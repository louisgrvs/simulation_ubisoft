// ============================================================================
//  simlab : fumée 2D interactive (jour 5)
// ============================================================================
//  Boucle de l'application, dans l'ordre à chaque image :
//    1. événements (clavier, souris -> pinceau en coordonnées physiques) ;
//    2. un pas de simulation si on n'est pas en pause (ou si « Un pas ») ;
//    3. rendu du champ choisi ;
//    4. panneau ImGui (paramètres, diagnostics, courbes).
//
//  Souris (dans la zone du champ) :
//    clic gauche : pousse le fluide dans le sens du geste et injecte de la fumée ;
//    clic droit  : injecte seulement de la fumée.
//  Clavier : Espace = pause, N ou flèche droite = un pas, R = réinitialiser.
//
//  Le texte de l'interface est sans accents : la police par défaut d'ImGui
//  ne contient que l'ASCII.
// ============================================================================
#include <glad/glad.h>   // toujours AVANT GLFW
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <vector>

#include "render/FieldRenderer.h"
#include "solvers/StableFluidsSolver.h"

#ifndef SIMLAB_SHADER_DIR
#define SIMLAB_SHADER_DIR "shaders/"
#endif

namespace {

void glfwErrorCallback(int code, const char* message) {
    std::fprintf(stderr, "GLFW %d : %s\n", code, message);
}

// Sortie de débogage OpenGL (cœur depuis 4.3) : chaque erreur est signalée ici,
// au moment même où elle se produit.
void APIENTRY glDebugCallback(GLenum, GLenum type, GLuint, GLenum severity, GLsizei,
                              const GLchar* message, const void*) {
    if (severity == GL_DEBUG_SEVERITY_NOTIFICATION) return;
    std::fprintf(stderr, "[GL %s] %s\n", type == GL_DEBUG_TYPE_ERROR ? "ERREUR" : "info", message);
}

// ImGui travaille en float, la simulation en double.
bool sliderReal(const char* label, real& value, float vmin, float vmax, const char* fmt,
                ImGuiSliderFlags flags = 0) {
    float f = float(value);
    if (ImGui::SliderFloat(label, &f, vmin, vmax, fmt, flags)) {
        value = real(f);
        return true;
    }
    return false;
}

// Une ligne de mesures par pas : exportée en CSV pour les tracés Python.
struct Sample {
    double t, dt, energy, smoke, divBefore, divAfter;
    int cgPressure;
    double ms;
};

void exportCsv(const std::vector<Sample>& log, const char* path) {
    std::ofstream out(path);
    out << "t,dt,energie,masse_fumee,div_avant,div_apres,iter_cg_pression,ms\n";
    for (const Sample& s : log)
        out << s.t << ',' << s.dt << ',' << s.energy << ',' << s.smoke << ',' << s.divBefore << ','
            << s.divAfter << ',' << s.cgPressure << ',' << s.ms << '\n';
}

} // namespace

int main() {
    // --- 1. Fenêtre et contexte OpenGL 4.3 Core ------------------------------------
    glfwSetErrorCallback(glfwErrorCallback);
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GLFW_TRUE);

    GLFWwindow* window = glfwCreateWindow(1280, 720, "simlab - Stable Fluids (jour 5)", nullptr, nullptr);
    if (!window) {
        std::fprintf(stderr, "Echec de la creation de la fenetre (OpenGL 4.3 indisponible ?)\n");
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);   // synchronisation verticale : au plus une image par rafraîchissement

    if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress))) {
        std::fprintf(stderr, "Echec de l'initialisation de glad\n");
        return 1;
    }
    std::printf("OpenGL %s | %s\n", glGetString(GL_VERSION), glGetString(GL_RENDERER));
    if (GLAD_GL_VERSION_4_3) {
        glEnable(GL_DEBUG_OUTPUT);
        glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
        glDebugMessageCallback(glDebugCallback, nullptr);
    }

    // --- 2. ImGui ---------------------------------------------------------------
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGuiIO& io = ImGui::GetIO();
    io.IniFilename = nullptr;               // pas de imgui.ini dans le dossier courant
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 430");

    // --- 3. Simulation et rendu ---------------------------------------------------
    int exitCode = 0;
    try {
        StableFluidsSolver solver;
        FieldRenderer renderer(SIMLAB_SHADER_DIR);

        bool running = true, stepOnce = false;
        int fieldIndex = 0;
        bool smooth = true, autoRange = true;
        float vmin = 0, vmax = 1;
        const int gridSizes[] = {32, 48, 64, 96, 128};
        int gridChoice = 2;                         // 64 × 64 par défaut
        std::vector<Sample> log;
        std::vector<float> energyPlot, smokePlot;   // dernières valeurs, pour les courbes
        const size_t plotLength = 600;

        FieldRenderer::Rect fieldRect;              // où le champ a été dessiné à l'image précédente
        Vec2 prevMouse;
        bool prevMouseValid = false;
        double prevTime = glfwGetTime();

        auto resetAll = [&]() {
            solver.reset();
            log.clear();
            energyPlot.clear();
            smokePlot.clear();
        };

        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            const double now = glfwGetTime();
            const double frameDt = std::max(now - prevTime, 1e-4);
            prevTime = now;

            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();

            // --- Clavier ---------------------------------------------------------
            if (!io.WantCaptureKeyboard) {
                if (ImGui::IsKeyPressed(ImGuiKey_Escape)) glfwSetWindowShouldClose(window, true);
                if (ImGui::IsKeyPressed(ImGuiKey_Space)) running = !running;
                if (ImGui::IsKeyPressed(ImGuiKey_N) || ImGui::IsKeyPressed(ImGuiKey_RightArrow)) stepOnce = true;
                if (ImGui::IsKeyPressed(ImGuiKey_R)) resetAll();
            }

            // --- Souris -> pinceau ------------------------------------------------
            // GLFW donne la position en pixels « fenêtre », origine en HAUT à gauche ;
            // le champ est dessiné en pixels « framebuffer », origine en BAS à gauche.
            int fbW = 0, fbH = 0, winW = 1, winH = 1;
            glfwGetFramebufferSize(window, &fbW, &fbH);
            glfwGetWindowSize(window, &winW, &winH);
            const double sx = double(fbW) / std::max(winW, 1), sy = double(fbH) / std::max(winH, 1);
            double mx = 0, my = 0;
            glfwGetCursorPos(window, &mx, &my);
            const double px = mx * sx, py = (winH - my) * sy;

            Brush brush;
            const Vec2 L = solver.grid().size();
            const bool inField = fieldRect.w > 0 && px >= fieldRect.x && px < fieldRect.x + fieldRect.w &&
                                 py >= fieldRect.y && py < fieldRect.y + fieldRect.h;
            if (inField) {
                const Vec2 pos{(px - fieldRect.x) / fieldRect.w * L.x, (py - fieldRect.y) / fieldRect.h * L.y};
                const bool left  = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
                const bool right = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
                if ((left || right) && !io.WantCaptureMouse) {
                    Vec2 vel = prevMouseValid ? (pos - prevMouse) / real(frameDt) : Vec2{};
                    const real speed = norm(vel), vMaxBrush = 5;      // garde-fou (m/s)
                    if (speed > vMaxBrush) vel = (vMaxBrush / speed) * vel;
                    brush.active   = true;
                    brush.pos      = pos;
                    brush.velocity = vel;
                    brush.push     = left;
                    brush.inject   = true;
                }
                prevMouse = pos;
                prevMouseValid = true;
            } else {
                prevMouseValid = false;
            }
            solver.setBrush(brush);

            // --- Pas de simulation ---------------------------------------------------
            if (running || stepOnce) {
                solver.step(solver.chosenDt());
                stepOnce = false;
                const StepStats& st = solver.lastStats();
                log.push_back({solver.time(), st.dt, solver.kineticEnergy(), solver.smokeAmount(),
                               st.divBefore, st.divAfter, st.pressureIterations, st.msTotal()});
                energyPlot.push_back(float(solver.kineticEnergy()));
                smokePlot.push_back(float(solver.smokeAmount()));
                if (energyPlot.size() > plotLength) energyPlot.erase(energyPlot.begin());
                if (smokePlot.size() > plotLength) smokePlot.erase(smokePlot.begin());
            }

            // --- Panneau --------------------------------------------------------------
            const float panelW = 380;
            ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImGui::SetNextWindowSize(ImVec2(panelW, io.DisplaySize.y));
            ImGui::Begin("Stable Fluids", nullptr,
                         ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
            ImGui::PushItemWidth(170);            // place laissée aux libellés, à droite
            StableFluidsParams& p = solver.params();
            ImGui::Text("%.0f img/s   |   %.2f ms par pas", io.Framerate, solver.lastStats().msTotal());

            ImGui::SeparatorText("Simulation");
            bool paused = !running;
            if (ImGui::Checkbox("Pause (Espace)", &paused)) running = !paused;
            ImGui::SameLine();
            if (ImGui::Button("Un pas (N)")) stepOnce = true;
            ImGui::SameLine();
            if (ImGui::Button("Reinit. (R)")) resetAll();

            int scene = int(solver.scene());
            if (ImGui::Combo("Scene", &scene, "Panache\0Bulle\0Tourbillon (mesure)\0\0")) {
                solver.setScene(SmokeScene(scene));
                resetAll();
            }
            if (ImGui::Combo("Grille", &gridChoice, "32 x 32\0" "48 x 48\0" "64 x 64\0" "96 x 96\0" "128 x 128\0\0")) {
                p.nx = p.ny = gridSizes[gridChoice];
                resetAll();
            }

            ImGui::SeparatorText("Pas de temps");
            sliderReal("dt (s)", p.dt, 0.001f, 1.0f, "%.4f", ImGuiSliderFlags_Logarithmic);
            ImGui::Checkbox("dt auto : min(dt, CFL cible)", &p.autoDt);
            if (p.autoDt) sliderReal("CFL cible", p.cflTarget, 0.1f, 20.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
            ImGui::Text("dt utilise %.4f s   CFL du dernier pas %.2f", solver.chosenDt(), solver.lastStats().cfl);

            ImGui::SeparatorText("Physique");
            const ImGuiSliderFlags logf = ImGuiSliderFlags_Logarithmic;
            sliderReal("nu viscosite", p.viscosity, 0.0f, 0.05f, "%.2e m2/s", logf);
            sliderReal("kappa diffusion", p.diffusivity, 0.0f, 0.05f, "%.2e m2/s", logf);
            sliderReal("alpha dissipation", p.dissipation, 0.0f, 5.0f, "%.3f /s", logf);
            sliderReal("beta flottabilite", p.buoyancy, -10.0f, 20.0f, "%.2f m/s2");

            if (ImGui::CollapsingHeader("Souris et emetteur")) {
                ImGui::TextDisabled("clic gauche : pousser + fumee");
                ImGui::TextDisabled("clic droit  : fumee seule");
                sliderReal("rayon pinceau", p.brushRadius, 0.01f, 0.2f, "%.3f m");
                sliderReal("gain de force", p.brushGain, 0.0f, 100.0f, "%.1f /s");
                sliderReal("debit pinceau", p.brushSmokeRate, 0.0f, 30.0f, "%.1f /s");
                if (solver.scene() == SmokeScene::Plume) {
                    ImGui::Checkbox("emetteur actif", &p.emitterOn);
                    sliderReal("debit emetteur", p.emitterRate, 0.0f, 20.0f, "%.1f /s");
                }
            }

            ImGui::SeparatorText("Affichage");
            const std::vector<FieldView> views = solver.fields();
            fieldIndex = std::clamp(fieldIndex, 0, int(views.size()) - 1);
            if (ImGui::BeginCombo("Champ", views[fieldIndex].name.c_str())) {
                for (int k = 0; k < int(views.size()); ++k)
                    if (ImGui::Selectable(views[k].name.c_str(), k == fieldIndex)) fieldIndex = k;
                ImGui::EndCombo();
            }
            ImGui::Checkbox("Lissage", &smooth);
            ImGui::SameLine();
            ImGui::Checkbox("Echelle auto", &autoRange);
            const FieldView& view = views[fieldIndex];
            if (autoRange) {
                const float lo = float(view.field->minValue()), hi = float(view.field->maxValue());
                if (view.signedField) {
                    const float m = std::max({std::abs(lo), std::abs(hi), 1e-12f});
                    vmin = -m;
                    vmax = m;
                } else {
                    vmin = std::min(lo, 0.0f);
                    vmax = std::max(hi, vmin + 1e-12f);
                }
            } else {
                ImGui::DragFloat("min", &vmin, 0.01f);
                ImGui::DragFloat("max", &vmax, 0.01f);
            }
            ImGui::Text("echelle [%.3g, %.3g]", vmin, vmax);

            ImGui::SeparatorText("Diagnostics");
            for (const Diagnostic& d : solver.diagnostics())
                ImGui::Text("%-22s %.4g", d.name.c_str(), d.value);
            if (!energyPlot.empty()) {
                ImGui::PlotLines("##E", energyPlot.data(), int(energyPlot.size()), 0, "energie cinetique",
                                 0.0f, FLT_MAX, ImVec2(-1, 70));
                ImGui::PlotLines("##M", smokePlot.data(), int(smokePlot.size()), 0, "masse de fumee",
                                 0.0f, FLT_MAX, ImVec2(-1, 70));
            }
            if (ImGui::Button("Exporter CSV")) {
                exportCsv(log, "simlab_mesures.csv");
                std::printf("%zu pas exportes dans simlab_mesures.csv\n", log.size());
            }
            ImGui::SameLine();
            ImGui::TextDisabled("%zu pas enregistres", log.size());
            ImGui::PopItemWidth();
            ImGui::End();

            // --- Rendu ------------------------------------------------------------------
            glViewport(0, 0, fbW, fbH);
            glClearColor(0.06f, 0.06f, 0.07f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);

            const int panelPx = int(panelW * sx);
            const int margin = int(12 * sx);
            FieldRenderer::Rect area{panelPx + margin, margin, fbW - panelPx - 2 * margin, fbH - 2 * margin};
            renderer.setSmooth(smooth);
            renderer.upload(*view.field);
            fieldRect = renderer.draw(area, vmin, vmax,
                                      view.signedField ? FieldRenderer::Palette::Diverging
                                                       : FieldRenderer::Palette::Sequential);

            ImGui::Render();
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            glfwSwapBuffers(window);
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "%s\n", e.what());
        exitCode = 1;
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return exitCode;
}
