// ============================================================================
//  simlab_particles : particules PIC interactives (jour 7)
// ============================================================================
//  Même organisation que l'application de fumée (main.cpp, jour 5) :
//    1. événements (clavier, souris -> pinceau en coordonnées physiques) ;
//    2. un pas de simulation si on n'est pas en pause (ou si « Un pas ») ;
//    3. rendu : champ de grille en fond (facultatif), particules par-dessus ;
//    4. panneau ImGui.
//
//  Scènes :
//    Boite pleine - tourbillon : le tourbillon cellulaire du jour 5, porté par
//                                des particules ; comparer sa décroissance à Stam ;
//    Boite pleine - repos      : à remuer avec la souris (clic gauche) ;
//    Disque en rotation        : transferts PIC SEULS, particules immobiles.
//                                Avec l'échelle de couleur fixe, on voit la
//                                vitesse s'éteindre, d'abord au bord du disque.
//  Clavier : Espace = pause, N ou flèche droite = un pas, R = réinitialiser.
//
//  Les deux applications seront réunies quand FlipSolver aura son mode FLIP
//  (jour 9) ; d'ici là, chacune garde son panneau.
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
#include <vector>

#include "diagnostics/ParticleDiagnostics.h"
#include "render/FieldRenderer.h"
#include "render/ParticleRenderer.h"
#include "solvers/FlipSolver.h"

#ifndef SIMLAB_SHADER_DIR
#define SIMLAB_SHADER_DIR "shaders/"
#endif

namespace {

void glfwErrorCallback(int code, const char* message) {
    std::fprintf(stderr, "GLFW %d : %s\n", code, message);
}

void APIENTRY glDebugCallback(GLenum, GLenum type, GLuint, GLenum severity, GLsizei,
                              const GLchar* message, const void*) {
    if (severity == GL_DEBUG_SEVERITY_NOTIFICATION) return;
    std::fprintf(stderr, "[GL %s] %s\n", type == GL_DEBUG_TYPE_ERROR ? "ERREUR" : "info", message);
}

bool sliderReal(const char* label, real& value, float vmin, float vmax, const char* fmt,
                ImGuiSliderFlags flags = 0) {
    float f = float(value);
    if (ImGui::SliderFloat(label, &f, vmin, vmax, fmt, flags)) {
        value = real(f);
        return true;
    }
    return false;
}

// Plus grand rectangle aux proportions du domaine, centré dans la zone
// (même calcul que FieldRenderer::draw, pour que champ et particules coïncident).
FieldRenderer::Rect fitRect(FieldRenderer::Rect area, Vec2 domain) {
    FieldRenderer::Rect r;
    const float aspect = float(domain.x / domain.y);
    r.w = area.w;
    r.h = int(float(area.w) / aspect);
    if (r.h > area.h) {
        r.h = area.h;
        r.w = int(float(area.h) * aspect);
    }
    r.x = area.x + (area.w - r.w) / 2;
    r.y = area.y + (area.h - r.h) / 2;
    return r;
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

    GLFWwindow* window = glfwCreateWindow(1280, 720, "simlab - particules PIC (jour 7)", nullptr, nullptr);
    if (!window) {
        std::fprintf(stderr, "Echec de la creation de la fenetre (OpenGL 4.3 indisponible ?)\n");
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

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
    io.IniFilename = nullptr;
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 430");

    // --- 3. Simulation et rendu ---------------------------------------------------
    int exitCode = 0;
    try {
        FlipSolver solver;
        FieldRenderer fieldRenderer(SIMLAB_SHADER_DIR);
        ParticleRenderer particleRenderer(SIMLAB_SHADER_DIR);

        bool running = true, stepOnce = false;
        int background = 0;                          // 0 : aucun, sinon index + 1 dans fields()
        bool colorBySpeed = true, fixedScale = true;
        float pointScale = 1.0f;
        float speedScale = 1.0f;                     // haut de l'échelle de couleur (m/s)
        const int gridSizes[] = {32, 48, 64, 96, 128};
        int gridChoice = 2;
        int ppcChoice = 1;                           // 0 : 1, 1 : 4, 2 : 9 particules par cellule
        int seed = int(solver.params().seed);
        std::vector<float> energyPlot, angularPlot, speeds;
        const size_t plotLength = 600;

        FieldRenderer::Rect domainRect;
        Vec2 prevMouse;
        bool prevMouseValid = false;
        double prevTime = glfwGetTime();

        auto resetAll = [&]() {
            solver.params().seed = unsigned(std::max(seed, 0));
            solver.reset();
            energyPlot.clear();
            angularPlot.clear();
            const real vmax = particleMaxSpeed(solver.particles());
            speedScale = vmax > 0 ? float(vmax) : 1.0f;
        };
        resetAll();

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

            // --- Souris -> pinceau (comme au jour 5) -------------------------------
            int fbW = 0, fbH = 0, winW = 1, winH = 1;
            glfwGetFramebufferSize(window, &fbW, &fbH);
            glfwGetWindowSize(window, &winW, &winH);
            const double sx = double(fbW) / std::max(winW, 1), sy = double(fbH) / std::max(winH, 1);
            double mx = 0, my = 0;
            glfwGetCursorPos(window, &mx, &my);
            const double px = mx * sx, py = (winH - my) * sy;

            FlipBrush brush;
            const Vec2 L = solver.grid().size();
            const bool inDomain = domainRect.w > 0 && px >= domainRect.x && px < domainRect.x + domainRect.w &&
                                  py >= domainRect.y && py < domainRect.y + domainRect.h;
            if (inDomain) {
                const Vec2 pos{(px - domainRect.x) / domainRect.w * L.x, (py - domainRect.y) / domainRect.h * L.y};
                const bool left = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
                if (left && !io.WantCaptureMouse) {
                    Vec2 vel = prevMouseValid ? (pos - prevMouse) / real(frameDt) : Vec2{};
                    const real speed = norm(vel), vMaxBrush = 5;
                    if (speed > vMaxBrush) vel = (vMaxBrush / speed) * vel;
                    brush.active = true;
                    brush.pos = pos;
                    brush.velocity = vel;
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
                energyPlot.push_back(float(particleEnergy(solver.particles())));
                angularPlot.push_back(float(particleAngularMomentum(solver.particles(), solver.domainCenter())));
                if (energyPlot.size() > plotLength) energyPlot.erase(energyPlot.begin());
                if (angularPlot.size() > plotLength) angularPlot.erase(angularPlot.begin());
            }

            // --- Panneau --------------------------------------------------------------
            const float panelW = 380;
            ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImGui::SetNextWindowSize(ImVec2(panelW, io.DisplaySize.y));
            ImGui::Begin("Particules PIC", nullptr,
                         ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse);
            ImGui::PushItemWidth(170);
            FlipParams& p = solver.params();
            ImGui::Text("%.0f img/s   |   %.2f ms par pas", io.Framerate, solver.lastStats().msTotal());

            ImGui::SeparatorText("Simulation");
            bool paused = !running;
            if (ImGui::Checkbox("Pause (Espace)", &paused)) running = !paused;
            ImGui::SameLine();
            if (ImGui::Button("Un pas (N)")) stepOnce = true;
            ImGui::SameLine();
            if (ImGui::Button("Reinit. (R)")) resetAll();

            int scene = int(solver.scene());
            if (ImGui::Combo("Scene", &scene, "Boite pleine - tourbillon\0Boite pleine - repos\0"
                                              "Disque en rotation (transferts)\0\0")) {
                solver.setScene(ParticleScene(scene));
                resetAll();
            }
            if (ImGui::Combo("Grille", &gridChoice, "32 x 32\0" "48 x 48\0" "64 x 64\0" "96 x 96\0" "128 x 128\0\0")) {
                p.nx = p.ny = gridSizes[gridChoice];
                resetAll();
            }
            if (ImGui::Combo("Particules/cellule", &ppcChoice, "1\0" "4\0" "9\0\0")) {
                p.particlesPerAxis = ppcChoice + 1;
                resetAll();
            }
            sliderReal("desordre (jitter)", p.jitter, 0.0f, 1.0f, "%.2f");
            ImGui::InputInt("graine", &seed);
            ImGui::TextDisabled("desordre et graine : appliques a la reinit.");
            if (solver.transfersOnly())
                ImGui::TextColored(ImVec4(1.0f, 0.75f, 0.4f, 1.0f), "Transferts seuls : particules immobiles");

            ImGui::SeparatorText("Pas de temps");
            sliderReal("dt (s)", p.dt, 0.001f, 0.2f, "%.4f", ImGuiSliderFlags_Logarithmic);
            ImGui::Checkbox("dt auto : min(dt, CFL cible)", &p.autoDt);
            if (p.autoDt) sliderReal("CFL cible", p.cflTarget, 0.1f, 10.0f, "%.2f", ImGuiSliderFlags_Logarithmic);
            sliderReal("CFL des sous-pas", p.substepCfl, 0.25f, 5.0f, "%.2f");
            ImGui::Text("dt %.4f s   CFL %.2f   sous-pas %d", solver.chosenDt(), solver.lastStats().cfl,
                        solver.lastStats().substeps);

            if (ImGui::CollapsingHeader("Souris")) {
                ImGui::TextDisabled("clic gauche : pousser le fluide");
                sliderReal("rayon pinceau", p.brushRadius, 0.01f, 0.2f, "%.3f m");
                sliderReal("gain de force", p.brushGain, 0.0f, 100.0f, "%.1f /s");
            }

            ImGui::SeparatorText("Affichage");
            const std::vector<FieldView> views = solver.fields();
            const char* bgName = background == 0 ? "aucun" : views[background - 1].name.c_str();
            if (ImGui::BeginCombo("Fond", bgName)) {
                if (ImGui::Selectable("aucun", background == 0)) background = 0;
                for (int k = 0; k < int(views.size()); ++k)
                    if (ImGui::Selectable(views[k].name.c_str(), background == k + 1)) background = k + 1;
                ImGui::EndCombo();
            }
            ImGui::Checkbox("Couleur = |u_p|", &colorBySpeed);
            ImGui::SameLine();
            ImGui::Checkbox("Echelle fixe", &fixedScale);
            if (!fixedScale) speedScale = std::max(float(particleMaxSpeed(solver.particles())), 1e-6f);
            ImGui::DragFloat("|u| max (m/s)", &speedScale, 0.01f, 1e-3f, 10.0f, "%.3f");
            ImGui::SliderFloat("taille des points", &pointScale, 0.2f, 3.0f, "%.2f");

            ImGui::SeparatorText("Diagnostics");
            for (const Diagnostic& d : solver.diagnostics())
                ImGui::Text("%-22s %.4g", d.name.c_str(), d.value);
            if (!energyPlot.empty()) {
                ImGui::PlotLines("##E", energyPlot.data(), int(energyPlot.size()), 0, "energie (particules)",
                                 0.0f, FLT_MAX, ImVec2(-1, 60));
                ImGui::PlotLines("##L", angularPlot.data(), int(angularPlot.size()), 0, "moment cinetique",
                                 FLT_MAX, FLT_MAX, ImVec2(-1, 60));
            }
            ImGui::PopItemWidth();
            ImGui::End();

            // --- Rendu ------------------------------------------------------------------
            glViewport(0, 0, fbW, fbH);
            glClearColor(0.06f, 0.06f, 0.07f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);

            const int panelPx = int(panelW * sx);
            const int margin = int(12 * sx);
            const FieldRenderer::Rect area{panelPx + margin, margin, fbW - panelPx - 2 * margin, fbH - 2 * margin};
            domainRect = fitRect(area, L);

            // Cadre du domaine, légèrement plus clair que le fond.
            glEnable(GL_SCISSOR_TEST);
            glScissor(domainRect.x, domainRect.y, domainRect.w, domainRect.h);
            glClearColor(0.10f, 0.10f, 0.12f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            glDisable(GL_SCISSOR_TEST);

            if (background > 0) {
                const FieldView& view = views[background - 1];
                const float lo = float(view.field->minValue()), hi = float(view.field->maxValue());
                float vmin = std::min(lo, 0.0f), vmax = std::max(hi, vmin + 1e-12f);
                if (view.signedField) {
                    vmax = std::max({std::abs(lo), std::abs(hi), 1e-12f});
                    vmin = -vmax;
                }
                fieldRenderer.upload(*view.field);
                fieldRenderer.draw(domainRect, vmin, vmax,
                                   view.signedField ? FieldRenderer::Palette::Diverging
                                                    : FieldRenderer::Palette::Sequential);
            }

            const ParticleSet2D& ps = solver.particles();
            speeds.resize(ps.size());
            for (std::size_t k = 0; k < ps.size(); ++k) speeds[k] = float(norm(ps.v[k]));
            particleRenderer.upload(ps, speeds);
            // Diamètre ≈ espacement des particules (h / k) à l'écran, réglable.
            const float spacingPx = float(domainRect.w) / float(p.nx * p.particlesPerAxis);
            const float size = std::max(1.5f, pointScale * 0.9f * spacingPx);
            particleRenderer.draw(domainRect, L, size, 0.0f, speedScale, !colorBySpeed);

            ImGui::Render();
            glViewport(0, 0, fbW, fbH);
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
