// ============================================================================
//  Démo console du jour 8 : surface libre → CSV pour tools/plot_free_surface.py
// ============================================================================
//  Trois expériences, une par vérification du planning :
//    1. cuve au repos (10 s)      → rest.csv      : vitesses parasites, volume
//    2. profil de pression        → pressure.csv  : p contre ρ g (y_s − y)
//    3. cuve agitée (10 s), 32² et 64², avec et sans comblement des bulles
//                                 → stirred_n32_bulles8.csv, stirred_n32_bulles0.csv, …
//
//  Utilisation : free_surface_demo [dossier]   (par défaut : out_free_surface)
// ============================================================================
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include "diagnostics/FreeSurfaceDiagnostics.h"
#include "solvers/FlipSolver.h"

namespace fs = std::filesystem;

namespace {

constexpr real kG = real(9.81);

// Tourbillon à divergence nulle dans la moitié basse de la cuve (tangent aux parois et à la surface).
void stirTank(FlipSolver& s, real A) {
    ParticleSet2D& p = s.particles();
    for (std::size_t k = 0; k < p.size(); ++k) {
        const Vec2 x = p.x[k];
        p.v[k] = {A * std::sin(kPi * x.x) * std::cos(2 * kPi * x.y),
                  -A * real(0.5) * std::cos(kPi * x.x) * std::sin(2 * kPi * x.y)};
    }
}

// Fait tourner s pendant T secondes en écrivant une ligne CSV par pas.
FreeSurfaceStats run(FlipSolver& s, real T, const fs::path& csvPath) {
    std::ofstream csv(csvPath);
    writeFreeSurfaceCsvHeader(csv);
    FreeSurfaceStats st;
    real t = 0;
    while (t < T) {
        const real dt = s.chosenDt();
        s.step(dt);
        t += dt;
        st = measureFreeSurface(s.grid(), s.labels(), s.particles(), s.params().density, kG);
        writeFreeSurfaceCsvRow(csv, t, st);
    }
    return st;
}

} // namespace

int main(int argc, char** argv) {
    const fs::path dir = argc > 1 ? argv[1] : "out_free_surface";
    fs::create_directories(dir);

    // 1 et 2. Cuve au repos, 64 × 64, moitié pleine, ensemencement aléatoire.
    {
        FlipParams prm;
        prm.nx = prm.ny = 64;
        prm.tankHeight = 0.5;
        prm.jitter = 1;
        prm.density = 1000;                      // de l'eau : p en pascals
        prm.projectionTol = 1e-8;
        FlipSolver s(prm, ParticleScene::TankAtRest);
        const FreeSurfaceStats st = run(s, 10, dir / "rest.csv");
        std::printf("Cuve au repos 64x64, 10 s : max|u_p| = %.2e m/s, ecart a rho g h = %.2e, %d cellules fluides\n",
                    st.maxSpeed, st.hydroErr, st.fluidCells);
        std::ofstream prof(dir / "pressure.csv");
        writePressureProfileCsv(prof, s.grid(), s.labels(), prm.nx / 2, prm.density, kG);
    }

    // 3. Cuve agitée : effet du comblement des bulles sur le volume.
    for (int n : {32, 64})
        for (int bubbles : {8, 0}) {
            FlipParams prm;
            prm.nx = prm.ny = n;
            prm.tankHeight = 0.5;
            prm.autoDt = true;
            prm.maxBubbleCells = bubbles;
            FlipSolver s(prm, ParticleScene::TankAtRest);
            stirTank(s, real(0.5));
            s.step(s.chosenDt());
            const int a0 = s.lastStats().fluidCells;
            const std::string name = "stirred_n" + std::to_string(n) + "_bulles" + std::to_string(bubbles) + ".csv";
            const FreeSurfaceStats st = run(s, 10, dir / name);
            std::printf("Cuve agitee %dx%d, bulles comblees <= %d : volume final / initial = %.3f\n", n, n, bubbles,
                        double(st.fluidCells) / a0);
        }
    std::printf("CSV ecrits dans %s\n", dir.string().c_str());
    return 0;
}
