// Démo sans OpenGL (jour 5) : mesures d'énergie de la boucle Stable Fluids complète.
//
// Scène : tourbillon cellulaire, solution STATIONNAIRE d'Euler (voir SmokeScenes.h).
// Toute perte d'énergie avec ν = 0 est donc numérique ; on la convertit en
// viscosité numérique effective :
//      E(T) = E(0) exp(-2 ν_eff k² T)   =>   ν_eff = -ln(E(T)/E(0)) / (2 k² T)
//
// Trois séries, chacune écrite en CSV :
//   energy_nu.csv    ν ∈ {0, 5e-4, 1e-3, 2e-3}, n = 64, dt = 0,01 s
//   energy_dt.csv    ν = 0, dt de 0,005 à 0,32 s (CFL ≈ 0,3 à 20), n = 64
//   energy_n.csv     ν = 0, dt = 0,01 s, n ∈ {16, 32, 64, 128}
// plus summary.csv (une ligne par simulation) et un test de stabilité du panache.
//
// Usage : smoke_demo [dossier_de_sortie]
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include "scenes/SmokeScenes.h"
#include "solvers/StableFluidsSolver.h"

namespace fs = std::filesystem;

namespace {

struct RunResult {
    real energyRatio = 0;   // E(T)/E(0)
    real nuEff = 0;         // viscosité effective totale (physique + numérique)
    real cfl = 0;
    double msPerStep = 0;
    int pressureIterations = 0;   // moyenne par pas
};

// Fait tourner le tourbillon jusqu'à T et écrit (série, t, E/E0, exact) dans csv.
RunResult runVortex(int n, real nu, real dt, real T, std::ofstream& csv, const std::string& series) {
    StableFluidsParams p;
    p.nx = p.ny = n;
    p.viscosity = nu;
    p.buoyancy = 0;
    p.emitterOn = false;
    p.projectionTol = 1e-8;
    StableFluidsSolver s(p, SmokeScene::Vortex);

    const real k2 = smoke_scene::cellularVortexK2(s.grid());
    const real e0 = s.kineticEnergy();
    const int steps = int(std::lround(T / dt));
    RunResult r;
    long iters = 0;
    double ms = 0;

    csv << series << ',' << 0 << ',' << 1 << ',' << 1 << '\n';
    for (int k = 1; k <= steps; ++k) {
        s.step(dt);
        ms += s.lastStats().msTotal();
        iters += s.lastStats().pressureIterations;
        r.cfl = std::max(r.cfl, s.lastStats().cfl);
        const real t = k * dt;
        csv << series << ',' << t << ',' << s.kineticEnergy() / e0 << ','
            << std::exp(-2 * nu * k2 * t) << '\n';
    }
    r.energyRatio = s.kineticEnergy() / e0;
    r.nuEff = -std::log(r.energyRatio) / (2 * k2 * steps * dt);
    r.msPerStep = ms / steps;
    r.pressureIterations = int(iters / steps);
    return r;
}

} // namespace

int main(int argc, char** argv) {
    const fs::path dir = argc > 1 ? argv[1] : "out_smoke";
    fs::create_directories(dir);
    const real T = 3.0;

    std::ofstream summary(dir / "summary.csv");
    summary << "serie,n,nu,dt,cfl,E_ratio,nu_eff,nu_num,ms_per_step,cg_iter\n";
    std::printf("%-10s %5s %8s %7s %6s %8s %10s %10s %8s %6s\n",
                "serie", "n", "nu", "dt", "CFL", "E(T)/E0", "nu_eff", "nu_num", "ms/pas", "CG");

    auto report = [&](const char* serie, int n, real nu, real dt, const RunResult& r) {
        const real nuNum = r.nuEff - nu;
        std::printf("%-10s %5d %8.1e %7.3f %6.2f %8.4f %10.3e %10.3e %8.2f %6d\n",
                    serie, n, nu, dt, r.cfl, r.energyRatio, r.nuEff, nuNum, r.msPerStep, r.pressureIterations);
        summary << serie << ',' << n << ',' << nu << ',' << dt << ',' << r.cfl << ',' << r.energyRatio
                << ',' << r.nuEff << ',' << nuNum << ',' << r.msPerStep << ',' << r.pressureIterations << '\n';
    };

    // 1) Viscosité physique : la décroissance s'accélère avec ν.
    {
        std::ofstream csv(dir / "energy_nu.csv");
        csv << "serie,t,E_ratio,exact_viscous\n";
        for (real nu : {0.0, 5e-4, 1e-3, 2e-3}) {
            char tag[32];
            std::snprintf(tag, sizeof tag, "nu=%g", nu);
            report("viscosite", 64, nu, 0.01, runVortex(64, nu, 0.01, T, csv, tag));
        }
    }
    // 2) Pas de temps : stable partout, et MOINS dissipatif à grand pas
    //    (moins d'interpolations par seconde simulée, comme au jour 3).
    {
        std::ofstream csv(dir / "energy_dt.csv");
        csv << "serie,t,E_ratio,exact_viscous\n";
        for (real dt : {0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32}) {
            char tag[32];
            std::snprintf(tag, sizeof tag, "dt=%g", dt);
            report("pas", 64, 0.0, dt, runVortex(64, 0.0, dt, T, csv, tag));
        }
    }
    // 3) Résolution : la viscosité numérique diminue quand h diminue.
    {
        std::ofstream csv(dir / "energy_n.csv");
        csv << "serie,t,E_ratio,exact_viscous\n";
        for (int n : {16, 32, 64, 128}) {
            char tag[32];
            std::snprintf(tag, sizeof tag, "n=%d", n);
            report("resolution", n, 0.0, 0.01, runVortex(n, 0.0, 0.01, T, csv, tag));
        }
    }
    // 4) Panache forcé à très grand pas : valeurs finies et vitesse bornée ?
    {
        std::printf("\nPanache (flottabilite + emetteur), 64x64, 5 s simulees :\n");
        for (real dt : {1.0 / 60, 0.1, 0.5, 1.0}) {
            StableFluidsParams p;
            StableFluidsSolver s(p, SmokeScene::Plume);
            const int steps = int(std::ceil(5.0 / dt));
            real vmax = 0;
            bool finite = true;
            for (int k = 0; k < steps; ++k) {
                s.step(dt);
                vmax = std::max(vmax, s.grid().maxSpeed());
                finite = finite && std::isfinite(s.kineticEnergy()) && std::isfinite(s.smokeAmount());
            }
            std::printf("  dt = %5.3f s : %4d pas, max|u| = %6.3f m/s, CFL max = %5.1f, %s\n",
                        dt, steps, vmax, vmax * dt / s.grid().h(), finite ? "valeurs finies" : "NaN/inf !");
        }
    }
    std::printf("\nSorties dans %s\n", dir.string().c_str());
    return 0;
}
