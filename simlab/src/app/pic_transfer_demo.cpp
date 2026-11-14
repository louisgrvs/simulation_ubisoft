// Démo sans OpenGL (jour 7) : mesure de la dissipation des transferts PIC.
//
// Expérience principale (APIC §8, transferts seuls) : un disque de rayon 0,3
// centré dans [0,1]², en rotation rigide u = ω × (x - c). On répète P2G puis G2P
// SANS déplacer les particules : toute perte de moment cinétique L ou d'énergie
// E vient des transferts. Résultats écrits en CSV :
//
//   transfers_n.csv       L/L0 et E/E0 sur 100 allers-retours, n ∈ {32, 64, 128}
//   transfers_jitter.csv  idem, n = 64, jitter ∈ {0 ; 0,5 ; 1}
//   transfers_ppc.csv     idem, n = 64, 1, 4 ou 9 particules par cellule
//   radial.csv            après UN aller-retour : erreur |u_p' - u_p| / (ω R)
//                         en fonction de r/R (où se concentre la perte ?)
//   vortex_pic.csv        boîte pleine, tourbillon cellulaire, solveur PIC
//                         complet : E(t)/E0, pour comparer à smoke_demo (jour 5)
//
// Usage : pic_transfer_demo [dossier_de_sortie]
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include "diagnostics/ParticleDiagnostics.h"
#include "numerics/MacTransfer.h"
#include "scenes/Seeding.h"
#include "solvers/FlipSolver.h"

namespace fs = std::filesystem;

namespace {

const Vec2 kCenter{0.5, 0.5};
const real kRadius = 0.3, kOmega = 1;

ParticleSet2D rotatingDisk(int n, int perAxis, real jitter, unsigned seed = 1) {
    ParticleSet2D p;
    SeedOptions opt;
    opt.perAxis = perAxis;
    opt.jitter = jitter;
    opt.seed = seed;
    seedParticles(p, n, n, real(1) / n, region::disk(kCenter, kRadius), opt, [](Vec2 x) {
        const Vec2 r = x - kCenter;
        return Vec2{-kOmega * r.y, kOmega * r.x};
    });
    return p;
}

// Répète `trips` allers-retours ; écrit (série, k, L/L0, E/E0) ; renvoie L/L0 final.
real runTransfers(int n, int perAxis, real jitter, int trips, std::ofstream& csv, const std::string& series) {
    MacGrid2D g(n, n, real(1) / n);
    MacMasses m;
    ParticleSet2D p = rotatingDisk(n, perAxis, jitter);
    const real L0 = particleAngularMomentum(p, kCenter), E0 = particleEnergy(p);
    csv << series << ",0,1,1\n";
    real ratio = 1;
    for (int k = 1; k <= trips; ++k) {
        particlesToMac(p, g, m);
        macToParticlesPIC(g, p);
        ratio = particleAngularMomentum(p, kCenter) / L0;
        csv << series << ',' << k << ',' << ratio << ',' << particleEnergy(p) / E0 << '\n';
    }
    return ratio;
}

} // namespace

int main(int argc, char** argv) {
    const fs::path dir = argc > 1 ? argv[1] : "out_pic";
    fs::create_directories(dir);
    const int trips = 100;

    std::printf("Disque en rotation rigide (R = %.1f, omega = %.0f), transferts PIC seuls\n", kRadius, kOmega);
    std::printf("%-12s %5s %4s %7s %10s %14s %16s\n", "serie", "n", "ppc", "jitter", "particules",
                "L/L0 (1 a.-r.)", "L/L0 (100 a.-r.)");

    auto row = [&](const char* serie, int n, int k, real jitter, std::ofstream& csv, const std::string& tag) {
        {
            MacGrid2D g(n, n, real(1) / n);
            MacMasses m;
            ParticleSet2D p = rotatingDisk(n, k, jitter);
            const real L0 = particleAngularMomentum(p, kCenter);
            particlesToMac(p, g, m);
            macToParticlesPIC(g, p);
            const real one = particleAngularMomentum(p, kCenter) / L0;
            const real hundred = runTransfers(n, k, jitter, trips, csv, tag);
            std::printf("%-12s %5d %4d %7.2f %10zu %14.6f %16.6f\n", serie, n, k * k, jitter, p.size(), one, hundred);
        }
    };

    // 1) Résolution
    {
        std::ofstream csv(dir / "transfers_n.csv");
        csv << "serie,k,L_ratio,E_ratio\n";
        for (int n : {32, 64, 128}) row("resolution", n, 2, 1, csv, "n=" + std::to_string(n));
    }
    // 2) Désordre de l'ensemencement
    {
        std::ofstream csv(dir / "transfers_jitter.csv");
        csv << "serie,k,L_ratio,E_ratio\n";
        for (real j : {0.0, 0.5, 1.0}) {
            char tag[32];
            std::snprintf(tag, sizeof tag, "jitter=%g", j);
            row("jitter", 64, 2, j, csv, tag);
        }
    }
    // 3) Particules par cellule
    {
        std::ofstream csv(dir / "transfers_ppc.csv");
        csv << "serie,k,L_ratio,E_ratio\n";
        for (int k : {1, 2, 3}) row("ppc", 64, k, 1, csv, std::to_string(k * k) + " par cellule");
    }
    // 4) Où se perd la vitesse ? Erreur après un aller-retour selon le rayon,
    //    pour un ensemencement aléatoire (jitter 1) et un réseau régulier (jitter 0).
    {
        const int n = 64;
        std::ofstream csv(dir / "radial.csv");
        csv << "jitter,r_over_R,err_over_omegaR\n";
        std::printf("\nApres 1 aller-retour (n = 64), erreur moyenne |u' - u| / (omega R) :\n");
        for (real jitter : {1.0, 0.0}) {
            MacGrid2D g(n, n, real(1) / n);
            MacMasses m;
            ParticleSet2D p = rotatingDisk(n, 2, jitter);
            const ParticleSet2D p0 = p;
            particlesToMac(p, g, m);
            macToParticlesPIC(g, p);
            real edge = 0, inner = 0;
            int nEdge = 0, nInner = 0;
            for (std::size_t k = 0; k < p.size(); ++k) {
                const real r = norm(p0.x[k] - kCenter) / kRadius;
                const real e = norm(p.v[k] - p0.v[k]) / (kOmega * kRadius);
                csv << jitter << ',' << r << ',' << e << '\n';
                if (r > 1 - 2.0 / (n * kRadius)) { edge += e; ++nEdge; }    // à moins de 2h du bord
                else { inner += e; ++nInner; }
            }
            std::printf("  jitter = %.0f : %.2e dans les 2h du bord, %.2e a l'interieur\n", jitter,
                        edge / nEdge, inner / nInner);
        }
    }
    // 5) Solveur PIC complet sur le tourbillon cellulaire du jour 5.
    {
        std::ofstream csv(dir / "vortex_pic.csv");
        csv << "serie,t,E_ratio\n";
        const real T = 3, dt = 0.01, k2 = 2 * kPi * kPi;
        std::printf("\nTourbillon cellulaire en boite pleine, solveur PIC, dt = %.2f s, T = %.0f s :\n", dt, T);
        for (int n : {32, 64}) {
            FlipParams prm;
            prm.nx = prm.ny = n;
            prm.projectionTol = 1e-8;
            FlipSolver s(prm, ParticleScene::BoxVortex);
            const real E0 = particleEnergy(s.particles());
            const int steps = int(std::lround(T / dt));
            double ms = 0;
            csv << "n=" << n << ",0,1\n";
            for (int k = 1; k <= steps; ++k) {
                s.step(dt);
                ms += s.lastStats().msTotal();
                csv << "n=" << n << ',' << k * dt << ',' << particleEnergy(s.particles()) / E0 << '\n';
            }
            const real ratio = particleEnergy(s.particles()) / E0;
            std::printf("  n = %3d : %6zu particules, E(T)/E0 = %.4f, nu_num = %.2e m2/s, %.2f ms par pas\n", n,
                        s.particles().size(), ratio, -std::log(ratio) / (2 * k2 * T), ms / steps);
        }
    }
    std::printf("\nSorties dans %s\n", dir.string().c_str());
    return 0;
}
