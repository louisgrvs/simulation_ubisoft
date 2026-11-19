// Démo sans OpenGL (jours 9 et 10) : PIC contre FLIP contre APIC, énergie et bruit.
//
// Cinq schémas, avec exactement la même grille, le même pas de temps et la
// même graine : PIC, puis FLIP mélangé avec alpha = 0,95, 0,99 et 1 (FLIP pur),
// puis APIC (J10).
//
//   dam.csv      rupture de barrage (surface libre du jour 8, gravité) :
//                E_c, E_p, E_c + E_p, bruit, position du front, cellules occupées
//   vortex.csv   tourbillon en boîte pleine, sans gravité : E(t)/E0 et bruit.
//                L'énergie devrait rester constante : tout ce qui manque est
//                de la dissipation numérique (à comparer à vortex_pic.csv du jour 7)
//   summary.csv  valeurs finales, une ligne par expérience et par schéma
//
// Usage : flip_demo [dossier_de_sortie] [n] [durée_barrage_s]
//         (par défaut : sorties_flip 64 3)
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "diagnostics/ParticleDiagnostics.h"
#include "solvers/FlipSolver.h"

namespace fs = std::filesystem;

namespace {

struct Scheme {
    std::string name;
    TransferScheme scheme;
    real alpha;
};

const std::vector<Scheme> kSchemes = {
    {"PIC",       TransferScheme::PIC,  0},
    {"FLIP 0.95", TransferScheme::FLIP, real(0.95)},
    {"FLIP 0.99", TransferScheme::FLIP, real(0.99)},
    {"FLIP 1",    TransferScheme::FLIP, 1},
    {"APIC",      TransferScheme::APIC, 0},   // J10 : alpha ignoré
};

int occupiedCells(const FlipSolver& s) {
    Field2D c = s.grid().makeCellField();
    countParticlesPerCell(s.particles(), c);
    int n = 0;
    for (real v : c.data()) n += v > 0;
    return n;
}

real front(const FlipSolver& s) {
    real f = 0;
    for (const Vec2& x : s.particles().x) f = std::max(f, x.x);
    return f;
}

} // namespace

int main(int argc, char** argv) {
    const fs::path out = argc > 1 ? argv[1] : "sorties_flip";
    const int n = argc > 2 ? std::max(8, std::atoi(argv[2])) : 64;
    const real tDam = argc > 3 ? std::atof(argv[3]) : 3.0;
    fs::create_directories(out);

    std::ofstream summary(out / "summary.csv");
    summary << "experience,schema,E_final/E0,bruit_moyen,bruit_relatif_moyen,ms_par_pas\n";

    // ------------------------------------------------------------------------
    // 1) Rupture de barrage : colonne 0,25 × 0,5 contre la paroi gauche.
    // ------------------------------------------------------------------------
    {
        FlipParams prm;
        prm.nx = prm.ny = n;
        prm.dt = real(1) / 120;          // fixe, identique pour les 4 schémas
        const int steps = int(std::lround(tDam / prm.dt));
        std::ofstream csv(out / "dam.csv");
        csv << "schema,t,Ec,Ep,Etot,Etot/E0,bruit,bruit_relatif,front,cellules\n";

        for (const Scheme& sc : kSchemes) {
            FlipParams q = prm;
            q.scheme = sc.scheme;
            q.flipRatio = sc.alpha;
            FlipSolver s(q, ParticleScene::DamBreak);
            if (norm(s.activeGravity()) == 0) {
                std::printf("rupture de barrage ignoree : pas de gravite (surface libre du jour 8 absente)\n");
                break;
            }
            const real E0 = particleEnergy(s.particles()) + s.potentialEnergy();
            real noiseSum = 0, relSum = 0, msSum = 0;
            auto write = [&] {
                const real Ec = particleEnergy(s.particles()), Ep = s.potentialEnergy();
                csv << sc.name << ',' << s.time() << ',' << Ec << ',' << Ep << ',' << Ec + Ep << ','
                    << (Ec + Ep) / E0 << ',' << s.noise().meanVariance << ',' << s.noise().relative << ','
                    << front(s) << ',' << occupiedCells(s) << '\n';
            };
            write();
            for (int k = 1; k <= steps; ++k) {
                s.step(prm.dt);
                noiseSum += s.noise().meanVariance;
                relSum += s.noise().relative;
                msSum += s.lastStats().msTotal();
                write();
            }
            const real Ef = particleEnergy(s.particles()) + s.potentialEnergy();
            summary << "barrage," << sc.name << ',' << Ef / E0 << ',' << noiseSum / steps << ','
                    << relSum / steps << ',' << msSum / steps << '\n';
            std::printf("barrage  %-9s  (Ec+Ep)/E0 = %.4f   bruit moyen = %.3e m2/s2   %.2f ms/pas\n",
                        sc.name.c_str(), Ef / E0, noiseSum / steps, msSum / steps);
        }
    }

    // ------------------------------------------------------------------------
    // 2) Tourbillon cellulaire en boîte pleine (pas de gravité, pas de surface libre).
    // ------------------------------------------------------------------------
    {
        FlipParams prm;
        prm.nx = prm.ny = n;
        prm.dt = real(1) / 60;
        prm.projectionTol = 1e-8;
        const real T = 5;
        const int steps = int(std::lround(T / prm.dt));
        std::ofstream csv(out / "vortex.csv");
        csv << "schema,t,E/E0,bruit,bruit_relatif\n";

        for (const Scheme& sc : kSchemes) {
            FlipParams q = prm;
            q.scheme = sc.scheme;
            q.flipRatio = sc.alpha;
            FlipSolver s(q, ParticleScene::BoxVortex);
            const real E0 = particleEnergy(s.particles());
            real noiseSum = 0, relSum = 0, msSum = 0;
            csv << sc.name << ",0,1," << s.noise().meanVariance << ',' << s.noise().relative << '\n';
            for (int k = 1; k <= steps; ++k) {
                s.step(prm.dt);
                noiseSum += s.noise().meanVariance;
                relSum += s.noise().relative;
                msSum += s.lastStats().msTotal();
                csv << sc.name << ',' << s.time() << ',' << particleEnergy(s.particles()) / E0 << ','
                    << s.noise().meanVariance << ',' << s.noise().relative << '\n';
            }
            const real Ef = particleEnergy(s.particles()) / E0;
            summary << "tourbillon," << sc.name << ',' << Ef << ',' << noiseSum / steps << ','
                    << relSum / steps << ',' << msSum / steps << '\n';
            std::printf("vortex   %-9s  E/E0 = %.4f   bruit relatif moyen = %.3e   %.2f ms/pas\n",
                        sc.name.c_str(), Ef, relSum / steps, msSum / steps);
        }
    }
    std::printf("CSV ecrits dans %s\n", out.string().c_str());
    return 0;
}
