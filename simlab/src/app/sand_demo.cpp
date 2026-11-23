// Démo sans OpenGL (jour 13) : effondrement d'une colonne, eau contre sable (Zhu & Bridson §3).
//
// Colonne posée au centre du sol, lâchée à t = 0, schéma PIC pour tous les cas
// (Zhu & Bridson utilisent PIC pour le sable), même grille, même pas, même graine.
// Colonne par défaut 0,3 × 0,3 : moins haute que celle des tests (0,2 × 0,4), elle finit
// plus près d'un tas « au repos ». Avec la colonne haute et φ grand, il reste un pied de
// colonne tronqué et l'angle mesuré ne dit plus grand-chose.
//
//   sweep.csv     balayage de l'angle de frottement φ (15° à 50°), 3 graines par φ : forme du
//                 dépôt à t_fin (largeur, hauteur, angles des deux flancs), énergie cinétique finale
//   series.csv    évolution au cours du temps pour quelques cas : eau, sable φ = 20°, 30°, 45°,
//                 sable φ = 30° sans groupes rigides (il s'aplatit : rôle de l'étape 5) et
//                 sable φ = 30° sans deuxième projection (il se tasse : particules par cellule)
//   profiles.csv  profil de hauteur final h(x) des mêmes cas, plus la colonne initiale
//
// Usage : sand_demo [dossier_de_sortie] [n] [durée_s] [largeur] [hauteur]
//         (par défaut : sorties_sable 96 2 0.3 0.3 ; largeur et hauteur en fractions de Lx, Ly)
//         n = 96 : le tas fait une quinzaine de cellules de haut, assez pour mesurer ses
//         flancs. Compter 1 à 2 minutes ; n = 64 va deux fois plus vite, angles plus bruités.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "diagnostics/ParticleDiagnostics.h"
#include "diagnostics/SandDiagnostics.h"
#include "solvers/FlipSolver.h"

namespace fs = std::filesystem;

namespace {

struct Case {
    std::string name;
    Material material;
    real phiDeg;
    bool rigidGroups;
    bool reproject = true;   // deuxième projection après le frottement
};

real gColumnWidth = real(0.3), gColumnHeight = real(0.3);   // lus sur la ligne de commande

FlipParams baseParams(int n) {
    FlipParams prm;
    prm.nx = prm.ny = n;
    prm.scheme = TransferScheme::PIC;
    prm.dt = real(1) / 120;           // fixe, identique pour tous les cas
    prm.columnWidth = gColumnWidth;
    prm.columnHeight = gColumnHeight;
    return prm;
}

FlipParams caseParams(int n, const Case& c) {
    FlipParams prm = baseParams(n);
    prm.material = c.material;
    prm.frictionAngleDeg = c.phiDeg;
    prm.rigidGroups = c.rigidGroups;
    prm.reprojectAfterFriction = c.reproject;
    return prm;
}

int occupiedCells(const FlipSolver& s) {   // cellules contenant au moins une particule
    Field2D c = s.grid().makeCellField();
    countParticlesPerCell(s.particles(), c);
    int n = 0;
    for (real v : c.data()) n += v > 0;
    return n;
}

PileShape shapeOf(const FlipSolver& s) {
    return measurePile(s.particles(), s.grid().nx(), s.grid().h());
}

void writeProfile(std::ofstream& csv, const std::string& name, const FlipSolver& s) {
    const real h = s.grid().h();
    const std::vector<real> top = heightProfile(s.particles(), s.grid().nx(), h);
    for (std::size_t i = 0; i < top.size(); ++i)
        csv << name << ',' << (real(i) + real(0.5)) * h << ',' << top[i] << '\n';
}

} // namespace

int main(int argc, char** argv) {
    const fs::path out = argc > 1 ? argv[1] : "sorties_sable";
    const int n = argc > 2 ? std::max(16, std::atoi(argv[2])) : 96;
    const real tEnd = argc > 3 ? std::atof(argv[3]) : 2.0;
    if (argc > 4) gColumnWidth = real(std::atof(argv[4]));
    if (argc > 5) gColumnHeight = real(std::atof(argv[5]));
    std::printf("grille %d x %d, colonne %.2f x %.2f, t_fin = %.2f s\n", n, n, gColumnWidth, gColumnHeight, tEnd);
    fs::create_directories(out);
    const int steps = int(std::lround(tEnd / baseParams(n).dt));

    // ------------------------------------------------------------------------
    // 1) Séries temporelles et profils finaux.
    // ------------------------------------------------------------------------
    {
        const std::vector<Case> cases = {
            { "eau",                     Material::Water, 30, true },
            { "sable 20",                Material::Sand,  20, true },
            { "sable 30",                Material::Sand,  30, true },
            { "sable 45",                Material::Sand,  45, true },
            { "sable 30 sans groupes",   Material::Sand,  30, false },
            { "sable 30 sans 2e proj.",  Material::Sand,  30, true, false },
        };
        std::ofstream series(out / "series.csv");
        series << "cas,t,largeur,hauteur,angle,Ec,cellules_rigides,groupes,particules_par_cellule\n";
        std::ofstream profiles(out / "profiles.csv");
        profiles << "cas,x,h\n";

        bool initialWritten = false;
        for (const Case& c : cases) {
            const FlipParams prm = caseParams(n, c);
            FlipSolver s(prm, ParticleScene::ColumnCollapse);
            if (!initialWritten) { writeProfile(profiles, "initial", s); initialWritten = true; }
            auto write = [&] {
                const PileShape r = shapeOf(s);
                series << c.name << ',' << s.time() << ',' << r.spread << ',' << r.height << ','
                       << r.angleDeg << ',' << particleEnergy(s.particles()) << ','
                       << s.lastStats().rigidCells << ',' << s.lastStats().rigidGroups << ','
                       << double(s.particles().size()) / std::max(1, occupiedCells(s)) << '\n';
            };
            write();
            for (int k = 1; k <= steps; ++k) {
                s.step(prm.dt);
                if (k % 3 == 0 || k == steps) write();   // 40 points par seconde
            }
            writeProfile(profiles, c.name, s);
            const PileShape r = shapeOf(s);
            std::printf("%-22s  largeur %.3f   hauteur %.3f   angle %5.1f deg   Ec %.2e\n",
                        c.name.c_str(), r.spread, r.height, r.angleDeg, particleEnergy(s.particles()));
        }
    }

    // ------------------------------------------------------------------------
    // 2) Balayage de φ : forme du dépôt à t_fin, moyenne sur plusieurs graines
    //    (le dépôt dépend un peu du désordre initial des particules).
    // ------------------------------------------------------------------------
    {
        const unsigned seeds[] = { 1, 2, 3 };
        std::ofstream sweep(out / "sweep.csv");
        sweep << "phi,graine,largeur,hauteur,angle_gauche,angle_droit,angle,flancs,base,angle_equivalent,Ec_finale,"
                 "cellules_rigides,groupes,ms_par_pas\n";
        for (int phi = 15; phi <= 50; phi += 5) {
            real angleSum = 0, eqSum = 0, spreadSum = 0, msSum = 0;
            int measured = 0;
            for (unsigned seed : seeds) {
                const Case c{ "sable", Material::Sand, real(phi), true };
                FlipParams prm = caseParams(n, c);
                prm.seed = seed;
                FlipSolver s(prm, ParticleScene::ColumnCollapse);
                real ms = 0;
                for (int k = 1; k <= steps; ++k) {
                    s.step(prm.dt);
                    ms += s.lastStats().msTotal();
                }
                const PileShape r = shapeOf(s);
                sweep << phi << ',' << seed << ',' << r.spread << ',' << r.height << ',' << r.angleLeftDeg
                      << ',' << r.angleRightDeg << ',' << r.angleDeg << ',' << r.flanks << ','
                      << r.baseWidth << ',' << r.angleEquivalentDeg << ','
                      << particleEnergy(s.particles()) << ',' << s.lastStats().rigidCells << ','
                      << s.lastStats().rigidGroups << ',' << ms / steps << '\n';
                if (r.flanks > 0) { angleSum += r.angleDeg; ++measured; }
                spreadSum += r.spread;
                eqSum += r.angleEquivalentDeg;
                msSum += ms / steps;
            }
            const int ns = int(sizeof(seeds) / sizeof(seeds[0]));
            std::printf("phi = %2d deg : angle des flancs %5.1f deg   angle equivalent %5.1f deg   "
                        "largeur %.3f   %.2f ms/pas\n",
                        phi, measured ? angleSum / measured : 0.0, eqSum / ns, spreadSum / ns, msSum / ns);
        }
    }
    std::printf("CSV ecrits dans %s\n", out.string().c_str());
    return 0;
}
