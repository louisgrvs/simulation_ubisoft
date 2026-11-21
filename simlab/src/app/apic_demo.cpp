// ============================================================================
//  apic_demo : jour 11 (APIC, §8, fig. 4 et 13), PIC contre FLIP contre APIC
// ============================================================================
//  Mêmes schémas, même grille, même pas de temps et même graine pour tous :
//  PIC, FLIP mélangé avec alpha = 0,95, FLIP pur (alpha = 1) et APIC.
//
//  1) disk_transfers.csv : disque en rotation, TRANSFERTS SEULS (particules immobiles),
//     100 allers-retours P2G + G2P : L/L0 et E/E0.
//     Attendu : PIC perd à chaque aller-retour ; APIC exact (rotation = champ affine) ;
//     FLIP exact lui aussi, mais trivialement : sans projection ni force, la grille
//     ne change pas entre P2G et G2P, donc l'incrément ajouté aux particules est nul.
//
//  2) disk_free.csv : le même disque en simulation complète (RotatingDiskFree,
//     surface libre, sans gravité), sur plusieurs tours. Par pas :
//       L/L0      moment cinétique des particules autour du centre de masse
//       Lg/Lg0    celui de la grille juste après P2G (la mesure commune aux trois
//                 schémas : en APIC, une partie de L est portée par cx, cy)
//       E/E0, bruit relatif, R_g/R_g0 (forme), |P|, r_max/R (doigts), cellules occupées
//     Attention : sans tension de surface, ce disque est instable (Rayleigh-Taylor).
//     Il ne reste un disque que pendant environ un quart de tour ; au-delà, on mesure
//     sa désagrégation, plus la conservation par les schémas.
//
//  3) stability.csv : balayage en pas de temps (dt = 1/240 à 1/4 s), pour la rupture
//     de barrage (1 s) et pour le disque libre (un quart de tour), avec et sans sous-pas d'advection
//     (sans : un seul déplacement RK2 par pas, quel que soit le CFL). Un calcul est
//     déclaré « explosé » dès qu'apparaît une valeur non finie, une vitesse au-delà
//     de 10 fois l'échelle de la scène, ou une énergie au-delà de 1,05 E0.
//
//  4) summary.csv : une ligne par expérience et par schéma.
//
//  Usage : apic_demo [dossier_de_sortie] [n] [tours_du_disque]
//          (par défaut : sorties_apic, 64, 3)
// ============================================================================
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
    {"FLIP 1",    TransferScheme::FLIP, 1},
    {"APIC",      TransferScheme::APIC, 0},   // alpha ignoré
};

FlipParams paramsFor(const FlipParams& base, const Scheme& sc) {
    FlipParams q = base;
    q.scheme = sc.scheme;
    q.flipRatio = sc.alpha;
    return q;
}

int occupiedCells(const FlipSolver& s) {
    Field2D c = s.grid().makeCellField();
    countParticlesPerCell(s.particles(), c);
    int n = 0;
    for (real v : c.data()) n += v > 0;
    return n;
}

// Plus grande distance d'une particule au centre de masse : les doigts de l'instabilité.
real maxRadius(const FlipSolver& s) {
    const Vec2 G = s.centerOfMass();
    real r = 0;
    for (const Vec2& x : s.particles().x) r = std::max(r, norm(x - G));
    return r;
}

bool allFinite(const ParticleSet2D& p) {
    for (std::size_t k = 0; k < p.size(); ++k)
        if (!std::isfinite(p.x[k].x) || !std::isfinite(p.x[k].y) ||
            !std::isfinite(p.v[k].x) || !std::isfinite(p.v[k].y))
            return false;
    return true;
}

// Énergie mécanique : cinétique + potentielle (nulle sans gravité).
real totalEnergy(const FlipSolver& s) { return particleEnergy(s.particles()) + s.potentialEnergy(); }

} // namespace

int main(int argc, char** argv) {
    const fs::path out = argc > 1 ? argv[1] : "sorties_apic";
    const int n = argc > 2 ? std::max(8, std::atoi(argv[2])) : 64;
    const real turns = argc > 3 ? std::max(0.25, std::atof(argv[3])) : 3.0;
    fs::create_directories(out);

    std::ofstream summary(out / "summary.csv");
    summary << "experience,schema,L/L0,Lg/Lg0,E/E0,Rg/Rg0,bruit_relatif,ms_par_pas\n";

    FlipParams base;
    base.nx = base.ny = n;
    base.dt = real(1) / 60;
    base.autoDt = false;              // le pas est imposé partout : c'est lui qu'on compare
    const real w = base.diskAngularVelocity, R = base.diskRadius;
    const real period = 2 * kPi / w;  // un tour : pi s pour omega = 2

    // ------------------------------------------------------------------------
    // 1) Disque en rotation, transferts seuls : 100 allers-retours.
    // ------------------------------------------------------------------------
    {
        const int trips = 100;
        std::ofstream csv(out / "disk_transfers.csv");
        csv << "schema,k,L/L0,E/E0\n";
        for (const Scheme& sc : kSchemes) {
            FlipSolver s(paramsFor(base, sc), ParticleScene::RotatingDisk);
            const Vec2 G = s.centerOfMass();   // les particules ne bougent pas : G est fixe
            const real L0 = particleAngularMomentum(s.particles(), G);
            const real E0 = particleEnergy(s.particles());
            csv << sc.name << ",0,1,1\n";
            for (int k = 1; k <= trips; ++k) {
                s.step(base.dt);
                csv << sc.name << ',' << k << ',' << particleAngularMomentum(s.particles(), G) / L0 << ','
                    << particleEnergy(s.particles()) / E0 << '\n';
            }
            const real Lf = particleAngularMomentum(s.particles(), G) / L0;
            const real Ef = particleEnergy(s.particles()) / E0;
            summary << "disque_transferts_100," << sc.name << ',' << Lf << ",," << Ef << ",,,"
                    << s.lastStats().msTotal() << '\n';
            std::printf("disque, transferts seuls  %-9s  apres %d allers-retours : L/L0 = %.6f   E/E0 = %.6f\n",
                        sc.name.c_str(), trips, Lf, Ef);
        }
    }

    // ------------------------------------------------------------------------
    // 2) Disque en rotation, simulation complète (surface libre, sans gravité).
    // ------------------------------------------------------------------------
    {
        const int steps = int(std::lround(turns * period / base.dt));
        const int quarter = int(std::lround(0.25 * period / base.dt));   // le disque tient jusque-là
        std::ofstream csv(out / "disk_free.csv");
        csv << "schema,t,L/L0,Lg/Lg0,E/E0,bruit_relatif,Rg/Rg0,P,rmax/R,cellules\n";
        for (const Scheme& sc : kSchemes) {
            FlipSolver s(paramsFor(base, sc), ParticleScene::RotatingDiskFree);
            const real L0 = particleAngularMomentum(s.particles(), s.centerOfMass());
            const real E0 = particleEnergy(s.particles());
            const real Rg0 = s.initialRadiusOfGyration();
            real Lg0 = 0;   // connu au premier P2G
            real relSum = 0, msSum = 0;
            auto row = [&](real Lg) {
                const Vec2 G = s.centerOfMass();
                csv << sc.name << ',' << s.time() << ',' << particleAngularMomentum(s.particles(), G) / L0 << ','
                    << Lg << ',' << particleEnergy(s.particles()) / E0 << ',' << s.noise().relative << ','
                    << particleRadiusOfGyration(s.particles(), G) / Rg0 << ','
                    << norm(particleMomentum(s.particles())) << ',' << maxRadius(s) / R << ','
                    << occupiedCells(s) << '\n';
            };
            auto summarize = [&](const char* name, int k) {
                const Vec2 G = s.centerOfMass();
                summary << name << ',' << sc.name << ',' << particleAngularMomentum(s.particles(), G) / L0 << ','
                        << s.lastStats().angularMomentumGrid / Lg0 << ',' << particleEnergy(s.particles()) / E0 << ','
                        << particleRadiusOfGyration(s.particles(), G) / Rg0 << ',' << relSum / k << ','
                        << msSum / k << '\n';
            };
            row(1);
            for (int k = 1; k <= steps; ++k) {
                s.step(base.dt);
                if (k == 1) Lg0 = s.lastStats().angularMomentumGrid;
                relSum += s.noise().relative;
                msSum += s.lastStats().msTotal();
                row(s.lastStats().angularMomentumGrid / Lg0);
                if (k == quarter) {
                    summarize("disque_libre_quart_de_tour", k);
                    std::printf("disque libre, 1/4 de tour %-9s  Lg/Lg0 = %.4f   E/E0 = %.4f   R_g/R_g0 = %.4f\n",
                                sc.name.c_str(), s.lastStats().angularMomentumGrid / Lg0,
                                particleEnergy(s.particles()) / E0,
                                particleRadiusOfGyration(s.particles(), s.centerOfMass()) / Rg0);
                }
            }
            summarize("disque_libre_fin", steps);
            std::printf("disque libre, %.2g tours   %-9s  Lg/Lg0 = %.4f   E/E0 = %.4f   R_g/R_g0 = %.4f   %.2f ms/pas\n",
                        double(turns), sc.name.c_str(), s.lastStats().angularMomentumGrid / Lg0,
                        particleEnergy(s.particles()) / E0,
                        particleRadiusOfGyration(s.particles(), s.centerOfMass()) / Rg0, msSum / steps);
        }
    }

    // ------------------------------------------------------------------------
    // 3) Stabilité selon le pas de temps.
    // ------------------------------------------------------------------------
    {
        struct SceneCase {
            const char* name;
            ParticleScene scene;
            real duration;   // s
            real speed;      // échelle de vitesse de la scène (m/s)
        };
        const real damSpeed = std::sqrt(2 * norm(base.gravity) * base.damHeight);   // chute de la colonne
        const std::vector<SceneCase> scenes = {
            {"barrage",      ParticleScene::DamBreak,         1,              damSpeed},
            {"disque_libre", ParticleScene::RotatingDiskFree, real(0.25) * period, w * R},
        };
        const std::vector<int> perSecond = {240, 120, 60, 30, 15, 8, 4};   // dt = 1 / perSecond

        std::ofstream csv(out / "stability.csv");
        csv << "scene,schema,sous_pas,dt,pas,cfl_max,sous_pas_max,explose,t_explosion,E_final/E0,"
               "bruit_relatif_moyen,ms_par_pas\n";
        for (const SceneCase& sce : scenes)
            for (bool substeps : {true, false})
                for (const Scheme& sc : kSchemes)
                    for (int ps : perSecond) {
                        FlipParams q = paramsFor(base, sc);
                        q.dt = real(1) / ps;
                        // Sans sous-pas : un CFL de sous-pas immense donne un seul déplacement par pas.
                        q.substepCfl = substeps ? real(1) : real(1e9);
                        FlipSolver s(q, sce.scene);
                        const real E0 = totalEnergy(s);
                        const int steps = int(std::lround(sce.duration / q.dt));
                        real cflMax = 0, relSum = 0, msSum = 0, tBlow = -1;
                        int subMax = 0, done = 0;
                        for (int k = 1; k <= steps; ++k) {
                            s.step(q.dt);
                            ++done;
                            cflMax = std::max(cflMax, s.lastStats().cfl);
                            subMax = std::max(subMax, s.lastStats().substeps);
                            relSum += s.noise().relative;
                            msSum += s.lastStats().msTotal();
                            const bool blown = !allFinite(s.particles()) ||
                                               particleMaxSpeed(s.particles()) > 10 * sce.speed ||
                                               totalEnergy(s) > real(1.05) * E0;
                            if (blown) {
                                tBlow = s.time();
                                break;
                            }
                        }
                        const real Ef = totalEnergy(s) / E0;
                        csv << sce.name << ',' << sc.name << ',' << (substeps ? 1 : 0) << ',' << q.dt << ','
                            << done << ',' << cflMax << ',' << subMax << ',' << (tBlow >= 0 ? 1 : 0) << ','
                            << tBlow << ',' << Ef << ',' << relSum / done << ',' << msSum / done << '\n';
                        std::printf("stabilite %-12s %-9s %s dt = 1/%-3d CFL max = %6.2f  %s  E/E0 = %.4f\n",
                                    sce.name, sc.name.c_str(), substeps ? "avec sous-pas" : "sans sous-pas",
                                    ps, double(cflMax), tBlow >= 0 ? "EXPLOSE" : "stable ", double(Ef));
                    }
    }

    std::printf("CSV ecrits dans %s\n", out.string().c_str());
    return 0;
}
