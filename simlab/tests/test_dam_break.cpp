// Tests du jour 9 : rupture de barrage. Demande la surface libre du jour 8
// (gravité, étiquettes Fluid/Air/Solid, p = 0 dans l'air).
#include "doctest.h"
#include <algorithm>
#include <cmath>
#include <utility>
#include "diagnostics/ParticleDiagnostics.h"
#include "solvers/FlipSolver.h"

namespace {
// Cellules contenant au moins une particule (≈ cellules Fluid, sans comblement des bulles).
real occupiedCells(const FlipSolver& s) {
    Field2D c = s.grid().makeCellField();
    countParticlesPerCell(s.particles(), c);
    real n = 0;
    for (real v : c.data()) n += v > 0 ? 1 : 0;
    return n;
}

struct DamResult { real eTot0, eTot, area0, area, front; };

DamResult runDam(TransferScheme s, real alpha, real tEnd) {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    prm.scheme = s;
    prm.flipRatio = alpha;
    prm.dt = real(1) / 120;
    FlipSolver solver(prm, ParticleScene::DamBreak);
    auto eTot = [&] { return particleEnergy(solver.particles()) + solver.potentialEnergy(); };
    DamResult r{};
    r.eTot0 = eTot();
    r.area0 = occupiedCells(solver);
    const int n = int(std::lround(tEnd / prm.dt));
    for (int k = 0; k < n; ++k) solver.step(prm.dt);
    r.eTot = eTot();
    r.area = occupiedCells(solver);
    for (const Vec2& x : solver.particles().x) r.front = std::max(r.front, x.x);
    return r;
}
} // namespace

TEST_CASE("Rupture de barrage : l'eau avance, le volume tient, l'énergie totale ne croît pas") {
    for (const auto& [s, alpha] : {std::pair{TransferScheme::PIC, real(0)},
                                   std::pair{TransferScheme::FLIP, real(0.95)}}) {
        const DamResult r = runDam(s, alpha, 0.3);
        CHECK(r.front > 0.5);                              // le front a dépassé 2 fois la largeur initiale
        CHECK(std::abs(r.area / r.area0 - 1) < 0.25);
        CHECK(r.eTot <= r.eTot0 * 1.01);
    }
}

TEST_CASE("Rupture de barrage : PIC perd plus d'énergie mécanique que FLIP") {
    const DamResult pic = runDam(TransferScheme::PIC, 0, 1.0);
    const DamResult flip = runDam(TransferScheme::FLIP, 0.99, 1.0);
    MESSAGE("E(1 s)/E0 : PIC " << pic.eTot / pic.eTot0 << ", FLIP 0.99 " << flip.eTot / flip.eTot0);
    CHECK(pic.eTot < flip.eTot);
}
