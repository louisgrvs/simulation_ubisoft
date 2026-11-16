// Tests du jour 9 : transfert FLIP et mélange PIC/FLIP, bruit, rupture de barrage.
#include "doctest.h"
#include <algorithm>
#include <cmath>
#include <random>
#include "diagnostics/ParticleDiagnostics.h"
#include "numerics/MacTransfer.h"
#include "scenes/Seeding.h"
#include "solvers/FlipSolver.h"

namespace {

void randomVelocities(ParticleSet2D& p, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<real> U(-1, 1);
    for (Vec2& v : p.v) v = {U(rng), U(rng)};
}

// Particules partout (pochoirs tronqués aux parois compris), vitesses aléatoires.
ParticleSet2D fullBox(int n, unsigned seed = 3) {
    ParticleSet2D p;
    SeedOptions opt;
    opt.seed = seed;
    seedParticles(p, n, n, real(1) / n, region::everywhere(), opt);
    randomVelocities(p, seed + 1);
    return p;
}

real maxVelocityGap(const ParticleSet2D& a, const ParticleSet2D& b) {
    real e = 0;
    for (std::size_t k = 0; k < a.size(); ++k) e = std::max(e, norm(a.v[k] - b.v[k]));
    return e;
}

// Grille « après les forces » : la grille de P2G modifiée par un champ quelconque.
void perturbGrid(MacGrid2D& g, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<real> U(-0.5, 0.5);
    for (real& u : g.u.data()) u += U(rng);
    for (real& v : g.v.data()) v += U(rng);
}

} // namespace

// ---------------------------------------------------------------------------
//  Transfert FLIP seul
// ---------------------------------------------------------------------------

TEST_CASE("FLIP avec alpha = 0 : exactement PIC") {
    const int n = 16;
    MacGrid2D g(n, n, real(1) / n);
    MacMasses m;
    ParticleSet2D p = fullBox(n);
    particlesToMac(p, g, m);
    const Field2D uOld = g.u, vOld = g.v;
    perturbGrid(g, 11);

    ParticleSet2D pic = p, flip = p;
    macToParticlesPIC(g, pic);
    macToParticlesFLIP(g, uOld, vOld, flip, 0);
    CHECK(maxVelocityGap(pic, flip) < 1e-14);
}

TEST_CASE("FLIP pur sans variation de la grille : les vitesses des particules sont inchangées") {
    // C'est toute la différence avec PIC : PIC remplacerait u_p par la moyenne de grille.
    const int n = 16;
    MacGrid2D g(n, n, real(1) / n);
    MacMasses m;
    ParticleSet2D p = fullBox(n);
    particlesToMac(p, g, m);

    ParticleSet2D flip = p, pic = p;
    macToParticlesFLIP(g, g.u, g.v, flip, 1);
    macToParticlesPIC(g, pic);
    CHECK(maxVelocityGap(p, flip) == 0.0);
    CHECK(maxVelocityGap(p, pic) > 0.1);       // PIC, lui, lisse les vitesses aléatoires
}

TEST_CASE("FLIP pur : une variation uniforme de la grille s'ajoute telle quelle, parois comprises") {
    // Partition de l'unité des poids, y compris les pochoirs tronqués (renormalisés).
    const int n = 12;
    MacGrid2D g(n, n, real(1) / n);
    MacMasses m;
    ParticleSet2D p = fullBox(n);
    particlesToMac(p, g, m);
    const Field2D uOld = g.u, vOld = g.v;
    const Vec2 dv{0.3, -0.7};
    for (real& u : g.u.data()) u += dv.x;
    for (real& v : g.v.data()) v += dv.y;

    ParticleSet2D flip = p;
    macToParticlesFLIP(g, uOld, vOld, flip, 1);
    real err = 0;
    for (std::size_t k = 0; k < p.size(); ++k) err = std::max(err, norm(flip.v[k] - (p.v[k] + dv)));
    CHECK(err < 1e-13);
}

TEST_CASE("Mélange : linéaire en alpha entre PIC et FLIP") {
    const int n = 16;
    MacGrid2D g(n, n, real(1) / n);
    MacMasses m;
    ParticleSet2D p = fullBox(n);
    particlesToMac(p, g, m);
    const Field2D uOld = g.u, vOld = g.v;
    perturbGrid(g, 5);

    ParticleSet2D pic = p, flip = p, mix = p;
    macToParticlesFLIP(g, uOld, vOld, pic, 0);
    macToParticlesFLIP(g, uOld, vOld, flip, 1);
    macToParticlesFLIP(g, uOld, vOld, mix, 0.95);
    real err = 0;
    for (std::size_t k = 0; k < p.size(); ++k)
        err = std::max(err, norm(mix.v[k] - (0.95 * flip.v[k] + 0.05 * pic.v[k])));
    CHECK(err < 1e-14);
}

// ---------------------------------------------------------------------------
//  Bruit : variance des vitesses par cellule
// ---------------------------------------------------------------------------

TEST_CASE("Bruit : nul pour un champ uniforme, exact pour deux vitesses opposées") {
    const int n = 8;
    const real h = real(1) / n;
    Field2D var(n, n, h, 0.5, 0.5);

    ParticleSet2D p;
    SeedOptions opt;
    seedParticles(p, n, n, h, region::everywhere(), opt, [](Vec2) { return Vec2{1.5, -2}; });
    NoiseStats st = cellVelocityVariance(p, var);
    CHECK(st.meanVariance < 1e-28);
    CHECK(st.cells == n * n);

    // Dans chaque cellule, moitié des particules à +a, moitié à -a : sigma² = |a|².
    const Vec2 a{0.6, 0.8};   // |a| = 1
    for (std::size_t k = 0; k < p.size(); ++k) p.v[k] = (k % 2 == 0) ? a : -a;
    st = cellVelocityVariance(p, var);
    CHECK(st.meanVariance == doctest::Approx(1.0).epsilon(1e-12));
    CHECK(st.relative == doctest::Approx(1.0).epsilon(1e-12));   // <|u|²> = 1 aussi
    CHECK(var(3, 5) == doctest::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("Bruit : une cellule à une seule particule est ignorée") {
    const real h = 0.25;
    Field2D var(4, 4, h, 0.5, 0.5);
    ParticleSet2D p;
    p.add({0.1, 0.1}, {5, 0}, 1);
    const NoiseStats st = cellVelocityVariance(p, var);
    CHECK(st.cells == 0);
    CHECK(st.meanVariance == 0.0);
}

// ---------------------------------------------------------------------------
//  Solveur complet, sans surface libre
// ---------------------------------------------------------------------------

TEST_CASE("Disque en rotation, transferts seuls : FLIP pur ne dissipe rien, PIC si") {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    prm.scheme = TransferScheme::FLIP;
    prm.flipRatio = 1;
    FlipSolver flip(prm, ParticleScene::RotatingDisk);
    prm.scheme = TransferScheme::PIC;
    FlipSolver pic(prm, ParticleScene::RotatingDisk);
    const real E0 = particleEnergy(flip.particles());
    for (int k = 0; k < 50; ++k) { flip.step(0.01); pic.step(0.01); }
    CHECK(particleEnergy(flip.particles()) == doctest::Approx(E0).epsilon(1e-12));
    CHECK(particleEnergy(pic.particles()) < 0.9 * E0);
}

TEST_CASE("Tourbillon en boîte pleine : FLIP garde plus d'énergie que PIC, mais il est plus bruité") {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    prm.projectionTol = 1e-8;
    auto run = [&](TransferScheme s, real alpha) {
        FlipParams q = prm;
        q.scheme = s;
        q.flipRatio = alpha;
        FlipSolver solver(q, ParticleScene::BoxVortex);
        const real E0 = particleEnergy(solver.particles());
        for (int k = 0; k < 60; ++k) solver.step(0.02);
        return std::pair<real, real>{particleEnergy(solver.particles()) / E0, solver.noise().relative};
    };
    const auto [ePic, nPic] = run(TransferScheme::PIC, 0);
    const auto [e95, n95] = run(TransferScheme::FLIP, 0.95);
    const auto [e99, n99] = run(TransferScheme::FLIP, 0.99);
    const auto [e1, n1] = run(TransferScheme::FLIP, 1);
    MESSAGE("E/E0 apres 1,2 s : PIC " << ePic << ", FLIP 0.95 " << e95 << ", 0.99 " << e99 << ", 1 " << e1);
    MESSAGE("bruit relatif   : PIC " << nPic << ", FLIP 0.95 " << n95 << ", 0.99 " << n99 << ", 1 " << n1);
    CHECK(ePic < e95);
    CHECK(e95 < e99);
    CHECK(e99 <= e1 * (1 + 1e-9));
    CHECK(e1 <= 1 + 1e-6);          // FLIP ne crée pas d'énergie dans une boîte fermée
    CHECK(nPic < n95);
    CHECK(n95 < n1);
}

// La rupture de barrage (surface libre du jour 8) est testée dans test_dam_break.cpp.
