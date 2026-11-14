// Tests du jour 7 : ensemencement, transferts PIC sur grille MAC, solveur PIC.
#include "doctest.h"
#include <cmath>
#include <random>
#include "core/Kernels.h"
#include "diagnostics/ParticleDiagnostics.h"
#include "numerics/MacOperators.h"
#include "numerics/MacTransfer.h"
#include "scenes/Seeding.h"
#include "solvers/FlipSolver.h"

namespace {

// Vitesses aléatoires (graine fixe) sur les particules.
void randomVelocities(ParticleSet2D& p, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<real> U(-1, 1);
    for (Vec2& v : p.v) v = {U(rng), U(rng)};
}

// Disque de particules centré, loin des parois : tous les pochoirs sont complets.
ParticleSet2D interiorDisk(int n, int perAxis = 2, real jitter = 1, unsigned seed = 7) {
    ParticleSet2D p;
    SeedOptions opt;
    opt.perAxis = perAxis;
    opt.jitter = jitter;
    opt.seed = seed;
    seedParticles(p, n, n, real(1) / n, region::disk({0.5, 0.5}, 0.3), opt);
    return p;
}

bool allFinite(const ParticleSet2D& p) {
    for (std::size_t k = 0; k < p.size(); ++k)
        if (!std::isfinite(p.x[k].x) || !std::isfinite(p.x[k].y) || !std::isfinite(p.v[k].x) ||
            !std::isfinite(p.v[k].y))
            return false;
    return true;
}

} // namespace

// ---------------------------------------------------------------------------
//  Ensemencement
// ---------------------------------------------------------------------------

TEST_CASE("Ensemencement : k² particules par cellule, chacune dans sa sous-cellule, masse = ρ × aire") {
    const int nx = 20, ny = 12;
    const real h = 0.05;
    for (int k : {1, 2, 3}) {
        CAPTURE(k);
        ParticleSet2D p;
        SeedOptions opt;
        opt.perAxis = k;
        opt.density = 3;
        // Rectangle aligné sur la grille : cellules 4..13 en x, 2..7 en y (10 × 6 cellules).
        const int added = seedParticles(p, nx, ny, h, region::box({4 * h, 2 * h}, {14 * h, 8 * h}), opt);
        CHECK(added == 10 * 6 * k * k);
        CHECK(p.size() == std::size_t(added));
        CHECK(particleMass(p) == doctest::Approx(3 * (10 * h) * (6 * h)).epsilon(1e-12));

        // Chaque sous-cellule (h/k de côté) contient exactement une particule.
        std::vector<int> hits(std::size_t(nx * k) * (ny * k), 0);
        for (const Vec2& x : p.x) {
            const int a = int(std::floor(x.x / (h / k))), b = int(std::floor(x.y / (h / k)));
            REQUIRE(a >= 0);
            REQUIRE(b >= 0);
            ++hits[std::size_t(a) + std::size_t(nx * k) * b];
        }
        for (int c : hits) CHECK(c <= 1);
    }
}

TEST_CASE("Ensemencement : même graine, mêmes particules ; graine différente, autres particules") {
    auto seedWith = [](unsigned s) {
        ParticleSet2D p;
        SeedOptions opt;
        opt.seed = s;
        seedParticles(p, 16, 16, 1.0 / 16, region::disk({0.5, 0.5}, 0.3), opt);
        return p;
    };
    const ParticleSet2D a = seedWith(42), b = seedWith(42), c = seedWith(43);
    REQUIRE(a.size() == b.size());
    REQUIRE(a.size() == c.size());          // la région ne dépend pas du tirage
    bool same = true, differ = false;
    for (std::size_t k = 0; k < a.size(); ++k) {
        same = same && a.x[k].x == b.x[k].x && a.x[k].y == b.x[k].y;   // égalité au bit près
        differ = differ || a.x[k].x != c.x[k].x;
    }
    CHECK(same);
    CHECK(differ);
}

TEST_CASE("Ensemencement : jitter = 0 donne le réseau régulier des centres de sous-cellules") {
    ParticleSet2D p;
    SeedOptions opt;
    opt.jitter = 0;
    seedParticles(p, 4, 4, 0.25, region::everywhere(), opt);
    for (const Vec2& x : p.x) {
        // centres en (m + ½) h/2
        CHECK(std::fmod(x.x / 0.125 - 0.5, 1.0) == doctest::Approx(0.0));
        CHECK(std::fmod(x.y / 0.125 - 0.5, 1.0) == doctest::Approx(0.0));
    }
}

// ---------------------------------------------------------------------------
//  Poids d'interpolation
// ---------------------------------------------------------------------------

TEST_CASE("Pochoir linéaire : poids positifs, partition de l'unité, sur u comme sur v") {
    MacGrid2D g(16, 10, 1.0 / 16);
    std::mt19937 rng(1);
    std::uniform_real_distribution<real> X(0.1, 0.9), Y(0.1, 0.52);
    for (int t = 0; t < 200; ++t) {
        const Vec2 x{X(rng), Y(rng)};
        for (const Field2D* f : {&g.u, &g.v}) {
            const Stencil2D s = linearStencil(*f, x);
            real sum = 0;
            for (int b = 0; b < 2; ++b)
                for (int a = 0; a < 2; ++a) {
                    const real w = s.wx[a] * s.wy[b];
                    CHECK(w >= 0);
                    CHECK(inField(*f, s.i0 + a, s.j0 + b));      // point intérieur
                    // le poids est bien N((x - x_f)/h) N((y - y_f)/h)
                    const Vec2 xf = f->position(s.i0 + a, s.j0 + b);
                    CHECK(w == doctest::Approx(LinearKernel::N((x.x - xf.x) / f->h()) *
                                               LinearKernel::N((x.y - xf.y) / f->h())));
                    sum += w;
                }
            CHECK(sum == doctest::Approx(1.0).epsilon(1e-14));
        }
    }
}

TEST_CASE("G2P coïncide avec l'interpolation bilinéaire du jour 3 (Field2D::sample) à l'intérieur") {
    MacGrid2D g(24, 24, 1.0 / 24);
    std::mt19937 rng(2);
    std::uniform_real_distribution<real> U(-1, 1), X(0.05, 0.95);
    for (real& a : g.u.data()) a = U(rng);
    for (real& a : g.v.data()) a = U(rng);
    for (int t = 0; t < 200; ++t) {
        const Vec2 x{X(rng), X(rng)};
        const Vec2 a = interpolateMac(g, x), b = g.velocity(x);
        CHECK(a.x == doctest::Approx(b.x).epsilon(1e-13));
        CHECK(a.y == doctest::Approx(b.y).epsilon(1e-13));
    }
}

// ---------------------------------------------------------------------------
//  Transferts : ce qui est exact
// ---------------------------------------------------------------------------

TEST_CASE("Champ uniforme : l'aller-retour P2G puis G2P est exact, parois comprises") {
    const int n = 32;
    MacGrid2D g(n, n, 1.0 / n);
    ParticleSet2D p;
    seedParticles(p, n, n, g.h(), region::everywhere());     // toute la boîte, y compris contre les parois
    const Vec2 U{0.7, -1.3};
    for (Vec2& v : p.v) v = U;

    MacMasses m;
    particlesToMac(p, g, m);
    CHECK(m.emptyU == 0);
    CHECK(m.emptyV == 0);
    for (real a : g.u.data()) CHECK(a == doctest::Approx(U.x).epsilon(1e-13));
    for (real a : g.v.data()) CHECK(a == doctest::Approx(U.y).epsilon(1e-13));

    macToParticlesPIC(g, p);
    real err = 0;
    for (const Vec2& v : p.v) err = std::max(err, norm(v - U));
    CHECK(err < 1e-12);
}

TEST_CASE("Quantité de mouvement conservée par P2G et par G2P") {
    const int n = 32;
    MacGrid2D g(n, n, 1.0 / n);
    ParticleSet2D p = interiorDisk(n);
    randomVelocities(p, 3);

    const Vec2 P0 = particleMomentum(p);
    MacMasses m;
    particlesToMac(p, g, m);
    const Vec2 Pg = gridMomentum(g, m);
    macToParticlesPIC(g, p);
    const Vec2 P1 = particleMomentum(p);

    const real scale = particleMass(p);          // |u| <= 1, donc |P| <= M
    CHECK(norm(Pg - P0) < 1e-12 * scale);
    CHECK(norm(P1 - Pg) < 1e-12 * scale);
}

TEST_CASE("Énergie : E_particules >= E_grille >= E_après (dissipation de PIC)") {
    const int n = 32;
    for (unsigned seed : {4u, 5u, 6u}) {
        CAPTURE(seed);
        MacGrid2D g(n, n, 1.0 / n);
        ParticleSet2D p = interiorDisk(n, 2, 1, seed);
        randomVelocities(p, seed);

        const real E0 = particleEnergy(p);
        MacMasses m;
        particlesToMac(p, g, m);
        const real Eg = gridEnergy(g, m);
        macToParticlesPIC(g, p);
        const real E1 = particleEnergy(p);

        CHECK(Eg <= E0 * (1 + 1e-12));
        CHECK(E1 <= Eg * (1 + 1e-12));
        CHECK(E1 < 0.9 * E0);          // un champ aléatoire perd beaucoup : il oscille à l'échelle h
        MESSAGE("champ aleatoire : E_grille/E0 = " << Eg / E0 << ", E_apres/E0 = " << E1 / E0);
    }
}

// ---------------------------------------------------------------------------
//  Transferts : ce que PIC perd
// ---------------------------------------------------------------------------

TEST_CASE("Rotation rigide : PIC perd du moment cinétique à chaque aller-retour") {
    const int n = 32;
    const Vec2 c{0.5, 0.5};
    MacGrid2D g(n, n, 1.0 / n);
    ParticleSet2D p = interiorDisk(n);
    for (std::size_t k = 0; k < p.size(); ++k) {
        const Vec2 r = p.x[k] - c;
        p.v[k] = {-r.y, r.x};                    // ω = 1
    }
    const real L0 = particleAngularMomentum(p, c), E0 = particleEnergy(p);
    const Vec2 P0 = particleMomentum(p);         // ≈ 0, pas exactement : le tirage n'est pas symétrique
    MacMasses m;
    real Lprev = L0;
    for (int it = 1; it <= 10; ++it) {
        particlesToMac(p, g, m);
        macToParticlesPIC(g, p);
        const real L = particleAngularMomentum(p, c);
        CHECK(L < Lprev);                        // strictement
        Lprev = L;
    }
    CHECK(particleEnergy(p) < E0);
    CHECK(norm(particleMomentum(p) - P0) < 1e-12 * particleMass(p));   // P, elle, est conservée
    MESSAGE("rotation rigide, n = 32 : L/L0 apres 10 allers-retours = " << Lprev / L0);
}

TEST_CASE("Champ affine : l'erreur de PIC est exactement A · Σ_f w_fp (x̄_f - x_f)") {
    // Pour u(x) = A x + b, P2G donne u_f = u(x̄_f), x̄_f barycentre pondéré des
    // particules autour de la face ; l'interpolation bilinéaire reproduisant les
    // fonctions affines, G2P donne u_p' = u(x_p) + A · Σ_f w_fp (x̄_f - x_f).
    const int n = 32;
    MacGrid2D g(n, n, 1.0 / n);
    ParticleSet2D p = interiorDisk(n);
    const real A[2][2] = {{0.3, -1.1}, {0.8, 0.4}};      // quelconque, pas seulement une rotation
    const Vec2 b{0.2, -0.5};
    auto field = [&](Vec2 x) { return Vec2{A[0][0] * x.x + A[0][1] * x.y + b.x, A[1][0] * x.x + A[1][1] * x.y + b.y}; };
    for (std::size_t k = 0; k < p.size(); ++k) p.v[k] = field(p.x[k]);

    // Barycentres pondérés x̄_f, sur les faces u et sur les faces v.
    auto barycenters = [&](const Field2D& f, std::vector<Vec2>& xbar) {
        std::vector<real> wsum(f.data().size(), 0);
        xbar.assign(f.data().size(), Vec2{});
        for (std::size_t k = 0; k < p.size(); ++k) {
            const Stencil2D s = linearStencil(f, p.x[k]);
            for (int bb = 0; bb < 2; ++bb)
                for (int a = 0; a < 2; ++a) {
                    const std::size_t id = std::size_t(s.i0 + a) + std::size_t(f.nx()) * (s.j0 + bb);
                    const real w = s.wx[a] * s.wy[bb] * p.m[k];
                    wsum[id] += w;
                    xbar[id] += w * p.x[k];
                }
        }
        for (std::size_t id = 0; id < xbar.size(); ++id)
            if (wsum[id] > 0) xbar[id] = xbar[id] / wsum[id];
    };
    std::vector<Vec2> xbarU, xbarV;
    barycenters(g.u, xbarU);
    barycenters(g.v, xbarV);

    const ParticleSet2D before = p;
    MacMasses m;
    particlesToMac(p, g, m);
    macToParticlesPIC(g, p);

    real errMax = 0, errPred = 0;
    for (std::size_t k = 0; k < p.size(); ++k) {
        // Décalage moyen δ = Σ_f w_fp (x̄_f - x_f), séparément pour les faces u et v.
        auto delta = [&](const Field2D& f, const std::vector<Vec2>& xbar) {
            const Stencil2D s = linearStencil(f, before.x[k]);
            Vec2 d;
            for (int bb = 0; bb < 2; ++bb)
                for (int a = 0; a < 2; ++a) {
                    const int i = s.i0 + a, j = s.j0 + bb;
                    d += (s.wx[a] * s.wy[bb]) * (xbar[std::size_t(i) + std::size_t(f.nx()) * j] - f.position(i, j));
                }
            return d;
        };
        const Vec2 dU = delta(g.u, xbarU), dV = delta(g.v, xbarV);
        const Vec2 predicted{field(before.x[k]).x + A[0][0] * dU.x + A[0][1] * dU.y,
                             field(before.x[k]).y + A[1][0] * dV.x + A[1][1] * dV.y};
        errMax  = std::max(errMax, norm(p.v[k] - before.v[k]));
        errPred = std::max(errPred, norm(p.v[k] - predicted));
    }
    MESSAGE("champ affine : erreur PIC max = " << errMax << ", ecart a la formule = " << errPred);
    CHECK(errMax > 1e-4);            // PIC ne reproduit PAS les champs affines (APIC le fera, jour 10)
    CHECK(errPred < 1e-12);          // mais son erreur est exactement celle prévue
}

// ---------------------------------------------------------------------------
//  Solveur PIC complet (boîte pleine)
// ---------------------------------------------------------------------------

TEST_CASE("Solveur PIC : un fluide au repos reste au repos") {
    FlipParams prm;
    prm.nx = prm.ny = 24;
    FlipSolver s(prm, ParticleScene::BoxAtRest);
    const ParticleSet2D p0 = s.particles();
    for (int k = 0; k < 20; ++k) s.step(0.02);
    real moved = 0;
    for (std::size_t k = 0; k < p0.size(); ++k) moved = std::max(moved, norm(s.particles().x[k] - p0.x[k]));
    // Rien ne bouge, sauf les rares particules tirées à moins de 10⁻³ h d'une paroi,
    // que l'advection ramène à cette marge.
    CHECK(moved <= 1e-3 * s.grid().h());
    CHECK(particleMaxSpeed(s.particles()) == 0.0);
}

TEST_CASE("Solveur PIC : tourbillon en boîte, divergence projetée, particules dans le domaine") {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    prm.projectionTol = 1e-8;
    FlipSolver s(prm, ParticleScene::BoxVortex);
    const real E0 = particleEnergy(s.particles());
    const real M0 = particleMass(s.particles());
    for (int k = 0; k < 60; ++k) {
        s.step(0.02);
        CHECK(s.lastStats().divAfter <= 2 * prm.projectionTol);
    }
    const Vec2 L = s.grid().size();
    bool inside = true;
    for (const Vec2& x : s.particles().x) inside = inside && x.x > 0 && x.x < L.x && x.y > 0 && x.y < L.y;
    CHECK(inside);
    CHECK(allFinite(s.particles()));
    CHECK(particleMass(s.particles()) == M0);          // aucune particule perdue ni créée
    const real E = particleEnergy(s.particles());
    CHECK(E < E0);                                     // la dissipation de PIC
    CHECK(E > 0.2 * E0);
    MESSAGE("tourbillon PIC 32x32 : E(1,2 s)/E0 = " << E / E0);
}

TEST_CASE("Solveur PIC : la souris pousse les particules dans le sens du geste") {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    FlipSolver s(prm, ParticleScene::BoxAtRest);
    FlipBrush b;
    b.active = true;
    b.pos = {0.5, 0.5};
    b.velocity = {2, 0};
    s.setBrush(b);
    s.step(1.0 / 60);
    // vitesse moyenne des particules proches du pinceau
    Vec2 mean;
    int count = 0;
    for (std::size_t k = 0; k < s.particles().size(); ++k)
        if (norm(s.particles().x[k] - b.pos) < 0.03) { mean += s.particles().v[k]; ++count; }
    REQUIRE(count > 0);
    mean = mean / real(count);
    CHECK(mean.x > 0.1);
    CHECK(std::abs(mean.y) < 0.2 * mean.x);
}

TEST_CASE("Scène disque en rotation : transferts seuls, particules immobiles, L décroissant") {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    FlipSolver s(prm, ParticleScene::RotatingDisk);
    const ParticleSet2D p0 = s.particles();
    const real L0 = particleAngularMomentum(p0, s.domainCenter());
    for (int k = 0; k < 5; ++k) s.step(0.01);
    bool still = true;
    for (std::size_t k = 0; k < p0.size(); ++k) still = still && s.particles().x[k].x == p0.x[k].x;
    CHECK(still);
    CHECK(particleAngularMomentum(s.particles(), s.domainCenter()) < L0);
}
