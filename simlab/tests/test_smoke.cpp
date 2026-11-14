// Tests du jour 5 : diffusion implicite, forces, boucle Stable Fluids complète.
#include "doctest.h"
#include <cmath>
#include <random>
#include <vector>
#include "diagnostics/FluidDiagnostics.h"
#include "numerics/Diffusion.h"
#include "numerics/Forces.h"
#include "numerics/MacOperators.h"
#include "scenes/SmokeScenes.h"
#include "solvers/StableFluidsSolver.h"

namespace {

bool allFinite(const Field2D& f) {
    for (real a : f.data())
        if (!std::isfinite(a)) return false;
    return true;
}

// Solveur « nu » : aucune force, aucune source, pour isoler un mécanisme.
StableFluidsParams quietParams(int n) {
    StableFluidsParams p;
    p.nx = p.ny = n;
    p.buoyancy = 0;
    p.emitterOn = false;
    p.projectionTol = 1e-10;
    return p;
}

} // namespace

// ---------------------------------------------------------------------------
//  Diffusion implicite
// ---------------------------------------------------------------------------

TEST_CASE("Helmholtz : symétrique et défini positif, pour les 4 combinaisons de bords") {
    std::mt19937 rng(1);
    std::uniform_real_distribution<real> U(-1, 1);
    for (EdgeCondition ex : {EdgeCondition::Neumann, EdgeCondition::Dirichlet0})
        for (EdgeCondition ey : {EdgeCondition::Neumann, EdgeCondition::Dirichlet0}) {
            HelmholtzOperator A{13, 7, 2.5, ex, ey};
            Vector x(13 * 7), y(13 * 7), Ax, Ay;
            for (real& a : x) a = U(rng);
            for (real& a : y) a = U(rng);
            A(x, Ax);
            A(y, Ay);
            CHECK(dotProduct(x, Ay) == doctest::Approx(dotProduct(Ax, y)).epsilon(1e-12));
            CHECK(dotProduct(x, Ax) >= dotProduct(x, x));   // I - cΔ >= I
        }
}

TEST_CASE("Diffusion d'un scalaire : un mode de Fourier est amorti EXACTEMENT par 1/(1 + κ dt λ)") {
    // Mode cosinus aux centres, compatible Neumann : q = cos(π m (i+½)/n).
    // Valeur propre discrète de -Δ_h : λ = (4/h²) sin²(π m / (2n)).
    const int n = 32, m = 3;
    const real h = 1.0 / n, kappa = 0.01, dt = 0.5;
    Field2D q(n, n, h, 0.5, 0.5);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) q(i, j) = std::cos(kPi * m * (i + 0.5) / n);
    Field2D q0 = q;

    const real lambda = 4 / (h * h) * std::pow(std::sin(kPi * m / (2.0 * n)), 2);
    const real factor = 1 / (1 + kappa * dt * lambda);
    DiffusionStats st = diffuseCells(q, kappa, dt, 1e-13);
    CHECK(st.converged);
    for (std::size_t k = 0; k < q.data().size(); ++k)
        CHECK(q.data()[k] == doctest::Approx(factor * q0.data()[k]).epsilon(1e-10));
}

TEST_CASE("Diffusion d'un scalaire : quantité conservée, principe du maximum, stable à grand pas") {
    const int n = 40;
    Field2D q(n, n, 1.0 / n, 0.5, 0.5);
    std::mt19937 rng(2);
    std::uniform_real_distribution<real> U(0, 1);
    for (real& a : q.data()) a = U(rng);
    const real m0 = totalAmount(q), mx0 = q.maxValue(), mn0 = q.minValue();

    for (real dt : {1e-3, 1.0, 1e3}) {                      // κ dt / h² jusqu'à 1,6·10⁶
        CAPTURE(dt);
        Field2D r = q;
        DiffusionStats st = diffuseCells(r, 1e-3, dt);
        CHECK(st.converged);
        CHECK(allFinite(r));
        CHECK(totalAmount(r) == doctest::Approx(m0).epsilon(1e-9));     // Neumann : flux nul
        CHECK(r.maxValue() <= mx0 + 1e-9);                               // M-matrice
        CHECK(r.minValue() >= mn0 - 1e-9);
    }
}

TEST_CASE("Diffusion de la vitesse : le tourbillon cellulaire est un mode propre exact") {
    // u = sin(πx) cos(πy) aux faces u : Dirichlet en x, Neumann en y ;
    // valeur propre discrète λ = 2 · (4/h²) sin²(πh/2), la même pour v.
    const int n = 48;
    const real h = 1.0 / n, nu = 1e-3, dt = 0.2;
    MacGrid2D g(n, n, h);
    smoke_scene::setCellularVortex(g, 1);
    MacGrid2D g0 = g;

    const real lambda = 8 / (h * h) * std::pow(std::sin(kPi * h / 2), 2);
    const real factor = 1 / (1 + nu * dt * lambda);
    DiffusionStats st = diffuseVelocity(g, nu, dt, 1e-13);
    CHECK(st.converged);
    // Écart maximal relatif à l'amplitude (les faces proches de sin(π) valent ~1e-16).
    real err = 0;
    for (std::size_t k = 0; k < g.u.data().size(); ++k)
        err = std::max(err, std::abs(g.u.data()[k] - factor * g0.u.data()[k]));
    for (std::size_t k = 0; k < g.v.data().size(); ++k)
        err = std::max(err, std::abs(g.v.data()[k] - factor * g0.v.data()[k]));
    CHECK(err <= 1e-10 * g0.maxSpeed());
    CHECK(maxAbsDivergence(g) < 1e-10);       // la diffusion garde ce champ à divergence nulle
}

TEST_CASE("Diffusion de la vitesse : parois intactes et énergie décroissante") {
    const int n = 32;
    MacGrid2D g(n, n, 1.0 / n);
    std::mt19937 rng(3);
    std::uniform_real_distribution<real> U(-1, 1);
    for (real& a : g.u.data()) a = U(rng);
    for (real& a : g.v.data()) a = U(rng);
    enforceWallBC(g);
    const real e0 = kineticEnergy(g);

    diffuseVelocity(g, 1e-2, 0.1);
    for (int j = 0; j < n; ++j) { CHECK(g.u(0, j) == 0); CHECK(g.u(n, j) == 0); }
    for (int i = 0; i < n; ++i) { CHECK(g.v(i, 0) == 0); CHECK(g.v(i, n) == 0); }
    CHECK(kineticEnergy(g) < e0);
}

TEST_CASE("Dissipation : la masse est multipliée exactement par 1/(1 + α dt)") {
    Field2D q(16, 16, 1.0 / 16, 0.5, 0.5, 2.0);
    dissipate(q, 0.5, 0.2);
    CHECK(totalAmount(q) == doctest::Approx(2.0 / 1.1));
}

// ---------------------------------------------------------------------------
//  Forces
// ---------------------------------------------------------------------------

TEST_CASE("Flottabilité uniforme : c'est un gradient, la projection l'annule (fluide au repos)") {
    // f = β s ŷ avec s uniforme : f = ∇(β s y). Le fluide doit rester immobile,
    // la pression absorbe tout (équilibre hydrostatique). Test de la cohérence
    // force + parois + projection.
    StableFluidsParams p = quietParams(32);
    p.buoyancy = 5;
    StableFluidsSolver s(p, SmokeScene::Plume);
    s.smoke().fill(1);
    for (int k = 0; k < 20; ++k) s.step(0.05);
    CHECK(s.grid().maxSpeed() < 1e-8);
    // et la pression est bien hydrostatique : p(y) croissante vers le haut, pente β s
    const MacGrid2D& g = s.grid();
    const real slope = (g.p(10, 20) - g.p(10, 19)) / g.h();
    CHECK(slope == doctest::Approx(5.0).epsilon(1e-6));
}

TEST_CASE("Flottabilité : une bulle de fumée monte") {
    StableFluidsParams p = quietParams(48);
    p.buoyancy = 4;
    StableFluidsSolver s(p, SmokeScene::Bubble);
    auto centroidY = [&]() {
        const Field2D& q = s.smoke();
        real m = 0, my = 0;
        for (int j = 0; j < q.ny(); ++j)
            for (int i = 0; i < q.nx(); ++i) { m += q(i, j); my += q(i, j) * q.position(i, j).y; }
        return my / m;
    };
    const real y0 = centroidY();
    for (int k = 0; k < 30; ++k) s.step(1.0 / 60);
    CHECK(centroidY() > y0 + 2 * s.grid().h());
    CHECK(s.lastStats().divAfter <= 2 * p.projectionTol);
}

TEST_CASE("Souris : le pinceau pousse le fluide dans sa direction et injecte de la fumée") {
    StableFluidsParams p = quietParams(48);
    StableFluidsSolver s(p, SmokeScene::Plume);
    Brush b;
    b.active = true;
    b.pos = {0.5, 0.5};
    b.velocity = {2, 0};
    s.setBrush(b);
    s.step(1.0 / 60);
    CHECK(s.grid().velocity({0.5, 0.5}).x > 0.1);
    CHECK(std::abs(s.grid().velocity({0.5, 0.5}).y) < 1e-3);
    CHECK(s.smokeAmount() > 0);
    CHECK(s.lastStats().divAfter <= 2 * p.projectionTol);
}

// ---------------------------------------------------------------------------
//  Validations du planning (J5)
// ---------------------------------------------------------------------------

TEST_CASE("J5 : énergie décroissante avec ν = 0 (dissipation numérique mesurée)") {
    // Le tourbillon cellulaire est une solution stationnaire d'Euler :
    // l'énergie exacte est constante, toute perte est numérique.
    StableFluidsParams p = quietParams(32);
    StableFluidsSolver s(p, SmokeScene::Vortex);
    const real e0 = s.kineticEnergy();
    real prev = e0;
    bool monotone = true;
    for (int k = 0; k < 100; ++k) {
        s.step(0.01);
        const real e = s.kineticEnergy();
        monotone = monotone && e <= prev * (1 + 1e-12);
        prev = e;
    }
    MESSAGE("E(1 s)/E(0) = " << prev / e0 << " avec nu = 0 (perte purement numerique)");
    CHECK(monotone);
    CHECK(prev < e0);
    CHECK(prev > 0.5 * e0);      // dissipation présente, mais pas absurde
}

TEST_CASE("J5 : décroissance plus rapide avec ν > 0, proche du facteur visqueux exact") {
    const real T = 1.0, dt = 0.01, nu = 2e-3;
    auto run = [&](real viscosity) {
        StableFluidsParams p = quietParams(32);
        p.viscosity = viscosity;
        StableFluidsSolver s(p, SmokeScene::Vortex);
        const real e0 = s.kineticEnergy();
        for (int k = 0; k < int(std::lround(T / dt)); ++k) s.step(dt);
        return s.kineticEnergy() / e0;
    };
    const real r0 = run(0), rNu = run(nu);
    const real k2 = 2 * kPi * kPi;                       // domaine unité
    const real viscousFactor = std::exp(-2 * nu * k2 * T);
    MESSAGE("E/E0 : nu = 0 -> " << r0 << ", nu = " << nu << " -> " << rNu
            << " ; rapport " << rNu / r0 << " (exact visqueux : " << viscousFactor << ")");
    CHECK(rNu < r0);
    // Les deux dissipations se composent : le rapport doit approcher exp(-2 ν k² T).
    CHECK(rNu / r0 == doctest::Approx(viscousFactor).epsilon(0.05));
}

TEST_CASE("J5 : rien n'explose à grand pas de temps") {
    for (real dt : {0.05, 0.2, 1.0}) {
        CAPTURE(dt);
        // (1) Tourbillon libre : l'énergie ne doit jamais croître, quel que soit dt.
        {
            StableFluidsParams p = quietParams(32);
            StableFluidsSolver s(p, SmokeScene::Vortex);
            const real e0 = s.kineticEnergy();
            real emax = 0;
            for (int k = 0; k < 20; ++k) { s.step(dt); emax = std::max(emax, s.kineticEnergy()); }
            CHECK(s.lastStats().cfl > 0.5);              // on est bien à grand pas
            CHECK(allFinite(s.grid().u));
            CHECK(allFinite(s.grid().v));
            CHECK(emax <= e0 * (1 + 1e-12));
        }
        // (2) Panache forcé (flottabilité + émetteur + viscosité + diffusion) :
        //     valeurs finies, vitesse bornée, fumée positive.
        {
            StableFluidsParams p;
            p.nx = p.ny = 32;
            p.viscosity = 1e-4;
            p.diffusivity = 1e-4;
            p.dissipation = 0.1;
            StableFluidsSolver s(p, SmokeScene::Plume);
            for (int k = 0; k < 40; ++k) s.step(dt);
            CHECK(allFinite(s.grid().u));
            CHECK(allFinite(s.grid().v));
            CHECK(allFinite(s.smoke()));
            CHECK(s.grid().maxSpeed() < 10);             // ~ quelques m/s au plus
            CHECK(s.smoke().minValue() >= -1e-9);        // l'advection est convexe
        }
    }
}
