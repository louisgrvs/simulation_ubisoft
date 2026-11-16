// ============================================================================
//  Tests du jour 8 : surface libre (Zhu-Bridson §4.2.3 et §4.2.5)
// ============================================================================
//  Du plus simple au plus complet :
//    1. briques : étiquettes, faces fluides, distance signée, repoussement
//    2. projection : non-régression, colonne hydrostatique EXACTE sur la grille seule
//    3. extrapolation et déplacement en sous-pas
//    4. solveur complet : cuve au repos, cuve agitée (volume), obstacle
//  Si la colonne hydrostatique (2) passe, toute erreur dans la cuve (4) vient des
//  étiquettes, des transferts ou de l'extrapolation, pas de Poisson.
// ============================================================================
#include "doctest.h"
#include <cmath>
#include "diagnostics/FreeSurfaceDiagnostics.h"
#include "diagnostics/ParticleDiagnostics.h"
#include "geometry/SolidSdf2D.h"
#include "numerics/Extrapolation.h"
#include "numerics/Forces.h"
#include "numerics/FreeSurface.h"
#include "numerics/MacOperators.h"
#include "numerics/ParticleAdvection.h"
#include "numerics/Projection.h"
#include "solvers/FlipSolver.h"

namespace {

constexpr real kG = real(9.81);

// Tourbillon à divergence nulle dans la moitié basse [0, 1] × [0, 0.5] de la cuve,
// tangent au fond, aux murs et à la surface y = 0.5 : ψ = A/(2π) sin(πx) sin(2πy).
void stirTank(FlipSolver& s, real A) {
    ParticleSet2D& p = s.particles();
    for (std::size_t k = 0; k < p.size(); ++k) {
        const Vec2 x = p.x[k];
        p.v[k] = {A * std::sin(kPi * x.x) * std::cos(2 * kPi * x.y),
                  -A * real(0.5) * std::cos(kPi * x.x) * std::sin(2 * kPi * x.y)};
    }
}

bool particlesFinite(const ParticleSet2D& p) {
    for (std::size_t k = 0; k < p.size(); ++k)
        if (!std::isfinite(p.x[k].x) || !std::isfinite(p.x[k].y) || !std::isfinite(p.v[k].x) || !std::isfinite(p.v[k].y))
            return false;
    return true;
}

ParticleSet2D oneParticleAt(Vec2 x) {
    ParticleSet2D p;
    p.add(x, {0, 0}, 1);
    return p;
}

} // namespace

// ----------------------------------------------------------------------------
//  1. Briques
// ----------------------------------------------------------------------------

TEST_CASE("Etiquettes : particule -> Fluid, vide -> Air, solide -> Solid, hors domaine -> Solid") {
    const int n = 8;
    const real h = real(1) / n;
    CellLabels L(n, n);
    SolidSdf2D sdf(Vec2{1, 1});
    sdf.addBox(Vec2{0.5, 0.5}, Vec2{0.75, 0.75});   // cellules (4..5, 4..5)

    ParticleSet2D p = oneParticleAt({0.1, 0.1});      // cellule (0, 0)
    p.add({0.6, 0.6}, {0, 0}, 1);                     // dans l'obstacle : la cellule reste Solid
    labelCells(L, p, h, &sdf);

    CHECK(L(0, 0) == CellType::Fluid);
    CHECK(L(3, 3) == CellType::Air);
    CHECK(L(4, 4) == CellType::Solid);
    CHECK(L(5, 5) == CellType::Solid);
    CHECK(L(6, 6) == CellType::Air);
    CHECK(L(-1, 0) == CellType::Solid);
    CHECK(L(0, n) == CellType::Solid);
    CHECK(L.count(CellType::Fluid) == 1);
    CHECK(L.count(CellType::Solid) == 4);
}

TEST_CASE("Etiquettes : une petite bulle enfermee devient Fluid, la grande zone d'air reste Air") {
    const int n = 8;
    const real h = real(1) / n;
    ParticleSet2D p;
    for (int j = 0; j < 4; ++j)                        // 4 rangées d'eau…
        for (int i = 0; i < n; ++i)
            if (!(i == 3 && j == 1)) p.add({(i + real(0.5)) * h, (j + real(0.5)) * h}, {0, 0}, 1);   // …avec un trou en (3, 1)

    CellLabels L(n, n);
    labelCells(L, p, h, nullptr, 0);                   // sans comblement : le trou est de l'air
    CHECK(L(3, 1) == CellType::Air);
    labelCells(L, p, h, nullptr);                      // comblement par défaut (≤ 8 cellules)
    CHECK(L(3, 1) == CellType::Fluid);
    CHECK(L(3, 5) == CellType::Air);                   // l'air au-dessus (32 cellules) n'est pas une bulle
    CHECK(L.count(CellType::Fluid) == 4 * n);
}

TEST_CASE("Faces fluides : une cellule Fluid et aucune Solid") {
    CellLabels L(3, 1);
    L.set(0, 0, CellType::Fluid);
    L.set(1, 0, CellType::Air);
    L.set(2, 0, CellType::Solid);
    // Faces u(i, 0) : entre (i−1, 0) et (i, 0).
    CHECK_FALSE(uFaceIsFluid(L, 0, 0));   // mur / Fluid
    CHECK(uFaceIsFluid(L, 1, 0));         // Fluid / Air : la surface libre
    CHECK_FALSE(uFaceIsFluid(L, 2, 0));   // Air / Solid
    CHECK_FALSE(uFaceIsFluid(L, 3, 0));   // Solid / mur
    L.set(1, 0, CellType::Fluid);
    CHECK(uFaceIsFluid(L, 1, 0));         // Fluid / Fluid
    CHECK_FALSE(uFaceIsFluid(L, 2, 0));   // Fluid / Solid
}

TEST_CASE("Distance signee : murs, obstacle, normale") {
    SolidSdf2D s(Vec2{2, 1});
    CHECK(s.phi({1, 0.5}) == doctest::Approx(0.5));    // centre : à 0,5 du sol et du plafond
    CHECK(s.phi({0.1, 0.5}) == doctest::Approx(0.1));  // près du mur gauche
    CHECK(s.phi({-0.2, 0.5}) == doctest::Approx(-0.2)); // dans le mur
    const Vec2 n = s.normal({0.1, 0.5});
    CHECK(n.x == doctest::Approx(1));
    CHECK(n.y == doctest::Approx(0).epsilon(1e-9));

    s.addBox(Vec2{0.8, 0}, Vec2{1.2, 0.3});           // marche posée au fond
    CHECK(s.phi({1, 0.2}) == doctest::Approx(-0.1));   // dedans, à 0,1 du dessus
    CHECK(s.phi({1, 0.4}) == doctest::Approx(0.1));    // au-dessus, à 0,1
    CHECK(s.normal({1, 0.4}).y == doctest::Approx(1));
}

TEST_CASE("Repoussement : retour a la marge, vitesse entrante annulee, tangentielle gardee") {
    SolidSdf2D s(Vec2{1, 1});
    ParticleSet2D p;
    p.add({-0.05, 0.5}, {-2, 3}, 1);    // dans le mur gauche, va vers le mur
    p.add({0.5, 0.5}, {1, 1}, 1);       // loin des murs : intouchée
    const real margin = 1e-3;
    pushOutOfSolids(p, s, margin);
    CHECK(s.phi(p.x[0]) == doctest::Approx(margin));
    CHECK(p.v[0].x == doctest::Approx(0));   // composante normale entrante retirée
    CHECK(p.v[0].y == doctest::Approx(3));   // glissement le long du mur
    CHECK(p.x[1].x == 0.5);
    CHECK(p.v[1].x == 1);
}

// ----------------------------------------------------------------------------
//  2. Projection
// ----------------------------------------------------------------------------

TEST_CASE("Projection surface libre : toutes les cellules Fluid -> meme resultat que le J4") {
    const int n = 16;
    MacGrid2D a(n, n, real(1) / n);
    a.setVelocity([](Vec2 x) { return Vec2{std::sin(3 * x.y) + x.x * x.x, std::cos(2 * x.x) * x.y}; });
    MacGrid2D b = a;
    CellLabels L(n, n);
    L.fill(CellType::Fluid);

    ProjectionOptions opt;
    opt.tol = 1e-12;
    opt.warmStart = false;
    project(a, real(0.01), opt);
    const ProjectionStats st = project(b, L, real(0.01), opt);

    real du = 0, dv = 0;
    for (std::size_t k = 0; k < a.u.data().size(); ++k) du = std::max(du, std::abs(a.u.data()[k] - b.u.data()[k]));
    for (std::size_t k = 0; k < a.v.data().size(); ++k) dv = std::max(dv, std::abs(a.v.data()[k] - b.v.data()[k]));
    CHECK(du < 1e-10);
    CHECK(dv < 1e-10);
    CHECK(st.divAfter <= 2e-12);
}

TEST_CASE("Colonne hydrostatique sur la grille seule : p = rho g (y_s - y) et v = 0, exactement") {
    const int n = 16, k = 8;                 // 8 rangées d'eau, de l'air au-dessus
    const real h = real(1) / n, dt = real(0.01), rho = 1000;
    MacGrid2D g(n, n, h);
    CellLabels L(n, n);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) L.set(i, j, j < k ? CellType::Fluid : CellType::Air);

    addGravity(g, Vec2{0, -kG}, dt);         // v = −g dt partout
    ProjectionOptions opt;
    opt.tol = 1e-12;
    opt.density = rho;
    const ProjectionStats st = project(g, L, dt, opt);
    CHECK(st.converged);
    CHECK(st.divAfter <= 2e-12);

    // Surface effective au centre de la première cellule d'air : y_s − y_j = (k − j) h.
    real maxRelErr = 0;
    for (int j = 0; j < k; ++j)
        for (int i = 0; i < n; ++i) {
            const real pRef = rho * kG * (k - j) * h;
            maxRelErr = std::max(maxRelErr, std::abs(g.p(i, j) - pRef) / (rho * kG * k * h));
        }
    CHECK(maxRelErr < 1e-9);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) CHECK(g.p(i, j) == (j < k ? g.p(i, j) : 0.0));   // p = 0 dans l'air

    real vmax = 0;                           // l'eau ne bouge pas : la pression compense la gravité
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i < n; ++i)
            if (vFaceIsFluid(L, i, j)) vmax = std::max(vmax, std::abs(g.v(i, j)));
    CHECK(vmax < 1e-10);

    real pRef;
    REQUIRE(hydrostaticReference(L, h, 3, 0, rho, kG, pRef));
    CHECK(pRef == doctest::Approx(rho * kG * k * h));
}

// ----------------------------------------------------------------------------
//  3. Extrapolation et déplacement
// ----------------------------------------------------------------------------

TEST_CASE("Extrapolation couche par couche : l'exemple 3 x 3") {
    Field2D f(3, 3, 1, 0, 0);
    FaceMask valid(9, 0);
    for (int i = 0; i < 3; ++i) { f(i, 0) = real(i + 1); valid[i] = 1; }   // rangée du bas : 1, 2, 3

    FaceMask m1 = valid;
    Field2D f1 = f;
    CHECK(extrapolateLayers(f1, m1, 1) == 1);
    for (int i = 0; i < 3; ++i) {
        CHECK(f1(i, 1) == real(i + 1));      // chaque face copie son voisin du bas…
        CHECK(m1[i + 3] == 1);
        CHECK(m1[i + 6] == 0);               // …et la rangée du haut attend la couche 2 (style Jacobi)
    }
    CHECK(extrapolateLayers(f, valid, 10) == 2);   // tout est rempli en 2 couches
    for (int i = 0; i < 3; ++i) CHECK(f(i, 2) == real(i + 1));
}

TEST_CASE("Extrapolation : moyenne des voisins valides, champ constant reproduit") {
    Field2D f(3, 3, 1, 0, 0);
    FaceMask valid(9, 0);
    f(1, 0) = 2; valid[1] = 1;              // en dessous de (1, 1)
    f(2, 1) = 4; valid[5] = 1;              // à droite de (1, 1)
    extrapolateLayers(f, valid, 1);
    CHECK(f(1, 1) == 3);                     // (2 + 4) / 2

    Field2D c(5, 4, 1, 0, 0, real(7));
    FaceMask mc(20, 0);
    mc[0] = 1;
    extrapolateLayers(c, mc, 20);
    for (real x : c.data()) CHECK(x == 7);
}

TEST_CASE("Deplacement RK2 : champ uniforme exact, nombre de sous-pas CFL") {
    const int n = 32;
    const real h = real(1) / n;
    MacGrid2D g(n, n, h);
    g.setVelocity([](Vec2) { return Vec2{real(0.3), real(-0.1)}; });
    SolidSdf2D s(g.size());
    ParticleSet2D p = oneParticleAt({0.3, 0.6});

    const real dt = real(0.5);                       // 0,3 × 0,5 = 0,15 m ≈ 4,8 cellules
    const AdvectionStats st = advectParticlesRK2(p, g, dt, 1, s, 1e-3 * h);
    CHECK(st.substeps == int(std::ceil(real(0.3) * dt / h)));
    CHECK(p.x[0].x == doctest::Approx(0.3 + 0.15).epsilon(1e-12));
    CHECK(p.x[0].y == doctest::Approx(0.6 - 0.05).epsilon(1e-12));
}

TEST_CASE("Deplacement RK2 : une particule qui fonce dans le mur reste dans la cuve") {
    const int n = 16;
    const real h = real(1) / n;
    MacGrid2D g(n, n, h);
    g.setVelocity([](Vec2) { return Vec2{-3, 0}; });
    SolidSdf2D s(g.size());
    ParticleSet2D p = oneParticleAt({0.2, 0.5});
    advectParticlesRK2(p, g, real(0.5), 1, s, 1e-3 * h);
    CHECK(s.phi(p.x[0]) >= 1e-3 * h * (1 - 1e-9));
}

// ----------------------------------------------------------------------------
//  4. Solveur complet
// ----------------------------------------------------------------------------

TEST_CASE("Cuve au repos : vitesses parasites, pression rho g h et volume") {
    for (const real H : {real(0.5), real(0.5) + real(0.5) / 32})       // surface alignée, puis décalée d'une demi-cellule
        for (const real jitter : {real(0), real(1)}) {
            CAPTURE(H);
            CAPTURE(jitter);
            FlipParams prm;
            prm.nx = prm.ny = 32;
            prm.tankHeight = H;
            prm.jitter = jitter;
            prm.projectionTol = 1e-8;
            FlipSolver s(prm, ParticleScene::TankAtRest);
            s.step(prm.dt);
            const int cells0 = s.lastStats().fluidCells;

            FreeSurfaceStats st;
            real maxSpeed = 0;
            for (int k = 0; k < 300; ++k) {                 // 5 s
                s.step(prm.dt);
                st = measureFreeSurface(s.grid(), s.labels(), s.particles(), prm.density, kG);
                maxSpeed = std::max(maxSpeed, st.maxSpeed);
                REQUIRE(st.fluidCells == cells0);           // volume : pas une cellule de perdue
            }
            CHECK(maxSpeed < 1e-8);                          // sqrt(gH) ≈ 2,2 m/s
            CHECK(st.hydroErr < 1e-8);
            CHECK(st.divAfter <= 2 * prm.projectionTol);
            CHECK(particlesFinite(s.particles()));
        }
}

TEST_CASE("Cuve agitee : le volume de fluide reste stable, aucune particule dans les murs") {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    prm.tankHeight = 0.5;
    prm.autoDt = true;
    FlipSolver s(prm, ParticleScene::TankAtRest);
    stirTank(s, real(0.5));
    s.step(s.chosenDt());
    const real a0 = s.lastStats().fluidCells;

    real amin = a0, amax = a0, minPhi = 1, t = 0;
    while (t < 5) {
        const real dt = s.chosenDt();
        s.step(dt);
        t += dt;
        amin = std::min<real>(amin, s.lastStats().fluidCells);
        amax = std::max<real>(amax, s.lastStats().fluidCells);
        for (const Vec2& x : s.particles().x) minPhi = std::min(minPhi, s.solids().phi(x));
        REQUIRE(s.lastStats().divAfter <= 2 * prm.projectionTol);
    }
    MESSAGE("cuve agitee 32x32, 5 s : aire fluide dans [" << amin / a0 << ", " << amax / a0 << "] x a0");
    CHECK(amin > 0.95 * a0);
    CHECK(amax < 1.08 * a0);                         // la surface bouge : quelques cellules en plus
    CHECK(minPhi >= 0.99e-3 * s.grid().h());         // marge de repoussement respectée
    CHECK(particlesFinite(s.particles()));
}

TEST_CASE("Bulle isolee : avec p = 0, elle implose (pourquoi labelCells comble les bulles)") {
    // Une cellule d'air au fond de l'eau : p = 0 y est imposée, alors que la pression
    // hydrostatique autour vaut ρ g × profondeur. Le gradient est énorme : l'eau se
    // précipite dans la bulle, et l'effet croît quand h diminue (≈ 1/h).
    const int n = 32, depth = 24;
    const real h = real(1) / n, dt = real(1) / 60;
    auto fallSpeed = [&](bool bubble) {
        MacGrid2D g(n, n, h);
        CellLabels L(n, n);
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i) L.set(i, j, j < depth ? CellType::Fluid : CellType::Air);
        if (bubble) L.set(16, 2, CellType::Air);
        addGravity(g, Vec2{0, -kG}, dt);
        ProjectionOptions opt;
        opt.tol = 1e-10;
        project(g, L, dt, opt);
        return std::abs(g.v(16, 3));     // face entre la bulle (16, 2) et l'eau au-dessus (16, 3)
    };
    const real calm = fallSpeed(false), implosion = fallSpeed(true);
    MESSAGE("sans bulle : " << calm << " m/s ; avec bulle : " << implosion << " m/s ; g dt = " << kG * dt);
    CHECK(calm < 1e-8);
    CHECK(implosion > 3 * kG * dt);      // ≈ 0,8 m/s à 32 × 32 : près d'une demi-cellule par pas
}

TEST_CASE("Cuve avec obstacle : cellules Solid, aucune particule dedans, volume stable") {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    prm.autoDt = true;
    prm.tankObstacle = true;                 // marche [0,4 ; 0,6] × [0 ; 0,25]
    FlipSolver s(prm, ParticleScene::TankAtRest);
    stirTank(s, real(0.5));
    s.step(s.chosenDt());
    // Centres dans la marche : i = 13…18 (6 colonnes), j = 0…7 (8 rangées).
    CHECK(s.labels().count(CellType::Solid) == 48);
    const real a0 = s.lastStats().fluidCells;

    real t = 0, minPhi = 1, amin = a0;
    while (t < 3) {
        const real dt = s.chosenDt();
        s.step(dt);
        t += dt;
        amin = std::min<real>(amin, s.lastStats().fluidCells);
        for (const Vec2& x : s.particles().x) minPhi = std::min(minPhi, s.solids().phi(x));
    }
    CHECK(minPhi >= 0.99e-3 * s.grid().h());   // ni dans les murs, ni dans la marche
    CHECK(amin > 0.95 * a0);
    CHECK(particlesFinite(s.particles()));
}
