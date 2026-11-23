// Tests du jour 13 : sable (Zhu & Bridson §3.1 à 3.3).
//   1. briques de SandFriction (tenseurs, éq. 2 sur une cellule, facteur de l'éq. 4) ;
//   2. éq. 4, éq. 2 et 3 et groupes rigides (étape 5) sur une grille, sans le solveur ;
//   3. effondrement d'une colonne avec FlipSolver (scène ColumnCollapse).
#include "doctest.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include "diagnostics/SandDiagnostics.h"
#include "grid/CellLabels.h"
#include "grid/MacGrid2D.h"
#include "numerics/FreeSurface.h"
#include "numerics/SandFriction.h"
#include "solvers/FlipSolver.h"

namespace {

// Boîte n × n de pas 1/n, toutes les cellules Fluid (hors domaine = Solid).
struct Box {
    MacGrid2D g;
    CellLabels L;
    explicit Box(int n) : g(n, n, real(1) / n), L(n, n) { L.fill(CellType::Fluid); }
};

// Nombre dans [0, 1[ tiré directement de std::mt19937 (CONVENTIONS.md : pas de
// std::uniform_real_distribution, dont l'algorithme dépend du compilateur).
real uniform01(std::mt19937& rng) { return real(rng()) / (real(std::mt19937::max()) + 1); }

real maxAbsDiff(const Field2D& a, const Field2D& b) {
    real m = 0;
    for (std::size_t k = 0; k < a.data().size(); ++k)
        m = std::max(m, std::abs(a.data()[k] - b.data()[k]));
    return m;
}

// Remplit u et v directement aux positions des faces (CONVENTIONS.md) :
// u(i, j) en (i h, (j+1/2) h), v(i, j) en ((i+1/2) h, j h).
template <class F>
void setFaces(MacGrid2D& g, F f) {
    const real h = g.h();
    for (int j = 0; j < g.ny(); ++j)
        for (int i = 0; i <= g.nx(); ++i) g.u(i, j) = f(Vec2{ i * h, (j + real(0.5)) * h }).x;
    for (int j = 0; j <= g.ny(); ++j)
        for (int i = 0; i < g.nx(); ++i) g.v(i, j) = f(Vec2{ (i + real(0.5)) * h, j * h }).y;
}

} // namespace

// ============================================================================
//  1. Briques
// ============================================================================
TEST_CASE("Sable, briques : norme de Frobenius et déviateur 2D") {
    const Sym2 S{ 1, 2, 3 };
    CHECK(frobeniusNorm(S) == doctest::Approx(std::sqrt(real(18))));   // 1 + 2·4 + 9
    const Sym2 d = deviator(S);
    CHECK(d.xx == doctest::Approx(-1));
    CHECK(d.xy == doctest::Approx(2));
    CHECK(d.yy == doctest::Approx(1));
    CHECK(d.xx + d.yy == doctest::Approx(0));                        // trace nulle en 2D : tr/2, pas tr/3
}

TEST_CASE("Sable : coefficient du seuil en 2D, k = √2 sin φ") {
    CHECK(frictionCoefficient2D(30) == doctest::Approx(std::sqrt(real(2)) / 2));
    CHECK(frictionCoefficient2D(90) == doctest::Approx(std::sqrt(real(2))));
    CHECK(frictionCoefficient2D(0) == doctest::Approx(0));
}

TEST_CASE("Sable, éq. 2 : cisaillement simple ⇒ σ_xy = sin φ · p, dans le sens de D") {
    const real p = 2, k = frictionCoefficient2D(30);   // sin 30° · 2 = 1
    for (real gamma : { real(3), real(300) }) {         // le frottement ne dépend pas de la vitesse
        const Sym2 D{ 0, gamma / 2, 0 };
        const Sym2 s = frictionalStressAtCell(D, p, k, /*etaMax*/ 0, real(1e-8));
        CHECK(s.xy == doctest::Approx(1));               // > 0 : même signe que D, donc dissipatif
        CHECK(s.xx == doctest::Approx(0));
        CHECK(s.yy == doctest::Approx(0));
    }
}

TEST_CASE("Sable, éq. 2 : σ est déviatorique et posé sur le seuil |σ|_F = k p") {
    const real p = real(1.7), k = frictionCoefficient2D(35);
    const Sym2 D{ real(1), real(0.3), real(0.2) };   // trace non nulle : seul le déviateur compte
    const Sym2 s = frictionalStressAtCell(D, p, k, 0, real(1e-8));
    CHECK(s.xx + s.yy == doctest::Approx(0).epsilon(1e-12));
    CHECK(frobeniusNorm(s) == doctest::Approx(k * p));
}

TEST_CASE("Sable, éq. 2 : |D| ≈ 0 ou p = 0 ⇒ σ = 0, sans NaN") {
    const real k = frictionCoefficient2D(30);
    bool clamped = true;
    const Sym2 s0 = frictionalStressAtCell(Sym2{}, 1, k, 0, real(1e-8), &clamped);
    CHECK(std::isfinite(s0.xx));
    CHECK(frobeniusNorm(s0) == 0);
    CHECK_FALSE(clamped);
    const Sym2 s1 = frictionalStressAtCell(Sym2{ 0, 1, 0 }, 0, k, 0, real(1e-8));
    CHECK(frobeniusNorm(s1) == 0);
}

TEST_CASE("Sable, éq. 3 : sous le seuil σ = σ_rigid, au-dessus σ = σ_f") {
    const real k = frictionCoefficient2D(30), p = 1;
    bool rigid = false;
    // |D| petit : σ_rigid = 0,1 · dev D est bien sous le seuil k p ≈ 0,71.
    const Sym2 Dsmall{ 0, real(1e-3), 0 };
    const Sym2 s = frictionalStressAtCell(Dsmall, p, k, real(0.1), real(1e-8), &rigid);
    CHECK(rigid);
    CHECK(s.xy == doctest::Approx(real(0.1) * Dsmall.xy));
    // |D| grand : σ_rigid dépasserait le seuil, on prend σ_f, posé sur le seuil.
    const Sym2 Dbig{ 0, real(50), 0 };
    const Sym2 s2 = frictionalStressAtCell(Dbig, p, k, real(0.1), real(1e-8), &rigid);
    CHECK_FALSE(rigid);
    CHECK(frobeniusNorm(s2) == doctest::Approx(k * p));
    // rigidScale ≤ 0 : éq. 3 désactivée, toujours σ_f.
    const Sym2 s3 = frictionalStressAtCell(Dsmall, p, k, 0, real(1e-8), &rigid);
    CHECK_FALSE(rigid);
    CHECK(frobeniusNorm(s3) == doctest::Approx(k * p));
}

TEST_CASE("Sable : taux de déformation exact pour un champ affine, nul pour un mouvement rigide") {
    Box b(10);
    const real a = real(0.7), bb = real(-1.3), c = real(0.4);
    setFaces(b.g, [=](Vec2 x) { return Vec2{ a * x.x + bb * x.y, c * x.x - a * x.y }; });
    for (int j = 1; j < 9; ++j)
        for (int i = 1; i < 9; ++i) {
            const Sym2 D = strainRateAtCell(b.g, i, j);
            CHECK(D.xx == doctest::Approx(a));
            CHECK(D.yy == doctest::Approx(-a));
            CHECK(D.xy == doctest::Approx((bb + c) / 2));
        }
    // Translation + rotation : D = 0 partout, y compris dans les cellules du bord.
    const real w = real(2.5);
    setFaces(b.g, [=](Vec2 x) { return Vec2{ real(0.3) - w * (x.y - real(0.5)), real(-0.2) + w * (x.x - real(0.5)) }; });
    for (int j = 0; j < 10; ++j)
        for (int i = 0; i < 10; ++i)
            CHECK(frobeniusNorm(strainRateAtCell(b.g, i, j)) == doctest::Approx(0).epsilon(1e-12));
}

TEST_CASE("Sable, éq. 4 : facteur de Coulomb max(0, 1 − μ|u_N|/|u_T|)") {
    CHECK(wallFrictionFactor(1, real(-0.5), real(0.6)) == doctest::Approx(0.7));
    CHECK(wallFrictionFactor(-1, real(-0.5), real(0.6)) == doctest::Approx(0.7));   // symétrique en u_T
    CHECK(wallFrictionFactor(real(0.2), real(-0.5), real(0.6)) == 0);               // arrêt, sans inversion
    CHECK(wallFrictionFactor(1, real(0.5), real(0.6)) == 1);                        // s'éloigne du mur
    CHECK(wallFrictionFactor(0, real(-0.5), real(0.6)) == 1);                       // rien à freiner, pas de NaN
    CHECK(wallFrictionFactor(1, real(-0.5), 0) == 1);                               // μ = 0 : glissement parfait
}

// ============================================================================
//  2. Sur une grille
// ============================================================================
TEST_CASE("Sable, éq. 4 sur la grille : le sol freine la rangée du bas, et elle seule") {
    const int n = 8;
    Box b(n);
    for (int j = 0; j < n; ++j)
        for (int i = 1; i < n; ++i) b.g.u(i, j) = 1;   // glissement vers la droite
    for (int i = 0; i < n; ++i) b.g.v(i, 0) = real(-0.5);  // poussé dans le sol (comme −g dt)
    const Field2D v0 = b.g.v;
    const WallFrictionStats st = applyWallFriction(b.g, b.L, real(0.6));
    for (int i = 1; i < n; ++i) {
        CHECK(b.g.u(i, 0) == doctest::Approx(0.7));
        CHECK(b.g.u(i, 1) == doctest::Approx(1));          // pas de mur : intact
        CHECK(b.g.u(i, n - 1) == doctest::Approx(1));      // plafond : u·n = 0, pas de frottement
    }
    CHECK(maxAbsDiff(b.g.v, v0) == 0);                     // l'éq. 4 ne touche pas aux normales
    CHECK(st.frictionFaces == n - 1);

    // Même chose, mais le fluide quitte le sol : aucun frottement.
    for (int i = 0; i < n; ++i) b.g.v(i, 0) = real(0.5);
    for (int i = 1; i < n; ++i) b.g.u(i, 0) = 1;
    applyWallFriction(b.g, b.L, real(0.6));
    for (int i = 1; i < n; ++i) CHECK(b.g.u(i, 0) == doctest::Approx(1));
}

TEST_CASE("Sable, éq. 4 sur la grille : mur de gauche (rôles de x et y échangés)") {
    const int n = 8;
    Box b(n);
    for (int j = 1; j < n; ++j)
        for (int i = 0; i < n; ++i) b.g.v(i, j) = 1;   // glissement vers le haut
    for (int j = 0; j < n; ++j) b.g.u(0, j) = real(-0.25);  // poussé dans le mur de gauche
    applyWallFriction(b.g, b.L, real(2));               // μ|u_N| = 0,5
    for (int j = 1; j < n; ++j) {
        CHECK(b.g.v(0, j) == doctest::Approx(0.5));
        CHECK(b.g.v(1, j) == doctest::Approx(1));
    }
}

TEST_CASE("Sable, éq. 2 sur la grille : un mouvement rigide n'est pas freiné") {
    Box b(12);
    b.g.p.fill(5);
    const real w = real(1.5);
    setFaces(b.g, [=](Vec2 x) { return Vec2{ real(0.4) - w * (x.y - real(0.5)), w * (x.x - real(0.5)) }; });
    const Field2D u0 = b.g.u, v0 = b.g.v;
    const FrictionalStressStats st = applyFrictionalStress(b.g, b.L, real(1) / 60);
    CHECK(st.stressCells == 0);
    CHECK(maxAbsDiff(b.g.u, u0) < 1e-12);
    CHECK(maxAbsDiff(b.g.v, v0) < 1e-12);
}

TEST_CASE("Sable, éq. 2 sur la grille : sans pression, aucun frottement") {
    Box b(12);
    b.g.p.fill(0);
    setFaces(b.g, [](Vec2 x) { return Vec2{ 3 * (x.y - real(0.5)), 0 }; });
    const Field2D u0 = b.g.u;
    applyFrictionalStress(b.g, b.L, real(1) / 60);
    CHECK(maxAbsDiff(b.g.u, u0) == 0);
}

TEST_CASE("Sable, éq. 2 sur la grille : pression négative traitée comme p = 0") {
    Box b(12);
    b.g.p.fill(-3);                                       // pas de traction dans le sable
    setFaces(b.g, [](Vec2 x) { return Vec2{ 3 * (x.y - real(0.5)), 0 }; });
    const Field2D u0 = b.g.u;
    applyFrictionalStress(b.g, b.L, real(1) / 60);
    CHECK(maxAbsDiff(b.g.u, u0) == 0);
}

TEST_CASE("Sable, éq. 2 sur la grille : cisaillement simple, freiné au bord, intact à l'intérieur") {
    // u = γ (y − 1/2) : D uniforme. Le frottement ne dépend pas de |D|, donc σ est le même
    // dans toutes les cellules, et sa divergence n'est non nulle qu'au bord.
    const int n = 16;
    Box b(n);
    b.g.p.fill(1);
    const real gamma = 2;
    setFaces(b.g, [=](Vec2 x) { return Vec2{ gamma * (x.y - real(0.5)), 0 }; });
    enforceSolidBC(b.g, b.L);
    const Field2D u0 = b.g.u;
    const FrictionalStressStats st = applyFrictionalStress(b.g, b.L, real(1) / 600);
    CHECK(st.rigidCells == 0);                            // σ_rigid > k p partout : on est sur le seuil
    const int i = n / 2;
    CHECK(b.g.u(i, n - 1) < u0(i, n - 1));               // rangée du haut (u > 0) : ralentie
    CHECK(b.g.u(i, 0) > u0(i, 0));                       // rangée du bas (u < 0) : ralentie aussi
    CHECK(b.g.u(i, n / 2) == doctest::Approx(u0(i, n / 2)).epsilon(1e-12));   // intérieur : intact
}

TEST_CASE("Sable, éq. 2 sur la grille : les faces dans l'air ne bougent pas") {
    const int n = 12;
    Box b(n);
    for (int j = n / 2; j < n; ++j)
        for (int i = 0; i < n; ++i) b.L.set(i, j, CellType::Air);   // moitié haute : air
    b.g.p.fill(1);
    setFaces(b.g, [](Vec2 x) { return Vec2{ 3 * x.y, -x.x }; });
    const Field2D u0 = b.g.u, v0 = b.g.v;
    applyFrictionalStress(b.g, b.L, real(1) / 60);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i <= n; ++i)
            if (!uFaceIsFluid(b.L, i, j)) CHECK(b.g.u(i, j) == u0(i, j));
    for (int j = 0; j <= n; ++j)
        for (int i = 0; i < n; ++i)
            if (!vFaceIsFluid(b.L, i, j)) CHECK(b.g.v(i, j) == v0(i, j));
}

TEST_CASE("Sable, éq. 2 et 3 sur la grille : le frottement dissipe l'énergie") {
    const int n = 16;
    SUBCASE("champ aléatoire, p ~ 10") {
        std::mt19937 rng(7);
        Box b(n);
        for (real& u : b.g.u.data()) u = 2 * uniform01(rng) - 1;
        for (real& v : b.g.v.data()) v = 2 * uniform01(rng) - 1;
        for (real& p : b.g.p.data()) p = 10 * uniform01(rng);
        enforceSolidBC(b.g, b.L);
        const FrictionalStressStats st = applyFrictionalStress(b.g, b.L, real(1) / 60);
        MESSAGE("aléatoire : E après / E avant = " << st.energyAfter / st.energyBefore
                << ", cellules rigides " << st.rigidCells << " / " << st.stressCells);
        CHECK(st.energyBefore > 0);
        CHECK(st.energyAfter < st.energyBefore);
        CHECK(st.stressCells == n * n);
        CHECK(st.rigidCells > 0);
        CHECK(st.rigidCells < n * n);                      // les deux régimes sont présents
    }
    SUBCASE("tourbillon lisse, p énorme : tout est sous le seuil, σ_rigid évite le dépassement") {
        Box b(n);
        setFaces(b.g, [](Vec2 x) {
            return Vec2{ std::sin(kPi * x.x) * std::cos(kPi * x.y), -std::cos(kPi * x.x) * std::sin(kPi * x.y) };
        });
        b.g.p.fill(real(1e6));
        enforceSolidBC(b.g, b.L);
        const FrictionalStressStats st = applyFrictionalStress(b.g, b.L, real(1) / 60);
        MESSAGE("tourbillon : E après / E avant = " << st.energyAfter / st.energyBefore);
        CHECK(st.rigidCells == n * n);
        CHECK(st.energyAfter < st.energyBefore);
    }
    SUBCASE("sans l'éq. 3, σ_f seule dépasse et CRÉE de l'énergie (ce que l'éq. 3 corrige)") {
        Box b(n);
        setFaces(b.g, [](Vec2 x) {
            return Vec2{ std::sin(kPi * x.x) * std::cos(kPi * x.y), -std::cos(kPi * x.x) * std::sin(kPi * x.y) };
        });
        b.g.p.fill(100);
        enforceSolidBC(b.g, b.L);
        FrictionalStressOptions opt;
        opt.rigidStress = false;
        const FrictionalStressStats st = applyFrictionalStress(b.g, b.L, real(1) / 60, opt);
        MESSAGE("σ_f seule : E après / E avant = " << st.energyAfter / st.energyBefore);
        CHECK(st.energyAfter > st.energyBefore);
    }
}

TEST_CASE("Sable, étape 5 : composantes connexes en 4-voisinage") {
    // 5 × 4 :   . X X . X        (j = 3, en haut)
    //           . X . . X
    //           . . . X .        X(3,1) ne touche X(4,2) qu'en diagonale : groupe à part
    //           X . . . .        (j = 0)
    const int nx = 5, ny = 4;
    std::vector<char> m(nx * ny, 0);
    auto set = [&](int i, int j) { m[i + nx * j] = 1; };
    set(1, 3); set(2, 3); set(1, 2);   // groupe A
    set(4, 3); set(4, 2);              // groupe B
    set(3, 1);                         // groupe C (diagonale de B seulement)
    set(0, 0);                         // groupe D
    std::vector<int> grp;
    CHECK(labelRigidGroups(m, nx, ny, grp) == 4);
    CHECK(grp[1 + nx * 3] == grp[2 + nx * 3]);
    CHECK(grp[1 + nx * 3] == grp[1 + nx * 2]);
    CHECK(grp[4 + nx * 3] == grp[4 + nx * 2]);
    CHECK(grp[3 + nx * 1] != grp[4 + nx * 2]);
    CHECK(grp[0] >= 0);
    CHECK(grp[2 + nx * 2] == -1);
}

TEST_CASE("Sable, étape 5 : un mouvement rigide est retrouvé exactement") {
    const int n = 10;
    Box b(n);
    const Vec2 U{ real(0.3), real(-0.7) };
    const real w = real(1.7);
    setFaces(b.g, [=](Vec2 x) { return Vec2{ U.x - w * (x.y - real(0.5)), U.y + w * (x.x - real(0.5)) }; });
    const Field2D u0 = b.g.u, v0 = b.g.v;
    std::vector<int> grp(n * n, -1);
    for (int j = 2; j < 7; ++j)
        for (int i = 3; i < 8; ++i) grp[i + n * j] = 0;   // un bloc 5 × 5
    int faces = 0;
    const std::vector<RigidMotion2D> m = projectRigidGroups(b.g, b.L, grp, 1, &faces);
    CHECK(faces == 6 * 5 * 2);                          // 6 faces u par rangée × 5 rangées, idem en v
    CHECK(m[0].omega == doctest::Approx(w));
    CHECK(maxAbsDiff(b.g.u, u0) < 1e-12);
    CHECK(maxAbsDiff(b.g.v, v0) < 1e-12);
}

TEST_CASE("Sable, étape 5 : la projection conserve P et L du groupe, et ne crée pas d'énergie") {
    const int n = 12;
    Box b(n);
    std::mt19937 rng(3);
    for (real& u : b.g.u.data()) u = 2 * uniform01(rng) - 1;
    for (real& v : b.g.v.data()) v = 2 * uniform01(rng) - 1;
    enforceSolidBC(b.g, b.L);
    std::vector<int> grp(n * n, -1);
    for (int j = 1; j < 6; ++j)
        for (int i = 2; i < 9; ++i) grp[i + n * j] = 0;   // groupe 0 : bloc 7 × 5
    grp[10 + n * 9] = 1;                                   // groupe 1 : une seule cellule
    const real h = b.g.h();
    // Quantités du groupe 0, sur ses faces (celles qui touchent une de ses cellules).
    auto inGroup0U = [&](int i, int j) {
        return (i - 1 >= 0 && grp[(i - 1) + n * j] == 0) || (i < n && grp[i + n * j] == 0);
    };
    auto inGroup0V = [&](int i, int j) {
        return (j - 1 >= 0 && grp[i + n * (j - 1)] == 0) || (j < n && grp[i + n * j] == 0);
    };
    const Vec2 c{ real(5.5) * h, real(3.5) * h };          // centre des cellules du bloc
    auto measure = [&](real& Px, real& Py, real& Lz, real& E) {
        Px = Py = Lz = E = 0;
        for (int j = 0; j < n; ++j)
            for (int i = 0; i <= n; ++i)
                if (uFaceIsFluid(b.L, i, j) && inGroup0U(i, j)) {
                    const real u = b.g.u(i, j);
                    Px += u; Lz -= ((j + real(0.5)) * h - c.y) * u; E += u * u;
                }
        for (int j = 0; j <= n; ++j)
            for (int i = 0; i < n; ++i)
                if (vFaceIsFluid(b.L, i, j) && inGroup0V(i, j)) {
                    const real v = b.g.v(i, j);
                    Py += v; Lz += ((i + real(0.5)) * h - c.x) * v; E += v * v;
                }
    };
    real Px0, Py0, L0, E0, Px1, Py1, L1, E1;
    measure(Px0, Py0, L0, E0);
    const std::vector<RigidMotion2D> m = projectRigidGroups(b.g, b.L, grp, 2);
    measure(Px1, Py1, L1, E1);
    CHECK(Px1 == doctest::Approx(Px0).epsilon(1e-12));
    CHECK(Py1 == doctest::Approx(Py0).epsilon(1e-12));
    CHECK(L1 == doctest::Approx(L0).epsilon(1e-12));
    CHECK(E1 <= E0);
    // Groupe d'une cellule : rotation non mesurable, ω = 0, U = moyenne des faces.
    CHECK(m[1].omega == 0);
    CHECK(b.g.u(10, 9) == doctest::Approx(b.g.u(11, 9)));
}

TEST_CASE("Sable, étape 5 sur la grille : sans groupes, rien ne change par rapport à l'éq. 3 seule") {
    const int n = 16;
    std::mt19937 rng(11);
    Box a(n);
    for (real& u : a.g.u.data()) u = 2 * uniform01(rng) - 1;
    for (real& v : a.g.v.data()) v = 2 * uniform01(rng) - 1;
    for (real& p : a.g.p.data()) p = 10 * uniform01(rng);
    enforceSolidBC(a.g, a.L);
    Box b(n);
    b.g = a.g;
    FrictionalStressOptions noGroups;
    noGroups.rigidGroups = false;
    const FrictionalStressStats sa = applyFrictionalStress(a.g, a.L, real(1) / 60, noGroups);
    const FrictionalStressStats sb = applyFrictionalStress(b.g, b.L, real(1) / 60);
    CHECK(sa.rigidGroups == 0);
    CHECK(sb.rigidGroups > 0);
    CHECK(sb.rigidCells == sa.rigidCells);              // même classification (éq. 3)
    CHECK(sb.energyAfter <= sa.energyAfter);            // les groupes dissipent davantage
}

// ============================================================================
//  3. Effondrement d'une colonne (FlipSolver, scène ColumnCollapse)
// ============================================================================
namespace {

FlipParams columnParams(Material m, real phiDeg) {
    FlipParams prm;
    prm.nx = prm.ny = 48;
    prm.scheme = TransferScheme::PIC;   // PIC pour les deux matériaux : comparaison à schéma égal
    prm.dt = real(1) / 120;
    prm.material = m;
    prm.frictionAngleDeg = phiDeg;
    prm.columnWidth = real(0.2);
    prm.columnHeight = real(0.4);
    return prm;
}

PileShape runColumn(Material m, real phiDeg, real tEnd, real width = real(0.2), real height = real(0.4)) {
    FlipParams prm = columnParams(m, phiDeg);
    prm.columnWidth = width;
    prm.columnHeight = height;
    FlipSolver s(prm, ParticleScene::ColumnCollapse);
    const int steps = int(std::lround(tEnd / prm.dt));
    for (int k = 0; k < steps; ++k) s.step(prm.dt);
    return measurePile(s.particles(), s.grid().nx(), s.grid().h());
}

} // namespace

TEST_CASE("Colonne : la scène est posée au centre du sol, au repos") {
    FlipSolver s(columnParams(Material::Sand, 30), ParticleScene::ColumnCollapse);
    REQUIRE(s.particles().size() > 0);
    CHECK(s.usesGravity());
    CHECK(s.isSand());
    const Vec2 G = s.centerOfMass();
    const real tol = s.grid().h() / 2;                   // une sous-cellule est gardée si son CENTRE est dedans
    CHECK(G.x == doctest::Approx(0.5).epsilon(0.01));
    CHECK(G.y == doctest::Approx(0.2).epsilon(0.02));    // milieu de la hauteur 0,4
    for (const Vec2& x : s.particles().x) {
        CHECK(x.x > real(0.4) - tol);
        CHECK(x.x < real(0.6) + tol);
        CHECK(x.y < real(0.4) + tol);
    }
}

TEST_CASE("Colonne : en eau, les étapes du sable ne s'exécutent pas") {
    FlipSolver s(columnParams(Material::Water, 30), ParticleScene::ColumnCollapse);
    for (int k = 0; k < 5; ++k) s.step(real(1) / 120);
    CHECK(s.lastStats().frictionFaces == 0);
    CHECK(s.lastStats().stressCells == 0);
}

TEST_CASE("Colonne : en sable, la friction ne crée jamais d'énergie") {
    FlipParams prm = columnParams(Material::Sand, 30);
    FlipSolver s(prm, ParticleScene::ColumnCollapse);
    int active = 0;
    for (int k = 0; k < 120; ++k) {
        s.step(prm.dt);
        const FlipStepStats& st = s.lastStats();
        if (st.stressCells > 0) ++active;
        CHECK(st.energyAfterFriction <= st.energyBeforeFriction * (1 + 1e-12));
    }
    CHECK(active > 100);                                  // l'éq. 2 agit bien à presque tous les pas
}

TEST_CASE("Colonne : au début de l'effondrement, le sable s'étale moins vite que l'eau") {
    const PileShape water = runColumn(Material::Water, 30, real(0.25));
    const PileShape sand = runColumn(Material::Sand, 30, real(0.25));
    MESSAGE("t = 0,25 s, eau  : largeur " << water.spread << ", hauteur " << water.height);
    MESSAGE("t = 0,25 s, sable: largeur " << sand.spread << ", hauteur " << sand.height);
    CHECK(sand.spread < water.spread - real(0.1));
}

TEST_CASE("Colonne : le sable forme un tas, d'autant plus étroit et pentu que φ est grand") {
    // Colonne 0,3 × 0,3 (celle de sand_demo) : moins haute, elle finit plus près du repos.
    const real W = real(0.3), H = real(0.3);
    const PileShape water = runColumn(Material::Water, 30, 2.0, W, H);
    const real phis[] = { 20, 30, 45 };
    PileShape sand[3];
    for (int k = 0; k < 3; ++k) {
        sand[k] = runColumn(Material::Sand, phis[k], 2.0, W, H);
        MESSAGE("φ = " << phis[k] << "° : largeur " << sand[k].spread << ", hauteur " << sand[k].height
                       << ", angle mesuré " << sand[k].angleDeg << "°");
    }
    CHECK(water.spread > real(0.9));                     // l'eau touche les deux murs
    CHECK(sand[2].spread < real(0.75));                  // φ = 45° : un tas, loin des murs
    CHECK(sand[0].spread > sand[1].spread);
    CHECK(sand[1].spread > sand[2].spread);
    CHECK(sand[0].angleDeg < sand[1].angleDeg);
    CHECK(sand[1].angleDeg < sand[2].angleDeg);
    CHECK(std::abs(sand[0].angleDeg - 20) < 8);          // angle du tas ≈ φ (pour φ modéré)
    CHECK(std::abs(sand[1].angleDeg - 30) < 8);
}

TEST_CASE("Colonne : le sable garde son volume grâce à la deuxième projection") {
    // Particules par cellule fluide : 4 au départ. Si le sable se tasse, ce nombre monte.
    auto particlesPerCell = [](bool reproject) {
        FlipParams prm = columnParams(Material::Sand, 45);
        prm.columnWidth = prm.columnHeight = real(0.3);
        prm.reprojectAfterFriction = reproject;
        FlipSolver s(prm, ParticleScene::ColumnCollapse);
        for (int k = 0; k < 240; ++k) s.step(prm.dt);
        return double(s.particles().size()) / s.lastStats().fluidCells;
    };
    const double with = particlesPerCell(true), without = particlesPerCell(false);
    MESSAGE("particules par cellule fluide à t = 2 s : " << with << " avec, " << without << " sans");
    CHECK(with < 4.5);
    CHECK(without > 6);
}

TEST_CASE("Colonne : sans groupes rigides, le sable s'aplatit (ce que l'étape 5 corrige)") {
    FlipParams prm = columnParams(Material::Sand, 30);
    prm.rigidGroups = false;
    FlipSolver s(prm, ParticleScene::ColumnCollapse);
    for (int k = 0; k < 240; ++k) s.step(prm.dt);
    const PileShape r = measurePile(s.particles(), s.grid().nx(), s.grid().h());
    MESSAGE("sans groupes : largeur " << r.spread << ", hauteur " << r.height);
    CHECK(r.spread > real(0.9));
}
