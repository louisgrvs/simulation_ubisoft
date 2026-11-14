// Tests du jour 4 : opérateurs MAC, gradient conjugué, projection.
// (J5 : deux retouches de compilation, aucun changement de fond :
//   - doctest est inclus comme les tests du jour 3 : "doctest.h", fourni dans tests/ ;
//   - kPi vient de core/Vec2.h : la copie locale rendait l'appel ambigu.)
#include "doctest.h"
#include <cmath>
#include <random>
#include "grid/MacGrid2D.h"
#include "numerics/MacOperators.h"
#include "numerics/ConjugateGradient.h"
#include "numerics/Projection.h"

namespace {

// Vitesse aléatoire (graine fixe) sur toutes les faces, puis u·n = 0.
void randomVelocity(MacGrid2D& g, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<real> U(-1, 1);
    for (real& a : g.u.data()) a = U(rng);
    for (real& a : g.v.data()) a = U(rng);
    enforceWallBC(g);
}

// Vitesse à divergence DISCRÈTE nulle : dérivée d'une fonction de courant ψ
// posée aux coins des cellules, nulle sur le bord (donc u·n = 0) :
//   u(i,j) =  (ψ(i,j+1) - ψ(i,j)) / h,   v(i,j) = -(ψ(i+1,j) - ψ(i,j)) / h
// D u = 0 exactement car les différences croisées se compensent (rot∘grad = 0 discret).
void divergenceFreeVelocity(MacGrid2D& g) {
    const real h = g.h(), Lx = g.size().x, Ly = g.size().y;
    auto psi = [&](int i, int j) {
        real x = i * h, y = j * h;
        return std::sin(kPi * x / Lx) * std::sin(kPi * y / Ly) * std::cos(2 * kPi * x / Lx);
    };
    for (int j = 0; j < g.u.ny(); ++j)
        for (int i = 0; i < g.u.nx(); ++i) g.u(i, j) = (psi(i, j + 1) - psi(i, j)) / h;
    for (int j = 0; j < g.v.ny(); ++j)
        for (int i = 0; i < g.v.nx(); ++i) g.v(i, j) = -(psi(i + 1, j) - psi(i, j)) / h;
}

Vector randomVector(size_t n, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<real> U(-1, 1);
    Vector x(n);
    for (real& a : x) a = U(rng);
    return x;
}

real maxDiff(const Field2D& a, const Field2D& b) {
    real m = 0;
    for (size_t k = 0; k < a.data().size(); ++k) m = std::max(m, std::abs(a.data()[k] - b.data()[k]));
    return m;
}

real faceDot(const MacGrid2D& a, const MacGrid2D& b) {
    return dotProduct(a.u.data(), b.u.data()) + dotProduct(a.v.data(), b.v.data());
}

} // namespace

// ---------------------------------------------------------------------------
//  Opérateurs discrets
// ---------------------------------------------------------------------------

TEST_CASE("Parois : u·n = 0 implique sum(D u) = 0 (théorème de la divergence discret)") {
    MacGrid2D g(40, 25, 1.0 / 40);
    randomVelocity(g, 1);
    Field2D div = g.makeCellField();
    computeDivergence(g, div);
    // somme télescopique : il ne reste que les flux de bord, tous nuls
    CHECK(std::abs(div.sum()) * g.h() < 1e-12);
}

TEST_CASE("Compatibilité MAC : D(G q) = -A q exactement") {
    const int nx = 37, ny = 23;                 // non carré, pour attraper une inversion nx/ny
    MacGrid2D g(nx, ny, 0.1);
    Field2D q = g.makeCellField();
    q.data() = randomVector(q.data().size(), 2);

    computeGradient(q, g.u, g.v);               // g.u, g.v <- G q
    Field2D dg = g.makeCellField();
    computeDivergence(g, dg);                   // D G q

    Vector Aq;
    NegLaplacianMAC{nx, ny, g.h()}(q.data(), Aq);
    real err = 0, ref = normInf(Aq);
    for (size_t k = 0; k < Aq.size(); ++k) err = std::max(err, std::abs(dg.data()[k] + Aq[k]));
    CHECK(err <= 1e-12 * ref);
}

TEST_CASE("Adjonction : <D u, q> = -<u, G q> quand u·n = 0") {
    MacGrid2D g(30, 30, 1.0 / 30), gq(30, 30, 1.0 / 30);
    randomVelocity(g, 3);
    Field2D q = g.makeCellField();
    q.data() = randomVector(q.data().size(), 4);

    Field2D div = g.makeCellField();
    computeDivergence(g, div);
    computeGradient(q, gq.u, gq.v);

    real lhs = dotProduct(div.data(), q.data());
    real rhs = -faceDot(g, gq);
    CHECK(std::abs(lhs - rhs) <= 1e-12 * std::abs(lhs));
}

TEST_CASE("Laplacien sans matrice : symétrique, positif, noyau = constantes") {
    const int nx = 20, ny = 13;
    NegLaplacianMAC A{nx, ny, 0.05};
    const size_t n = size_t(nx) * ny;
    Vector x = randomVector(n, 5), y = randomVector(n, 6), Ax, Ay;
    A(x, Ax);
    A(y, Ay);
    real xAy = dotProduct(x, Ay), Axy = dotProduct(Ax, y);
    CHECK(std::abs(xAy - Axy) <= 1e-12 * std::abs(xAy));
    CHECK(dotProduct(x, Ax) > 0);

    Vector ones(n, 1.0), A1;
    A(ones, A1);
    CHECK(normInf(A1) == 0.0);                  // exactement : n*1 - n*1
}

// ---------------------------------------------------------------------------
//  Gradient conjugué générique
// ---------------------------------------------------------------------------

TEST_CASE("CG : converge en autant d'itérations que de valeurs propres distinctes") {
    // A diagonale avec 3 valeurs propres distinctes -> 3 itérations en arithmétique exacte.
    const size_t n = 300;
    Vector d(n);
    for (size_t k = 0; k < n; ++k) d[k] = (k % 3 == 0) ? 1.0 : (k % 3 == 1 ? 10.0 : 100.0);
    LinearOperator A = [&](const Vector& x, Vector& y) {
        y.resize(x.size());
        for (size_t k = 0; k < x.size(); ++k) y[k] = d[k] * x[k];
    };
    Vector b = randomVector(n, 7), x;
    CGOptions opt;
    opt.tol = 1e-12;
    CGResult r = conjugateGradient(A, b, x, opt);
    CHECK(r.converged);
    CHECK(r.iterations <= 4);                   // 3 + une marge pour l'arrondi
    for (size_t k = 0; k < n; ++k) CHECK(x[k] == doctest::Approx(b[k] / d[k]).epsilon(1e-10));
}

TEST_CASE("CG : laplacien 1D de Dirichlet, retrouve une solution connue") {
    const int n = 64;
    LinearOperator A = [&](const Vector& x, Vector& y) {   // tridiag(-1, 2, -1)
        y.resize(x.size());
        for (int k = 0; k < n; ++k)
            y[k] = 2 * x[k] - (k > 0 ? x[k - 1] : 0) - (k < n - 1 ? x[k + 1] : 0);
    };
    Vector xStar = randomVector(n, 8), b, x;
    A(xStar, b);
    CGOptions opt;
    opt.tol = 1e-12;
    opt.maxIter = 2 * n;
    CGResult r = conjugateGradient(A, b, x, opt);
    CHECK(r.converged);
    CHECK(r.iterations <= n + 2);               // ≤ n en arithmétique exacte
    real err = 0;
    for (int k = 0; k < n; ++k) err = std::max(err, std::abs(x[k] - xStar[k]));
    CHECK(err < 1e-8);
}

// ---------------------------------------------------------------------------
//  Projection
// ---------------------------------------------------------------------------

TEST_CASE("Projection : un champ à divergence nulle est inchangé") {
    MacGrid2D g(64, 48, 1.0 / 64);
    divergenceFreeVelocity(g);
    CHECK(maxAbsDivergence(g) < 1e-10);
    MacGrid2D g0 = g;

    ProjectionStats st = project(g, 0.01);
    CHECK(st.converged);
    CHECK(st.iterations == 0);                  // b ≈ 0 : rien à résoudre
    CHECK(maxDiff(g.u, g0.u) < 1e-12);
    CHECK(maxDiff(g.v, g0.v) < 1e-12);
}

TEST_CASE("Projection : un gradient pur est annulé") {
    MacGrid2D g(64, 64, 1.0 / 64);
    Field2D phi = g.makeCellField();
    for (int j = 0; j < g.ny(); ++j)
        for (int i = 0; i < g.nx(); ++i) {
            Vec2 x = phi.position(i, j);
            phi(i, j) = std::exp(-20 * dot(x - Vec2{0.3, 0.6}, x - Vec2{0.3, 0.6}));
        }
    computeGradient(phi, g.u, g.v);             // u = G phi, nul sur les parois
    const real speed0 = g.maxSpeed();

    ProjectionOptions opt;
    opt.tol = 1e-10;
    ProjectionStats st = project(g, 0.01, opt);
    CHECK(st.converged);
    CHECK(g.maxSpeed() < 1e-6 * speed0);
}

TEST_CASE("Projection : divergence finale <= tolérance, idempotence, orthogonalité de Hodge") {
    for (int n : {16, 32, 64}) {
        CAPTURE(n);
        MacGrid2D g(n, n, 1.0 / n);
        randomVelocity(g, 9);
        MacGrid2D w = g;                        // champ avant projection

        ProjectionOptions opt;
        opt.tol = 1e-9;
        ProjectionStats st = project(g, 0.02, opt);
        CHECK(st.converged);
        CHECK(st.divBefore > 1);                // le champ aléatoire était bien divergent
        CHECK(st.divAfter <= 2 * opt.tol);      // résidu de récurrence ≈ vrai résidu

        // Orthogonalité : w = P w + G q, avec <P w, G q> = 0 (adjonction + D(Pw) ≈ 0)
        MacGrid2D gradPart = w;
        for (size_t k = 0; k < gradPart.u.data().size(); ++k) gradPart.u.data()[k] -= g.u.data()[k];
        for (size_t k = 0; k < gradPart.v.data().size(); ++k) gradPart.v.data()[k] -= g.v.data()[k];
        real ortho = faceDot(g, gradPart);
        real scale = std::sqrt(faceDot(g, g) * faceDot(gradPart, gradPart));
        CHECK(std::abs(ortho) <= 1e-8 * scale);

        // Idempotence : P(P w) = P w, sans démarrage à chaud pour un vrai test
        MacGrid2D once = g;
        g.p.data().assign(g.p.data().size(), 0);
        ProjectionOptions cold = opt;
        cold.warmStart = false;
        project(g, 0.02, cold);
        CHECK(maxDiff(g.u, once.u) < 1e-9);
        CHECK(maxDiff(g.v, once.v) < 1e-9);
    }
}

TEST_CASE("Projection : démarrage à chaud exact quand la divergence est la même") {
    // Le démarrage à chaud paie quand le champ à corriger ressemble à celui du pas
    // précédent (pression proche). Cas limite : on reprojette le MÊME champ w.
    // Attention : sur un champ DÉJÀ projeté, la bonne réponse est q = 0 et partir
    // de l'ancienne pression serait au contraire un mauvais départ.
    MacGrid2D g(48, 48, 1.0 / 48);
    randomVelocity(g, 10);
    MacGrid2D w = g;
    ProjectionStats first = project(g, 0.01);
    MacGrid2D once = g;

    w.p = g.p;                                  // même w, pression du premier passage
    ProjectionStats second = project(w, 0.01);
    CHECK(first.iterations > 10);
    CHECK(second.iterations == 0);
    CHECK(maxDiff(w.u, once.u) < 1e-12);
}
