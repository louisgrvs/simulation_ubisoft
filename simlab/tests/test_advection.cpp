#include "doctest.h"
#include <cmath>
#include <random>
#include <utility>
#include "numerics/Advection.h"
#include "scenes/RotationScene.h"

using namespace rotation_scene;

TEST_CASE("Advection : vitesse nulle = identité") {
    MacGrid2D g(16, 16, 1.0 / 16);  // u = v = 0
    Field2D q = g.makeCellField(), out = g.makeCellField();
    setSlottedDisk(q, 1.0);
    advectScalar(g, q, out, 0.1);
    for (std::size_t k = 0; k < q.data().size(); ++k) CHECK(out.data()[k] == q.data()[k]);
}

TEST_CASE("Advection : translation d'exactement une cellule = décalage exact") {
    const int n = 32;
    const real h = 1.0 / n, dt = 0.01;
    MacGrid2D g(n, n, h);
    g.setVelocity([=](Vec2) { return Vec2{h / dt, 0}; });  // on avance de h en dt

    Field2D q = g.makeCellField(), out = g.makeCellField();
    std::mt19937 rng(3);
    std::uniform_real_distribution<real> U(0, 1);
    for (real& a : q.data()) a = U(rng);

    advectScalar(g, q, out, dt);
    for (int j = 0; j < n; ++j) {
        CHECK(out(0, j) == doctest::Approx(q(0, j)));  // entrée : valeur du bord (clamp)
        for (int i = 1; i < n; ++i) CHECK(out(i, j) == doctest::Approx(q(i - 1, j)).epsilon(1e-12));
    }
}

// La rotation étant affine, la vitesse interpolée est exacte : seule l'erreur
// d'intégration en temps reste. Erreur locale : O(dt^2) pour Euler, O(dt^3) pour RK2.
TEST_CASE("Remontée des caractéristiques : ordre 1 (Euler) et ordre 2 (RK2)") {
    const int n = 64;
    const real omega = 2 * kPi;
    MacGrid2D g(n, n, 1.0 / n);
    setRotation(g, omega);
    const Vec2 c{0.5, 0.5}, x{0.72, 0.61};

    auto localError = [&](real dt, Integrator s) {
        Vec2 exact = c + rotate(x - c, -omega * dt);  // vraie position il y a dt
        return norm(traceBack(g, x, dt, s) - exact);
    };
    for (real dt : {0.01, 0.005}) {
        real orderEuler = std::log2(localError(dt, Integrator::Euler) / localError(dt / 2, Integrator::Euler));
        real orderRK2   = std::log2(localError(dt, Integrator::RK2)   / localError(dt / 2, Integrator::RK2));
        CAPTURE(dt);
        CHECK(orderEuler == doctest::Approx(2.0).epsilon(0.05));  // local O(dt^2) -> global ordre 1
        CHECK(orderRK2   == doctest::Approx(3.0).epsilon(0.05));  // local O(dt^3) -> global ordre 2
    }
}

// Rotation d'un tour complet du disque de Zalesak, à petit et à grand CFL :
//  * stabilité            : valeurs finies, quel que soit le CFL
//  * principe du maximum  : le max ne croît jamais, le min ne décroît jamais
//  * cohérence            : après un tour, le barycentre revient à sa place
TEST_CASE("Advection RK2 dans une rotation : stable, sans dépassement, barycentre conservé") {
    for (real cfl : {0.5, 2.0, 10.0}) {
        CAPTURE(cfl);
        const int n = 64;
        const real h = 1.0 / n, omega = 2 * kPi, T = 1.0;
        MacGrid2D g(n, n, h);
        setRotation(g, omega);

        const int steps = int(std::ceil(T / cflTimeStep(g, cfl)));
        const real dt = T / steps;
        CHECK(g.maxSpeed() * dt / h <= cfl + 1e-12);

        Field2D q = g.makeCellField(), tmp = g.makeCellField();
        setSlottedDisk(q, 1.0);
        const Vec2 c0 = centroid(q);
        real prevMax = q.maxValue(), prevMin = q.minValue();
        bool finite = true, maxOk = true, minOk = true;

        for (int s = 0; s < steps; ++s) {
            advectScalar(g, q, tmp, dt);
            std::swap(q, tmp);
            for (real a : q.data()) finite = finite && std::isfinite(a);
            real mx = q.maxValue(), mn = q.minValue();
            maxOk = maxOk && mx <= prevMax + 1e-12;
            minOk = minOk && mn >= prevMin - 1e-12;
            prevMax = mx;
            prevMin = mn;
        }
        CHECK(finite);
        CHECK(maxOk);
        CHECK(minOk);
        CHECK(norm(centroid(q) - c0) < 2 * h);
    }
}
