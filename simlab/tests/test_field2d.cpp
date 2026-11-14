#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "grid/Field2D.h"
#include "solvers/TestPatternSolver.h"

using namespace sim;

TEST_CASE("Field2D : indice lineaire = i + nx*j") {
    Field2D<Real> f(4, 3, Real(0.5), Stagger::Center);
    CHECK(f.size() == 12);
    CHECK(f.index(0, 0) == 0);
    CHECK(f.index(1, 0) == 1);
    CHECK(f.index(0, 1) == 4);
    CHECK(f.index(3, 2) == 11);
}

TEST_CASE("Field2D : position des echantillons selon le decalage") {
    const Real h = Real(0.5);
    CHECK(Field2D<Real>(4, 3, h, Stagger::Center).position(1, 2) == Vec2(0.75f, 1.25f));
    CHECK(Field2D<Real>(5, 3, h, Stagger::FaceX).position(1, 2) == Vec2(0.50f, 1.25f));
    CHECK(Field2D<Real>(4, 4, h, Stagger::FaceY).position(1, 2) == Vec2(0.75f, 1.00f));
    CHECK(Field2D<Real>(5, 4, h, Stagger::Node).position(1, 2) == Vec2(0.50f, 1.00f));
}

TEST_CASE("Field2D : interpolation exacte des champs affines") {
    const auto affine = [](Vec2 p) { return Real(1) + Real(2) * p.x - Real(3) * p.y; };
    const Vec2 points[] = {{0.9f, 0.6f}, {0.37f, 0.81f}, {1.5f, 0.2f}};
    for (Stagger s : {Stagger::Center, Stagger::FaceX, Stagger::FaceY, Stagger::Node}) {
        Field2D<Real> f(8, 5, Real(0.25), s);
        f.setFrom(affine);
        for (Vec2 p : points)
            CHECK(f.sample(p) == doctest::Approx(affine(p)).epsilon(1e-5));
    }
}

TEST_CASE("Field2D : valeur bornee hors du domaine") {
    Field2D<Real> f(4, 3, Real(1), Stagger::Center);
    f.setFrom([](Vec2 p) { return p.x + 10 * p.y; });
    CHECK(f.sample({-5.f, -5.f}) == doctest::Approx(f(0, 0)));
    CHECK(f.sample({50.f, 50.f}) == doctest::Approx(f(3, 2)));
}

TEST_CASE("TestPatternSolver : marqueur en haut a gauche") {
    TestPatternSolver s(64, 32, Real(1) / 32);
    const Field2D<Real>& m = *s.fields()[2].field;
    CHECK(m(16, 24) == 1);  // (0.25 Lx, 0.75 Ly)
    CHECK(m(48, 24) == 0);  // en haut à droite
    CHECK(m(16, 8) == 0);   // en bas à gauche
}
