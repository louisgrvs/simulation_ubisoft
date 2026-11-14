#include "doctest.h"
#include <random>
#include "grid/MacGrid2D.h"
#include "scenes/RotationScene.h"

TEST_CASE("MAC : tailles et positions des inconnues") {
    //La on test juste si on arrive bien à cosntruire une MACGrid comme il faut :
    const int nx = 5, ny = 3;
    const real h = 0.1;
    MacGrid2D g(nx, ny, h);

    CHECK(g.p.nx() == nx);     CHECK(g.p.ny() == ny);
    CHECK(g.u.nx() == nx + 1); CHECK(g.u.ny() == ny);
    CHECK(g.v.nx() == nx);     CHECK(g.v.ny() == ny + 1);

    // p(1,2) au centre de la cellule (1,2)
    CHECK(g.p.position(1, 2).x == doctest::Approx(1.5 * h));
    CHECK(g.p.position(1, 2).y == doctest::Approx(2.5 * h));
    // u(1,2) sur la face GAUCHE de la cellule (1,2)
    CHECK(g.u.position(1, 2).x == doctest::Approx(1.0 * h));
    CHECK(g.u.position(1, 2).y == doctest::Approx(2.5 * h));
    // v(1,2) sur la face BASSE de la cellule (1,2)
    CHECK(g.v.position(1, 2).x == doctest::Approx(1.5 * h));
    CHECK(g.v.position(1, 2).y == doctest::Approx(2.0 * h));
    // les faces extrêmes sont sur les bords du domaine
    CHECK(g.u.position(nx, 0).x == doctest::Approx(nx * h));
    CHECK(g.v.position(0, ny).y == doctest::Approx(ny * h));
}

TEST_CASE("Bilinéaire : retrouve exactement les échantillons") {
    //On vérifie que l'interpolation bilinéaire fonctionne bien :
    Field2D f(7, 4, 0.25, 0.0, 0.5);
    std::mt19937 rng(1);
    std::uniform_real_distribution<real> U(-1, 1);
    for (real& a : f.data()) a = U(rng);
    for (int j = 0; j < f.ny(); ++j)
        for (int i = 0; i < f.nx(); ++i)
            CHECK(f.sample(f.position(i, j)) == doctest::Approx(f(i, j)).epsilon(1e-12));
}

// Le test qui attrape toute erreur de décalage : l'interpolation bilinéaire est
// exacte sur les fonctions affines. Si u était lu avec le décalage de v (ou de p),
// on serait faux d'environ h/2 * pente. 
// On interpole une fonction affine pcq l'interpolation bilinéiare nous fait retrouver exactement une fonction affine
TEST_CASE("Bilinéaire : exacte sur un champ affine, pour u, v et p") {
    const int n = 16;
    const real h = 1.0 / n;
    MacGrid2D g(n, n, h);
    //Définition de la fonction affine :
    auto affine = [](Vec2 x) { return Vec2{0.3 + 2.0 * x.x - 1.5 * x.y, -0.7 + 0.5 * x.x + 3.0 * x.y}; };
    g.setVelocity(affine);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) g.p(i, j) = affine(g.p.position(i, j)).x;

    std::mt19937 rng(42);
    // Points tirés dans [h/2, 1-h/2]^2 : à l'intérieur des rectangles d'échantillons des trois champs.
    std::uniform_real_distribution<real> U(0.5 * h, 1.0 - 0.5 * h);
    for (int k = 0; k < 200; ++k) {
        Vec2 x{U(rng), U(rng)};
        Vec2 ex = affine(x);
        Vec2 vel = g.velocity(x);
        CHECK(vel.x == doctest::Approx(ex.x).epsilon(1e-12));
        CHECK(vel.y == doctest::Approx(ex.y).epsilon(1e-12));
        CHECK(g.p.sample(x) == doctest::Approx(ex.x).epsilon(1e-12));
    }
}

//Là, on vérifie que l'interpolation fait bien un principe du maximum pour une fonction convexe : min(f_{ij}) \leq f(x,y) \leq max(f_{ij})
TEST_CASE("Bilinéaire : combinaison convexe, y compris hors du domaine (clamp)") {
    Field2D f(9, 9, 0.1, 0.5, 0.5);
    std::mt19937 rng(7);
    std::uniform_real_distribution<real> U(-3, 5);
    for (real& a : f.data()) a = U(rng);
    std::uniform_real_distribution<real> X(-0.5, 1.5);  // déborde largement du domaine [0, 0.9]
    for (int k = 0; k < 500; ++k) {
        real s = f.sample({X(rng), X(rng)});
        CHECK(s >= f.minValue() - 1e-12);
        CHECK(s <= f.maxValue() + 1e-12);
    }
}

TEST_CASE("MAC : la rotation imposée est à divergence discrète nulle") {
    //La on vérifie que pour une rotation (champs affine) 
    // on a bien la condition de divergence nulle qui est gardé (le 1e-12 c'est une erreur du au flottants)
    MacGrid2D g(32, 24, 1.0 / 32);
    rotation_scene::setRotation(g, 2 * kPi);
    CHECK(g.maxAbsDivergence() < 1e-12);
}
