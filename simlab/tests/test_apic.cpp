// Tests du jour 10 : APIC sur grille MAC (Jiang et al. 2015, #5.1, #5.3, #6).
//
//   1. pochoir gradient : identités qui rendent APIC exact pour les champs affines ;
//   2. transferts seuls : champ affine, rotation rigide, quantité de mouvement,
//      pochoir tronqué aux bords, régressions PIC ;
//   3. solveur complet : disque en rotation, repos, tourbillon, rupture de barrage.
//   4. jour 11 : disque en rotation libre (pas complet) et diagnostics de forme.
#include "doctest.h"
#include <algorithm>
#include <cmath>
#include <random>
#include "core/Kernels.h"
#include "diagnostics/ParticleDiagnostics.h"
#include "numerics/MacTransfer.h"
#include "scenes/Seeding.h"
#include "solvers/FlipSolver.h"

namespace {

// Disque de particules centré, loin des parois : tous les pochoirs sont complets.
ParticleSet2D interiorDiskApic(int n, unsigned seed = 7) {
    ParticleSet2D p;
    SeedOptions opt;
    opt.perAxis = 2;
    opt.jitter = 1;
    opt.seed = seed;
    seedParticles(p, n, n, real(1) / n, region::disk({0.5, 0.5}, 0.3), opt);
    return p;
}

// Champ affine quelconque u(x) = A x + b. A n'est PAS symétrique : une
// transposition de C (colonnes au lieu des lignes) serait détectée.
const real kA[2][2] = {{0.3, -1.1}, {0.8, 0.4}};
const Vec2 kB{0.2, -0.5};
Vec2 affineField(Vec2 x) {
    return {kA[0][0] * x.x + kA[0][1] * x.y + kB.x, kA[1][0] * x.x + kA[1][1] * x.y + kB.y};
}
const Vec2 kCx{kA[0][0], kA[0][1]};   // ligne 0 de A : (du/dx, du/dy)
const Vec2 kCy{kA[1][0], kA[1][1]};   // ligne 1 de A : (dv/dx, dv/dy)

// Remplit la grille avec le champ affine, évalué à la position de chaque face.
void fillAffine(MacGrid2D& g) {
    for (int j = 0; j < g.u.ny(); ++j)
        for (int i = 0; i < g.u.nx(); ++i) g.u(i, j) = affineField(g.u.position(i, j)).x;
    for (int j = 0; j < g.v.ny(); ++j)
        for (int i = 0; i < g.v.nx(); ++i) g.v(i, j) = affineField(g.v.position(i, j)).y;
}

bool allFiniteParticles(const ParticleSet2D& p) {
    for (std::size_t k = 0; k < p.size(); ++k)
        for (real r : {p.x[k].x, p.x[k].y, p.v[k].x, p.v[k].y, p.cx[k].x, p.cx[k].y, p.cy[k].x, p.cy[k].y})
            if (!std::isfinite(r)) return false;
    return true;
}

real maxDiff(const std::vector<Vec2>& a, const std::vector<Vec2>& b) {
    real e = 0;
    for (std::size_t k = 0; k < a.size(); ++k) e = std::max(e, norm(a[k] - b[k]));
    return e;
}

} // namespace

// ---------------------------------------------------------------------------
//  1. Pochoir gradient
// ---------------------------------------------------------------------------

TEST_CASE("APIC, pochoir : somme des gradients nulle, somme des gradients ⊗ x_f = identité") {
    // Ce sont les deux identités qui rendent G2P exact pour un champ affine :
    // \sum_f \nabla w_f = 0 (une constante a une pente nulle) et
    // \sum_f \nabla w_f x_f^T = I (la pente de x est 1).
    MacGrid2D g(16, 10, 1.0 / 16);
    std::mt19937 rng(3);
    std::uniform_real_distribution<real> X(0.1, 0.9), Y(0.1, 0.52);
    for (int t = 0; t < 200; ++t) {
        const Vec2 x{X(rng), Y(rng)};
        for (const Field2D* f : {&g.u, &g.v}) {
            const Stencil2D s = linearStencil(*f, x);
            Vec2 sumGrad{0, 0};
            real I[2][2] = {{0, 0}, {0, 0}};   // \sum_f \nabla w_f x_f^T
            for (int b = 0; b < 2; ++b)
                for (int a = 0; a < 2; ++a) {
                    const Vec2 gw{s.dwx[a] * s.wy[b], s.wx[a] * s.dwy[b]};
                    const Vec2 xf = f->position(s.i0 + a, s.j0 + b);
                    sumGrad += gw;
                    I[0][0] += gw.x * xf.x; I[0][1] += gw.x * xf.y;
                    I[1][0] += gw.y * xf.x; I[1][1] += gw.y * xf.y;
                }
            CHECK(std::abs(sumGrad.x) < 1e-12);
            CHECK(std::abs(sumGrad.y) < 1e-12);
            CHECK(I[0][0] == doctest::Approx(1.0).epsilon(1e-12));
            CHECK(I[1][1] == doctest::Approx(1.0).epsilon(1e-12));
            CHECK(std::abs(I[0][1]) < 1e-12);
            CHECK(std::abs(I[1][0]) < 1e-12);
        }
    }
}

TEST_CASE("APIC, pochoir : dwx, dwy sont les dérivées des poids (différences finies)") {
    MacGrid2D g(16, 16, 1.0 / 16);
    const real h = g.h(), eps = 1e-7 * h;
    std::mt19937 rng(4);
    std::uniform_real_distribution<real> X(0.1, 0.9);
    for (int t = 0; t < 100; ++t) {
        const Vec2 x{X(rng), X(rng)};
        for (const Field2D* f : {&g.u, &g.v}) {
            const Stencil2D s = linearStencil(*f, x);
            // On reste dans la même cellule (sinon le pochoir change) : on saute les points trop près d'un bord.
            const Stencil2D sx = linearStencil(*f, x + Vec2{eps, 0});
            const Stencil2D sy = linearStencil(*f, x + Vec2{0, eps});
            if (sx.i0 != s.i0 || sy.j0 != s.j0) continue;
            for (int b = 0; b < 2; ++b)
                for (int a = 0; a < 2; ++a) {
                    const real w = s.wx[a] * s.wy[b];
                    const real dwdx = (sx.wx[a] * sx.wy[b] - w) / eps;
                    const real dwdy = (sy.wx[a] * sy.wy[b] - w) / eps;
                    CHECK(s.dwx[a] * s.wy[b] == doctest::Approx(dwdx).epsilon(1e-6));
                    CHECK(s.wx[a] * s.dwy[b] == doctest::Approx(dwdy).epsilon(1e-6));
                }
        }
    }
}

// ---------------------------------------------------------------------------
//  2. Transferts seuls
// ---------------------------------------------------------------------------

TEST_CASE("APIC : l'aller-retour d'un champ affine est exact (grille, vitesses et c)") {
    // Pour u = A x + b et c tiré des lignes de A :
    //   P2G : chaque particule apporte u(x_p) + A (x_f - x_p) = u(x_f), donc u_f = u(x_f) ;
    //   G2P : le bilinéaire reproduit les fonctions affines, donc u_p et c reviennent.
    // PIC, lui, fait une erreur A · \sum w (x̄_f - x_f) (test du jour 7).
    const int n = 32;
    MacGrid2D g(n, n, 1.0 / n);
    ParticleSet2D p = interiorDiskApic(n);
    for (std::size_t k = 0; k < p.size(); ++k) {
        p.v[k] = affineField(p.x[k]);
        p.cx[k] = kCx;
        p.cy[k] = kCy;
    }
    const ParticleSet2D before = p;
    MacMasses m;
    particlesToMac(p, g, m, true);

    // Étape intermédiaire : si ce test échoue, le bug est dans P2G ; sinon, dans G2P.
    real errGrid = 0;
    for (int j = 0; j < g.u.ny(); ++j)
        for (int i = 0; i < g.u.nx(); ++i)
            if (m.u(i, j) > 0) errGrid = std::max(errGrid, std::abs(g.u(i, j) - affineField(g.u.position(i, j)).x));
    for (int j = 0; j < g.v.ny(); ++j)
        for (int i = 0; i < g.v.nx(); ++i)
            if (m.v(i, j) > 0) errGrid = std::max(errGrid, std::abs(g.v(i, j) - affineField(g.v.position(i, j)).y));

    macToParticlesAPIC(g, p);
    const real errV = maxDiff(p.v, before.v);
    const real errC = std::max(maxDiff(p.cx, before.cx), maxDiff(p.cy, before.cy));
    MESSAGE("champ affine APIC : erreur grille = " << errGrid << ", vitesses = " << errV << ", c = " << errC);
    CHECK(errGrid < 1e-12);
    CHECK(errV < 1e-12);
    CHECK(errC < 1e-12);
}

TEST_CASE("APIC : rotation rigide, le moment cinétique est conservé à chaque aller-retour") {
    // Une rotation rigide est un champ affine : APIC la transporte exactement,
    // alors que PIC perd du moment cinétique à chaque aller-retour (jour 7).
    const int n = 32;
    const Vec2 c{0.5, 0.5};
    MacGrid2D g(n, n, 1.0 / n);
    ParticleSet2D p = interiorDiskApic(n);
    for (std::size_t k = 0; k < p.size(); ++k) {
        const Vec2 r = p.x[k] - c;
        p.v[k] = {-r.y, r.x};       // ω = 1
        p.cx[k] = {0, -1};          // (du/dx, du/dy)
        p.cy[k] = {1, 0};           // (dv/dx, dv/dy)
    }
    const real L0 = particleAngularMomentum(p, c), E0 = particleEnergy(p);
    MacMasses m;
    for (int it = 0; it < 10; ++it) {
        particlesToMac(p, g, m, true);
        macToParticlesAPIC(g, p);
    }
    MESSAGE("rotation rigide APIC, 10 allers-retours : L/L0 - 1 = " << particleAngularMomentum(p, c) / L0 - 1);
    CHECK(std::abs(particleAngularMomentum(p, c) / L0 - 1) < 1e-12);
    CHECK(std::abs(particleEnergy(p) / E0 - 1) < 1e-12);
}

TEST_CASE("APIC : P2G conserve la quantité de mouvement (pochoirs complets)") {
    // \sum_f m_f u_f = \sum_p m_p u_p + \sum_p m_p c_p · \sum_f w_fp (x_f - x_p),
    // et le dernier terme est nul : le pochoir reproduit les fonctions affines.
    const int n = 32;
    MacGrid2D g(n, n, 1.0 / n);
    ParticleSet2D p = interiorDiskApic(n);
    std::mt19937 rng(11);
    std::uniform_real_distribution<real> U(-1, 1);
    for (std::size_t k = 0; k < p.size(); ++k) {
        p.v[k] = {U(rng), U(rng)};
        p.cx[k] = {U(rng), U(rng)};   // gradients quelconques : la conservation ne dépend pas de c
        p.cy[k] = {U(rng), U(rng)};
    }
    MacMasses m;
    particlesToMac(p, g, m, true);
    Vec2 Pgrid{0, 0};
    for (std::size_t n2 = 0; n2 < g.u.data().size(); ++n2) Pgrid.x += m.u.data()[n2] * g.u.data()[n2];
    for (std::size_t n2 = 0; n2 < g.v.data().size(); ++n2) Pgrid.y += m.v.data()[n2] * g.v.data()[n2];
    CHECK(norm(Pgrid - particleMomentum(p)) < 1e-12 * particleMass(p));
}

TEST_CASE("APIC : pochoir tronqué au bord, dérivée tangentielle nulle, l'autre exacte") {
    // Pour u, le pochoir sort du champ dans les demi-cellules du bas et du haut (axe y) ;
    // pour v, dans celles de gauche et de droite (axe x). La vitesse y est prolongée
    // par une constante : la dérivée selon l'axe tronqué vaut 0.
    const int n = 16;
    const real h = 1.0 / n;
    MacGrid2D g(n, n, h);
    fillAffine(g);
    ParticleSet2D p;
    p.add({0.37, 0.2 * h}, {0, 0}, 1);       // bas    : u tronqué en y
    p.add({0.37, 1 - 0.2 * h}, {0, 0}, 1);   // haut   : u tronqué en y
    p.add({0.2 * h, 0.61}, {0, 0}, 1);       // gauche : v tronqué en x
    p.add({1 - 0.2 * h, 0.61}, {0, 0}, 1);   // droite : v tronqué en x
    macToParticlesAPIC(g, p);
    for (int k : {0, 1}) {
        CHECK(p.cx[k].y == 0);                                               // du/dy coupé
        CHECK(p.cx[k].x == doctest::Approx(kA[0][0]).epsilon(1e-12));        // du/dx exact
        CHECK(p.cy[k].x == doctest::Approx(kA[1][0]).epsilon(1e-12));        // v complet ici
        CHECK(p.cy[k].y == doctest::Approx(kA[1][1]).epsilon(1e-12));
    }
    for (int k : {2, 3}) {
        CHECK(p.cy[k].x == 0);                                               // dv/dx coupé
        CHECK(p.cy[k].y == doctest::Approx(kA[1][1]).epsilon(1e-12));        // dv/dy exact
        CHECK(p.cx[k].x == doctest::Approx(kA[0][0]).epsilon(1e-12));        // u complet ici
        CHECK(p.cx[k].y == doctest::Approx(kA[0][1]).epsilon(1e-12));
    }
}

TEST_CASE("APIC, régression : avec c = 0, P2G puis G2P donne les mêmes vitesses que PIC") {
    const int n = 32;
    MacGrid2D g1(n, n, 1.0 / n), g2(n, n, 1.0 / n);
    ParticleSet2D p = interiorDiskApic(n);
    std::mt19937 rng(5);
    std::uniform_real_distribution<real> U(-1, 1);
    for (Vec2& v : p.v) v = {U(rng), U(rng)};
    ParticleSet2D pic = p, apic = p;
    MacMasses m1, m2;
    particlesToMac(pic, g1, m1);
    macToParticlesPIC(g1, pic);
    particlesToMac(apic, g2, m2, true);
    macToParticlesAPIC(g2, apic);
    CHECK(maxDiff(pic.v, apic.v) == 0.0);
}

TEST_CASE("APIC, régression : P2G sans terme affine ignore cx et cy") {
    // PIC et FLIP appellent particlesToMac avec affine = false : des c non nuls
    // (restés d'un passage en APIC) ne doivent pas changer la grille.
    const int n = 32;
    MacGrid2D g1(n, n, 1.0 / n), g2(n, n, 1.0 / n);
    ParticleSet2D p = interiorDiskApic(n);
    for (std::size_t k = 0; k < p.size(); ++k) p.v[k] = affineField(p.x[k]);
    ParticleSet2D withC = p;
    for (std::size_t k = 0; k < p.size(); ++k) { withC.cx[k] = kCx; withC.cy[k] = kCy; }
    MacMasses m1, m2;
    particlesToMac(p, g1, m1);
    particlesToMac(withC, g2, m2);
    CHECK(g1.u.data() == g2.u.data());
    CHECK(g1.v.data() == g2.v.data());
}

// ---------------------------------------------------------------------------
//  3. Solveur complet
// ---------------------------------------------------------------------------

TEST_CASE("Solveur APIC, disque en rotation (transferts seuls) : énergie et moment cinétique conservés") {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    prm.scheme = TransferScheme::APIC;
    FlipSolver s(prm, ParticleScene::RotatingDisk);
    // reset() doit avoir donné la partie affine exacte de la rotation.
    const real w = prm.diskAngularVelocity;
    REQUIRE(!s.particles().empty());
    CHECK(s.particles().cx[0].y == doctest::Approx(-w));
    CHECK(s.particles().cy[0].x == doctest::Approx(w));
    const real E0 = particleEnergy(s.particles());
    const real L0 = particleAngularMomentum(s.particles(), s.domainCenter());
    for (int k = 0; k < 50; ++k) s.step(0.01);
    CHECK(particleEnergy(s.particles()) == doctest::Approx(E0).epsilon(1e-12));
    CHECK(particleAngularMomentum(s.particles(), s.domainCenter()) == doctest::Approx(L0).epsilon(1e-12));
}

TEST_CASE("Solveur : c reste à 0 en PIC et en FLIP, y compris pour le disque en rotation") {
    for (TransferScheme sc : {TransferScheme::PIC, TransferScheme::FLIP}) {
        FlipParams prm;
        prm.nx = prm.ny = 24;
        prm.scheme = sc;
        FlipSolver s(prm, ParticleScene::RotatingDisk);
        s.step(0.01);
        real cmax = 0;
        for (std::size_t k = 0; k < s.particles().size(); ++k)
            cmax = std::max({cmax, norm(s.particles().cx[k]), norm(s.particles().cy[k])});
        CHECK(cmax == 0.0);
    }
}

TEST_CASE("Solveur APIC : un fluide au repos reste au repos") {
    FlipParams prm;
    prm.nx = prm.ny = 24;
    prm.scheme = TransferScheme::APIC;
    FlipSolver s(prm, ParticleScene::BoxAtRest);
    for (int k = 0; k < 20; ++k) s.step(0.02);
    CHECK(particleMaxSpeed(s.particles()) == 0.0);
}

TEST_CASE("Solveur APIC, tourbillon en boîte pleine : beaucoup moins dissipatif que PIC") {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    prm.projectionTol = 1e-8;
    auto run = [&](TransferScheme sc, real alpha) {
        FlipParams q = prm;
        q.scheme = sc;
        q.flipRatio = alpha;
        FlipSolver s(q, ParticleScene::BoxVortex);
        const real E0 = particleEnergy(s.particles());
        for (int k = 0; k < 60; ++k) s.step(0.02);
        CHECK(allFiniteParticles(s.particles()));
        return particleEnergy(s.particles()) / E0;
    };
    const real ePic = run(TransferScheme::PIC, 0);
    const real eApic = run(TransferScheme::APIC, 0);
    const real eFlip = run(TransferScheme::FLIP, real(0.95));
    MESSAGE("tourbillon E(1,2 s)/E0 : PIC " << ePic << ", APIC " << eApic << ", FLIP 0.95 " << eFlip);
    CHECK(ePic < eApic);
    CHECK(eApic <= 1 + 1e-6);    // APIC ne crée pas d'énergie dans une boîte fermée
}

TEST_CASE("Solveur APIC, rupture de barrage : l'eau avance, l'énergie ne croît pas, moins de perte que PIC") {
    auto run = [](TransferScheme sc, real& front) {
        FlipParams prm;
        prm.nx = prm.ny = 32;
        prm.scheme = sc;
        prm.dt = real(1) / 120;
        FlipSolver s(prm, ParticleScene::DamBreak);
        auto eTot = [&] { return particleEnergy(s.particles()) + s.potentialEnergy(); };
        const real E0 = eTot();
        for (int k = 0; k < 120; ++k) s.step(prm.dt);   // 1 s
        front = 0;
        for (const Vec2& x : s.particles().x) front = std::max(front, x.x);
        return eTot() / E0;
    };
    real fPic = 0, fApic = 0;
    const real ePic = run(TransferScheme::PIC, fPic);
    const real eApic = run(TransferScheme::APIC, fApic);
    MESSAGE("barrage E(1 s)/E0 : PIC " << ePic << ", APIC " << eApic << " ; front APIC = " << fApic);
    CHECK(fApic > 0.5);
    CHECK(eApic <= 1.01);
    CHECK(ePic < eApic);
}

// ---------------------------------------------------------------------------
//  4. Jour 11 : disque en rotation libre (pas complet, sans gravité)
// ---------------------------------------------------------------------------

TEST_CASE("J11 : disque libre et disque figé partent du même état, seul le comportement change") {
    for (TransferScheme sc : {TransferScheme::PIC, TransferScheme::APIC}) {
        FlipParams prm;
        prm.nx = prm.ny = 32;
        prm.scheme = sc;
        FlipSolver fige(prm, ParticleScene::RotatingDisk);
        FlipSolver libre(prm, ParticleScene::RotatingDiskFree);
        // Le case partagé de reset() : mêmes particules, mêmes vitesses, même partie affine.
        REQUIRE(libre.particles().size() == fige.particles().size());
        CHECK(maxDiff(libre.particles().x, fige.particles().x) == 0.0);
        CHECK(maxDiff(libre.particles().v, fige.particles().v) == 0.0);
        CHECK(maxDiff(libre.particles().cx, fige.particles().cx) == 0.0);
        CHECK(maxDiff(libre.particles().cy, fige.particles().cy) == 0.0);
        // Le comportement : transferts seuls d'un côté, pas complet avec surface libre de l'autre.
        CHECK(fige.transfersOnly());
        CHECK(!fige.hasFreeSurface());
        CHECK(!libre.transfersOnly());
        CHECK(libre.hasFreeSurface());
        CHECK(!libre.usesGravity());
    }
}

TEST_CASE("J11 : après un pas, le disque libre a bougé, le disque figé non") {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    FlipSolver fige(prm, ParticleScene::RotatingDisk);
    FlipSolver libre(prm, ParticleScene::RotatingDiskFree);
    const std::vector<Vec2> x0 = libre.particles().x;
    fige.step(prm.dt);
    libre.step(prm.dt);
    CHECK(maxDiff(fige.particles().x, x0) == 0.0);
    // Le bord tourne à omega R = 0,6 m/s : en dt = 1/60 s, il parcourt 0,01, soit un tiers de cellule.
    CHECK(maxDiff(libre.particles().x, x0) > 0.005);
    CHECK(libre.lastStats().pressureIterations > 0);   // la projection a bien tourné
}

TEST_CASE("J11, diagnostics : centre de masse, R_g, L0 et E0 d'un disque uniforme proches des valeurs exactes") {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    FlipSolver s(prm, ParticleScene::RotatingDiskFree);
    const real R = prm.diskRadius, w = prm.diskAngularVelocity, rho = prm.density;
    const Vec2 G = s.centerOfMass();
    // Semis stratifié aléatoire : le centre bouge un peu, bien moins qu'une cellule.
    CHECK(norm(G - s.domainCenter()) < 0.05 * s.grid().h());
    // Disque plein : R_g = R / sqrt(2) ; L = rho omega pi R^4 / 2 ; E = rho omega^2 pi R^4 / 4.
    // L'écart (quelques 1e-3) vient du bord en escalier du semis.
    CHECK(s.initialRadiusOfGyration() == doctest::Approx(R / std::sqrt(real(2))).epsilon(0.01));
    CHECK(particleAngularMomentum(s.particles(), G) == doctest::Approx(rho * w * kPi * std::pow(R, 4) / 2).epsilon(0.01));
    CHECK(particleEnergy(s.particles()) == doctest::Approx(rho * w * w * kPi * std::pow(R, 4) / 4).epsilon(0.01));
    CHECK(norm(particleMomentum(s.particles())) < 1e-3 * particleMass(s.particles()) * w * R);
}

TEST_CASE("J11, L de la grille après P2G : égal à celui des particules en PIC, plus grand en APIC (part affine)") {
    FlipParams prm;
    prm.nx = prm.ny = 32;
    // PIC : P2G conserve exactement L (pochoirs complets, disque loin des parois).
    prm.scheme = TransferScheme::PIC;
    FlipSolver pic(prm, ParticleScene::RotatingDisk);
    const real Lpic = particleAngularMomentum(pic.particles(), pic.centerOfMass());
    pic.step(0.01);
    CHECK(pic.lastStats().angularMomentumGrid == doctest::Approx(Lpic).epsilon(1e-10));
    // APIC : la grille reçoit aussi la rotation portée par cx, cy, que la somme
    // des x_p × m_p v_p ne voit pas. Elle est donc strictement plus grande ici.
    prm.scheme = TransferScheme::APIC;
    FlipSolver apic(prm, ParticleScene::RotatingDisk);
    const real Lapic = particleAngularMomentum(apic.particles(), apic.centerOfMass());
    apic.step(0.01);
    const real ratio = apic.lastStats().angularMomentumGrid / Lapic;
    MESSAGE("APIC, disque fige : L grille / L particules = " << ratio);
    CHECK(ratio > 1 + 1e-3);
    CHECK(ratio < 1.05);
}

TEST_CASE("J11, disque libre sur un quart de tour : la forme tient, APIC garde L bien mieux que PIC") {
    // Un disque de fluide libre en rotation est instable (Rayleigh-Taylor : la force
    // centrifuge pousse le fluide vers l'air, sans tension de surface pour le retenir).
    // On ne le teste donc que sur un quart de tour, avant que les doigts ne se forment.
    struct Result { real L, E, Emax, Rg, drift; };
    auto run = [](TransferScheme sc, real alpha) {
        FlipParams prm;
        prm.nx = prm.ny = 32;
        prm.scheme = sc;
        prm.flipRatio = alpha;
        prm.dt = real(1) / 60;
        FlipSolver s(prm, ParticleScene::RotatingDiskFree);
        const Vec2 G0 = s.centerOfMass();
        const real E0 = particleEnergy(s.particles());
        s.step(prm.dt);
        const real Lg0 = s.lastStats().angularMomentumGrid;   // L de la grille au premier P2G
        const int steps = int(std::lround(kPi / (2 * prm.diskAngularVelocity) / prm.dt));   // quart de tour
        Result r{0, 0, 0, 0, 0};
        for (int k = 1; k < steps; ++k) {
            s.step(prm.dt);
            r.Emax = std::max(r.Emax, particleEnergy(s.particles()) / E0);
        }
        CHECK(allFiniteParticles(s.particles()));
        r.L = s.lastStats().angularMomentumGrid / Lg0;
        r.E = particleEnergy(s.particles()) / E0;
        r.Rg = particleRadiusOfGyration(s.particles(), s.centerOfMass()) / s.initialRadiusOfGyration();
        r.drift = norm(s.centerOfMass() - G0) / s.grid().h();
        return r;
    };
    const Result pic = run(TransferScheme::PIC, 0);
    const Result flip = run(TransferScheme::FLIP, 1);
    const Result apic = run(TransferScheme::APIC, 0);
    MESSAGE("quart de tour, L/L0 : PIC " << pic.L << ", FLIP " << flip.L << ", APIC " << apic.L
            << " ; E/E0 : PIC " << pic.E << ", FLIP " << flip.E << ", APIC " << apic.E
            << " ; R_g/R_g0 APIC " << apic.Rg);
    for (const Result& r : {pic, flip, apic}) {
        CHECK(r.Emax <= 1 + 1e-3);              // pas de création d'énergie
        CHECK(std::abs(r.Rg - 1) < 0.03);       // encore un disque
        CHECK(r.drift < 0.5);                   // pas de dérive (P = 0)
    }
    CHECK(apic.L > 0.95);                       // mesuré : environ 0,97
    CHECK(apic.L > pic.L + 0.1);                // mesuré : PIC environ 0,73 en 32 x 32
    CHECK(flip.L > pic.L + 0.1);
}
