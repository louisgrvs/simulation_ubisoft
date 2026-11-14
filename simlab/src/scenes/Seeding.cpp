#include "scenes/Seeding.h"
#include <cassert>
#include <random>

namespace region {

Region box(Vec2 lo, Vec2 hi) {
    return [=](Vec2 p) { return p.x >= lo.x && p.x < hi.x && p.y >= lo.y && p.y < hi.y; };
}

Region disk(Vec2 center, real radius) {
    return [=](Vec2 p) { const Vec2 d = p - center; return dot(d, d) < radius * radius; };
}

Region everywhere() {
    return [](Vec2) { return true; };
}

} // namespace region

namespace {

// SImulation d'un réel uniforme dans [0, 1[, à partir des 32 bits bruts de mt19937 (portable).
real uniform01(std::mt19937& rng) {
    return real(rng()) * (real(1) / real(4294967296.0));
}

} // namespace

int seedParticles(ParticleSet2D& particles, int nx, int ny, real h, const Region& inside,
                  const SeedOptions& opt, const std::function<Vec2(Vec2)>& velocity) {
    assert(opt.perAxis >= 1 && opt.jitter >= 0 && opt.jitter <= 1 && h > 0);
    const int k = opt.perAxis;
    const real sub = h / k;                       // côté d'une sous-cellule
    const real mass = opt.density * sub * sub;    // \rho h^2 / k^2
    std::mt19937 rng(opt.seed);

    int added = 0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            for (int b = 0; b < k; ++b)
                for (int a = 0; a < k; ++a) {
                    // On tire tjr les deux nombres, même si la sous-cellule est
                    // rejetée : ainsi le tirage d'une sous-cellule ne dépend pas de
                    // la forme de la région (changer de scène ne décale pas la suite).
                    const real ux = uniform01(rng), uy = uniform01(rng); // tirage d'un nombre uniforme dans [0;1[
                    const Vec2 center{i * h + (a + real(0.5)) * sub, j * h + (b + real(0.5)) * sub}; //calcule le centre de la sous-cellule
                    if (!inside(center)) continue; // teste si ce centre est dans la région
                    const Vec2 pos{center.x + opt.jitter * (ux - real(0.5)) * sub, // décale le point à partir du centre
                                   center.y + opt.jitter * (uy - real(0.5)) * sub};
                    particles.add(pos, velocity ? velocity(pos) : Vec2{}, mass);
                    ++added; //ajoute la particule et calcul sa masse(cf particleset2d.h)
                }
    return added;
}
