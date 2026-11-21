#include "diagnostics/ParticleDiagnostics.h"
#include <algorithm>
#include <cmath>
#include <vector>

real particleMass(const ParticleSet2D& p) {
    real s = 0;
    for (real m : p.m) s += m;
    return s;
}

Vec2 particleMomentum(const ParticleSet2D& p) {
    Vec2 s;
    for (std::size_t k = 0; k < p.size(); ++k) s += p.m[k] * p.v[k];
    return s;
}

real particleEnergy(const ParticleSet2D& p) {
    real s = 0;
    for (std::size_t k = 0; k < p.size(); ++k) s += p.m[k] * dot(p.v[k], p.v[k]);
    return real(0.5) * s;
}

real particleAngularMomentum(const ParticleSet2D& p, Vec2 c) {
    real s = 0;
    for (std::size_t k = 0; k < p.size(); ++k) s += p.m[k] * cross(p.x[k] - c, p.v[k]);
    return s;
}

real particleMaxSpeed(const ParticleSet2D& p) {
    real s = 0;
    for (const Vec2& v : p.v) s = std::max(s, norm(v));
    return s;
}

Vec2 gridMomentum(const MacGrid2D& g, const MacMasses& m) {
    Vec2 s;
    for (std::size_t n = 0; n < g.u.data().size(); ++n) s.x += m.u.data()[n] * g.u.data()[n];
    for (std::size_t n = 0; n < g.v.data().size(); ++n) s.y += m.v.data()[n] * g.v.data()[n];
    return s;
}

real gridEnergy(const MacGrid2D& g, const MacMasses& m) {
    real s = 0;
    for (std::size_t n = 0; n < g.u.data().size(); ++n) s += m.u.data()[n] * g.u.data()[n] * g.u.data()[n];
    for (std::size_t n = 0; n < g.v.data().size(); ++n) s += m.v.data()[n] * g.v.data()[n] * g.v.data()[n];
    return real(0.5) * s;
}

real gridAngularMomentum(const MacGrid2D& g, const MacMasses& m, Vec2 c) {
    real s = 0;
    for (int j = 0; j < g.u.ny(); ++j)
        for (int i = 0; i < g.u.nx(); ++i)
            s -= m.u(i, j) * (g.u.position(i, j).y - c.y) * g.u(i, j);
    for (int j = 0; j < g.v.ny(); ++j)
        for (int i = 0; i < g.v.nx(); ++i)
            s += m.v(i, j) * (g.v.position(i, j).x - c.x) * g.v(i, j);
    return s;
}

void countParticlesPerCell(const ParticleSet2D& p, Field2D& count) {
    count.fill(0);
    const real h = count.h();
    for (const Vec2& x : p.x) {
        const int i = std::clamp(int(std::floor(x.x / h)), 0, count.nx() - 1);
        const int j = std::clamp(int(std::floor(x.y / h)), 0, count.ny() - 1);
        count(i, j) += 1;
    }
}

// J9 -----------------------------------------------------------------------------

real particlePotentialEnergy(const ParticleSet2D& p, Vec2 gravity) {
    real s = 0;
    for (std::size_t k = 0; k < p.size(); ++k) s -= p.m[k] * dot(gravity, p.x[k]);
    return s;
}

NoiseStats cellVelocityVariance(const ParticleSet2D& p, Field2D& variance, int minCount) {
    const int nx = variance.nx(), ny = variance.ny();
    const real h = variance.h();
    const std::size_t nc = std::size_t(nx) * ny;
    std::vector<int>  cell(p.size());
    std::vector<int>  count(nc, 0);
    std::vector<real> mass(nc, 0);
    std::vector<Vec2> mean(nc);

    // 1er passage : masse et vitesse moyenne de chaque cellule.
    for (std::size_t k = 0; k < p.size(); ++k) {
        const int i = std::clamp(int(std::floor(p.x[k].x / h)), 0, nx - 1);
        const int j = std::clamp(int(std::floor(p.x[k].y / h)), 0, ny - 1);
        const int c = i + nx * j;
        cell[k] = c;
        ++count[c];
        mass[c] += p.m[k];
        mean[c] += p.m[k] * p.v[k];
    }
    for (std::size_t c = 0; c < nc; ++c)
        if (mass[c] > 0) mean[c] = mean[c] / mass[c];

    // 2e passage : écarts à la moyenne (deux passages : pas de soustraction de
    // grands nombres voisins, contrairement à <|u|²> - |<u>|²).
    std::vector<real> sq(nc, 0);
    real m2 = 0, mtot = 0;
    for (std::size_t k = 0; k < p.size(); ++k) {
        const int c = cell[k];
        const Vec2 d = p.v[k] - mean[c];
        sq[c] += p.m[k] * dot(d, d);
        m2 += p.m[k] * dot(p.v[k], p.v[k]);
        mtot += p.m[k];
    }

    NoiseStats st;
    real num = 0, den = 0;
    variance.fill(0);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int c = i + nx * j;
            if (count[c] < minCount) continue;
            variance(i, j) = sq[c] / mass[c];
            num += sq[c];
            den += mass[c];
            ++st.cells;
        }
    st.meanVariance = den > 0 ? num / den : 0;
    const real meanSq = mtot > 0 ? m2 / mtot : 0;
    st.relative = meanSq > 0 ? st.meanVariance / meanSq : 0;
    return st;
}

// J11 ----------------------------------------------------------------------------

Vec2 particleCenterOfMass(const ParticleSet2D& p) {
    Vec2 s;
    real m = 0;
    for (std::size_t k = 0; k < p.size(); ++k) {
        s += p.m[k] * p.x[k];
        m += p.m[k];
    }
    return m > 0 ? s / m : Vec2{};
}

real particleRadiusOfGyration(const ParticleSet2D& p, Vec2 c) {
    real s = 0, m = 0;
    for (std::size_t k = 0; k < p.size(); ++k) {
        const Vec2 r = p.x[k] - c;
        s += p.m[k] * dot(r, r);
        m += p.m[k];
    }
    return m > 0 ? std::sqrt(s / m) : real(0);
}