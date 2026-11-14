#include "diagnostics/ParticleDiagnostics.h"
#include <algorithm>
#include <cmath>

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
