#include "numerics/Forces.h"
#include <cmath>

real brushWeight(Vec2 x, Vec2 center, real radius) {
    const Vec2 d = x - center;
    return std::exp(-dot(d, d) / (radius * radius));
}

void addBrushForce(MacGrid2D& g, Vec2 center, real radius, Vec2 accel, real dt) {
    const int nx = g.nx(), ny = g.ny();
    for (int j = 0; j < ny; ++j)
        for (int i = 1; i < nx; ++i)                 // faces u intérieures
            g.u(i, j) += dt * accel.x * brushWeight(g.u.position(i, j), center, radius);
    for (int j = 1; j < ny; ++j)
        for (int i = 0; i < nx; ++i)                 // faces v intérieures
            g.v(i, j) += dt * accel.y * brushWeight(g.v.position(i, j), center, radius);
}

void addBrushSource(Field2D& q, Vec2 center, real radius, real rate, real dt) {
    for (int j = 0; j < q.ny(); ++j)
        for (int i = 0; i < q.nx(); ++i)
            q(i, j) += dt * rate * brushWeight(q.position(i, j), center, radius);
}

void addBuoyancy(MacGrid2D& g, const Field2D& smoke, real beta, real dt) {
    assert(smoke.nx() == g.nx() && smoke.ny() == g.ny());
    for (int j = 1; j < g.ny(); ++j)
        for (int i = 0; i < g.nx(); ++i) {
            const real sFace = real(0.5) * (smoke(i, j - 1) + smoke(i, j));
            g.v(i, j) += dt * beta * sFace;
        }
}

void addGravity(MacGrid2D& g, Vec2 gravity, real dt) {
    if (gravity.x != 0) for (real& u : g.u.data()) u += dt * gravity.x;
    if (gravity.y != 0) for (real& v : g.v.data()) v += dt * gravity.y;
}