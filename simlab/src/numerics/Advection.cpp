#include "numerics/Advection.h"
#include <algorithm>
#include <limits>

Vec2 traceBack(const MacGrid2D& grid, Vec2 x, real dt, Integrator scheme) {
    Vec2 xBack;
    if (scheme == Integrator::Euler) {
        xBack = x - dt * grid.velocity(x);
    } else {
        Vec2 xMid = x - (real(0.5) * dt) * grid.velocity(x);
        xBack     = x - dt * grid.velocity(xMid);
    }
    // Une trajectoire qui sort du domaine lit la valeur au bord.
    // C'est ICI qu'il faudra raffiner quand viendront murs et obstacles.
    Vec2 L = grid.size();
    xBack.x = std::clamp(xBack.x, real(0), L.x);
    xBack.y = std::clamp(xBack.y, real(0), L.y);
    return xBack;
}

void advectScalar(const MacGrid2D& grid, const Field2D& q, Field2D& qOut, real dt,
                  Integrator scheme) {
    assert(q.sameLayout(qOut) && &q != &qOut);
    for (int j = 0; j < q.ny(); ++j)
        for (int i = 0; i < q.nx(); ++i)
            qOut(i, j) = q.sample(traceBack(grid, q.position(i, j), dt, scheme));
}

real cflTimeStep(const MacGrid2D& grid, real cfl) {
    real umax = grid.maxSpeed();
    if (umax <= 0) return std::numeric_limits<real>::infinity();
    return cfl * grid.h() / umax;
}
