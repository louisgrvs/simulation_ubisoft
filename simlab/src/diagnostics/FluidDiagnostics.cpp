#include "diagnostics/FluidDiagnostics.h"
#include <cmath>
#include <vector>

real kineticEnergy(const MacGrid2D& g) {
    real s = 0;
    for (real a : g.u.data()) s += a * a;
    for (real a : g.v.data()) s += a * a;
    return real(0.5) * g.h() * g.h() * s;
}

real totalAmount(const Field2D& q) { return q.h() * q.h() * q.sum(); }

void computeSpeed(const MacGrid2D& g, Field2D& out) {
    for (int j = 0; j < g.ny(); ++j)
        for (int i = 0; i < g.nx(); ++i) {
            const real uc = real(0.5) * (g.u(i, j) + g.u(i + 1, j));
            const real vc = real(0.5) * (g.v(i, j) + g.v(i, j + 1));
            out(i, j) = std::sqrt(uc * uc + vc * vc);
        }
}

void computeVorticity(const MacGrid2D& g, Field2D& out) {
    const int nx = g.nx(), ny = g.ny();
    const real invh = 1 / g.h();
    // Coins (i, j), 0 <= i <= nx, 0 <= j <= ny ; nuls sur le bord.
    std::vector<real> corner(size_t(nx + 1) * (ny + 1), 0);
    for (int j = 1; j < ny; ++j)
        for (int i = 1; i < nx; ++i)
            corner[size_t(i) + size_t(nx + 1) * j] =
                (g.v(i, j) - g.v(i - 1, j)) * invh - (g.u(i, j) - g.u(i, j - 1)) * invh;
    auto c = [&](int i, int j) { return corner[size_t(i) + size_t(nx + 1) * j]; };
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            out(i, j) = real(0.25) * (c(i, j) + c(i + 1, j) + c(i, j + 1) + c(i + 1, j + 1));
}
