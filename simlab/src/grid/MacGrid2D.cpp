#include "grid/MacGrid2D.h"
#include <algorithm>
#include <cmath>

MacGrid2D::MacGrid2D(int nx, int ny, real h)
    : u(nx + 1, ny, h, 0.0, 0.5), //Field2D verticale
      v(nx, ny + 1, h, 0.5, 0.0), //Field2D horizontale
      p(nx, ny, h, 0.5, 0.5),     //Field2D de pression (au centre)
      nx_(nx), ny_(ny), h_(h) {}

void MacGrid2D::setVelocity(const std::function<Vec2(Vec2)>& f) {
    for (int j = 0; j < u.ny(); ++j)
        for (int i = 0; i < u.nx(); ++i) u(i, j) = f(u.position(i, j)).x;
    for (int j = 0; j < v.ny(); ++j)
        for (int i = 0; i < v.nx(); ++i) v(i, j) = f(v.position(i, j)).y;
}

real MacGrid2D::maxSpeed() const {
    real m = 0;
    for (real a : u.data()) m = std::max(m, std::abs(a));
    for (real a : v.data()) m = std::max(m, std::abs(a));
    return m;
}

real MacGrid2D::divergence(int i, int j) const {
    return (u(i + 1, j) - u(i, j) + v(i, j + 1) - v(i, j)) / h_;
}

real MacGrid2D::maxAbsDivergence() const {
    real m = 0;
    for (int j = 0; j < ny_; ++j)
        for (int i = 0; i < nx_; ++i) m = std::max(m, std::abs(divergence(i, j)));
    return m;
}
