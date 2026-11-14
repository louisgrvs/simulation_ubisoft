#include "grid/Field2D.h"
#include <algorithm>
#include <cmath>

Field2D::Field2D(int nx, int ny, real h, real ox, real oy, real init)
    : nx_(nx), ny_(ny), h_(h), ox_(ox), oy_(oy),
      data_(std::size_t(nx) * std::size_t(ny), init) {
    assert(nx >= 1 && ny >= 1 && h > 0);
}

real Field2D::sample(Vec2 x) const {
    // 1) Coordonnées physiques -> coordonnées « indice » (réelles) de CE champ.
    //    C'est ici, et seulement ici, que le décalage (ox, oy) intervient.
    real gx = x.x / h_ - ox_;
    real gy = x.y / h_ - oy_;

    // 2) On borne au rectangle des échantillons [0, nx-1] x [0, ny-1].
    gx = std::clamp(gx, real(0), real(nx_ - 1));
    gy = std::clamp(gy, real(0), real(ny_ - 1));

    // 3) Cellule d'interpolation (i0, j0) et coordonnées locales (fx, fy) dans [0, 1].
    //    Le min(…, n-2) gère le bord droit/haut (gx = nx-1 tombe dans la dernière cellule).
    int i0 = std::clamp(int(std::floor(gx)), 0, std::max(nx_ - 2, 0));
    int j0 = std::clamp(int(std::floor(gy)), 0, std::max(ny_ - 2, 0));
    int i1 = std::min(i0 + 1, nx_ - 1);
    int j1 = std::min(j0 + 1, ny_ - 1);
    real fx = gx - i0;
    real fy = gy - j0;

    // 4) Deux interpolations linéaires en x, puis une en y.
    const Field2D& f = *this;
    real bottom = (1 - fx) * f(i0, j0) + fx * f(i1, j0);
    real top    = (1 - fx) * f(i0, j1) + fx * f(i1, j1);
    return (1 - fy) * bottom + fy * top;
}

void Field2D::fill(real value) { std::fill(data_.begin(), data_.end(), value); }

bool Field2D::sameLayout(const Field2D& o) const {
    return nx_ == o.nx_ && ny_ == o.ny_ && h_ == o.h_ && ox_ == o.ox_ && oy_ == o.oy_;
}

real Field2D::minValue() const { return *std::min_element(data_.begin(), data_.end()); }
real Field2D::maxValue() const { return *std::max_element(data_.begin(), data_.end()); }

real Field2D::sum() const {
    real s = 0;
    for (real a : data_) s += a;
    return s;
}
