#include "numerics/MacOperators.h"
#include <algorithm>
#include <cmath>

void enforceWallBC(MacGrid2D& g) {
    //Condition de Neumann sur le bord de la grille mac
    const int nx = g.nx(), ny = g.ny();
    for (int j = 0; j < ny; ++j) { g.u(0, j) = 0; g.u(nx, j) = 0; }
    for (int i = 0; i < nx; ++i) { g.v(i, 0) = 0; g.v(i, ny) = 0; }
}

void computeDivergence(const MacGrid2D& g, Field2D& div) {
    //Calcule la divergence sur la macgrid g : 
    // div(i, j) = (u(i + 1, j) - u(i, j) + v(i, j + 1) - v(i, j)) / h(div : champ de cellules)
    assert(div.nx() == g.nx() && div.ny() == g.ny());
    for (int j = 0; j < g.ny(); ++j)
        for (int i = 0; i < g.nx(); ++i)
            div(i, j) = g.divergence(i, j);
}

real maxAbsDivergence(const MacGrid2D& g) {
    //max de |div| sur chaque cellule
    real m = 0;
    for (int j = 0; j < g.ny(); ++j)
        for (int i = 0; i < g.nx(); ++i)
            m = std::max(m, std::abs(g.divergence(i, j)));
    return m;
}

void computeGradient(const Field2D& q, Field2D& gx, Field2D& gy) {
    // Gradient d'un champ de cellules q, rangé sur les faces :
    //   gx(i,j) = (q(i,j) - q(i-1,j)) / h   pour 1 <= i <= nx-1, et 0 sur les parois
    //   gy(i,j) = (q(i,j) - q(i,j-1)) / h   pour 1 <= j <= ny-1, et 0 sur les parois
    // gx a le layout de u, gy celui de v.
    const int nx = q.nx(), ny = q.ny();
    const real invh = 1 / q.h();
    assert(gx.nx() == nx + 1 && gx.ny() == ny);
    assert(gy.nx() == nx && gy.ny() == ny + 1);
    for (int j = 0; j < ny; ++j) {
        gx(0, j) = 0; gx(nx, j) = 0;
        for (int i = 1; i < nx; ++i) gx(i, j) = (q(i, j) - q(i - 1, j)) * invh;
    }
    for (int i = 0; i < nx; ++i) { gy(i, 0) = 0; gy(i, ny) = 0; }
    for (int j = 1; j < ny; ++j)
        for (int i = 0; i < nx; ++i) gy(i, j) = (q(i, j) - q(i, j - 1)) * invh;
}

// Code qui correspond à faire : u <- u - scale * G q  (faces intérieures seulement ; les parois gardent u·n = 0).
void subtractGradient(MacGrid2D& g, const Field2D& q, real scale) {
    const int nx = g.nx(), ny = g.ny();
    const real s = scale / g.h();
    for (int j = 0; j < ny; ++j)
        for (int i = 1; i < nx; ++i) g.u(i, j) -= s * (q(i, j) - q(i - 1, j));
    for (int j = 1; j < ny; ++j)
        for (int i = 0; i < nx; ++i) g.v(i, j) -= s * (q(i, j) - q(i, j - 1));
}

//Laplacien négatif 
void NegLaplacianMAC::operator()(const std::vector<real>& x, std::vector<real>& y) const {
    const real invh2 = 1 / (h * h);
    y.resize(x.size());
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const size_t k = size_t(i) + size_t(nx) * size_t(j);
            int  n   = 0;      // nombre de voisins fluides (les parois ne comptent pas)
            real sum = 0;
            if (i > 0)      { ++n; sum += x[k - 1]; }
            if (i < nx - 1) { ++n; sum += x[k + 1]; }
            if (j > 0)      { ++n; sum += x[k - nx]; }
            if (j < ny - 1) { ++n; sum += x[k + nx]; }
            y[k] = (n * x[k] - sum) * invh2;
        }
}
