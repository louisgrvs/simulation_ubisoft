#include "numerics/Diffusion.h"
#include <algorithm>
#include <cassert>

void HelmholtzOperator::operator()(const Vector& x, Vector& y) const {
    y.resize(x.size());
    const bool dirX = xEdges == EdgeCondition::Dirichlet0;
    const bool dirY = yEdges == EdgeCondition::Dirichlet0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const size_t k = size_t(i) + size_t(nx) * size_t(j);
            int  n   = 0;
            real sum = 0;
            // Voisin présent : il compte et on l'ajoute. Voisin absent : il compte
            // seulement en Dirichlet (valeur 0 à distance h), et n'ajoute rien.
            if (i > 0)      { ++n; sum += x[k - 1]; }  else if (dirX) ++n;
            if (i < nx - 1) { ++n; sum += x[k + 1]; }  else if (dirX) ++n;
            if (j > 0)      { ++n; sum += x[k - nx]; } else if (dirY) ++n;
            if (j < ny - 1) { ++n; sum += x[k + nx]; } else if (dirY) ++n;
            y[k] = x[k] + c * (n * x[k] - sum);
        }
}

namespace {

// Résout (I - c Δ_h) x = b, en partant de x = b (bonne estimation quand c est petit).
DiffusionStats solveHelmholtz(const HelmholtzOperator& A, Vector& x, real tolRel) {
    DiffusionStats st;
    const Vector b = x;
    const real scale = normInf(b);
    if (A.c <= 0 || scale == 0 || b.empty()) return st;   // rien à faire

    CGOptions opt;
    opt.tol     = tolRel * scale;
    opt.maxIter = int(b.size());
    CGResult r  = conjugateGradient(A, b, x, opt);
    st.iterations = r.iterations;
    st.residual   = r.residual / scale;
    st.converged  = r.converged;
    return st;
}

} // namespace

DiffusionStats diffuseCells(Field2D& q, real kappa, real dt, real tol) {
    if (kappa <= 0 || dt <= 0) return {};
    HelmholtzOperator A;
    A.nx = q.nx();
    A.ny = q.ny();
    A.c  = kappa * dt / (q.h() * q.h());
    A.xEdges = A.yEdges = EdgeCondition::Neumann;
    return solveHelmholtz(A, q.data(), tol);
}

DiffusionStats diffuseVelocity(MacGrid2D& g, real nu, real dt, real tol) {
    if (nu <= 0 || dt <= 0) return {};
    const int nx = g.nx(), ny = g.ny();
    const real c = nu * dt / (g.h() * g.h());
    DiffusionStats total;

    // --- u : inconnues sur les faces intérieures i = 1..nx-1, j = 0..ny-1 ---------
    if (nx > 1) {
        HelmholtzOperator A{nx - 1, ny, c, EdgeCondition::Dirichlet0, EdgeCondition::Neumann};
        Vector x(size_t(nx - 1) * ny);
        for (int j = 0; j < ny; ++j)
            for (int i = 1; i < nx; ++i) x[size_t(i - 1) + size_t(nx - 1) * j] = g.u(i, j);
        DiffusionStats s = solveHelmholtz(A, x, tol);
        for (int j = 0; j < ny; ++j)
            for (int i = 1; i < nx; ++i) g.u(i, j) = x[size_t(i - 1) + size_t(nx - 1) * j];
        total.iterations += s.iterations;
        total.residual = std::max(total.residual, s.residual);
        total.converged = total.converged && s.converged;
    }
    // --- v : inconnues sur les faces intérieures i = 0..nx-1, j = 1..ny-1 ---------
    if (ny > 1) {
        HelmholtzOperator A{nx, ny - 1, c, EdgeCondition::Neumann, EdgeCondition::Dirichlet0};
        Vector x(size_t(nx) * (ny - 1));
        for (int j = 1; j < ny; ++j)
            for (int i = 0; i < nx; ++i) x[size_t(i) + size_t(nx) * (j - 1)] = g.v(i, j);
        DiffusionStats s = solveHelmholtz(A, x, tol);
        for (int j = 1; j < ny; ++j)
            for (int i = 0; i < nx; ++i) g.v(i, j) = x[size_t(i) + size_t(nx) * (j - 1)];
        total.iterations += s.iterations;
        total.residual = std::max(total.residual, s.residual);
        total.converged = total.converged && s.converged;
    }
    return total;
}

void dissipate(Field2D& q, real alpha, real dt) {
    if (alpha <= 0 || dt <= 0) return;
    const real factor = 1 / (1 + alpha * dt);
    for (real& a : q.data()) a *= factor;
}
