#include "numerics/ConjugateGradient.h"
#include <algorithm>
#include <cassert>
#include <cmath>

//produit scalaire entre deux vecteurs
real dotProduct(const Vector& a, const Vector& b) {
    assert(a.size() == b.size());
    real s = 0;
    for (size_t k = 0; k < a.size(); ++k) s += a[k] * b[k];
    return s;
}

//norme l infty de a
real normInf(const Vector& a) {
    real m = 0;
    for (real v : a) m = std::max(m, std::abs(v));
    return m;
}

// fait l'opération a <- a - mean(a)
void subtractMean(Vector& a) {
    if (a.empty()) return;
    real mean = 0;
    for (real v : a) mean += v;
    mean /= real(a.size());
    for (real& v : a) v -= mean;
}

//Un algo pour résoudre Ax = y  (faut trouver x). Pour ce faire, on minimise la formulation variationnelle associée. 
//C'est utile pour obtenir la pression et pour la diffusion implicite (cf le papier stable fluid).
CGResult conjugateGradient(const LinearOperator& A, const Vector& bIn, Vector& x,
                           const CGOptions& opt) {
    const size_t n = bIn.size();
    x.resize(n, 0);

    Vector b = bIn;
    if (opt.removeMean) subtractMean(b);        // b ∈ Im(A) = (constantes)^\perp

    // r0 = b - A x0
    Vector r(n), Ap(n);
    A(x, Ap);
    for (size_t k = 0; k < n; ++k) r[k] = b[k] - Ap[k];
    if (opt.removeMean) subtractMean(r);

    CGResult res;
    res.residual = normInf(r);
    if (opt.history) { opt.history->clear(); opt.history->push_back(res.residual); }
    if (res.residual <= opt.tol) { res.converged = true; return res; }

    Vector p = r;
    real rr = dotProduct(r, r);

    for (int it = 1; it <= opt.maxIter; ++it) {
        A(p, Ap);
        const real pAp = dotProduct(p, Ap);
        if (pAp <= 0) break;                    // direction dans le noyau : rien à gagner
        const real alpha = rr / pAp;

        for (size_t k = 0; k < n; ++k) { x[k] += alpha * p[k]; r[k] -= alpha * Ap[k]; }

        res.iterations = it;
        res.residual   = normInf(r);
        if (opt.history) opt.history->push_back(res.residual);
        if (res.residual <= opt.tol) { res.converged = true; break; }

        const real rrNew = dotProduct(r, r);
        const real beta  = rrNew / rr;
        rr = rrNew;
        for (size_t k = 0; k < n; ++k) p[k] = r[k] + beta * p[k];
    }

    if (opt.removeMean) subtractMean(x);        // on fixe la constante : moyenne nulle
    return res;
}
