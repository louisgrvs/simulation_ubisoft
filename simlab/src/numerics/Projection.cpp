#include "numerics/Projection.h"
#include "numerics/MacOperators.h"
#include <cassert>

//Algorithme de porjection comme dans le STAM, on résout \Delta q = div(w_3) et w_4 = w_3 - \nabla q
ProjectionStats project(MacGrid2D& g, real dt, const ProjectionOptions& opt,
                        Vector* residualHistory) {
    assert(dt > 0 && opt.density > 0);
    const int nx = g.nx(), ny = g.ny();
    const real scale = dt / opt.density;        // q = (dt/\rho) p

    // 1) u.n = 0 AVANT de calculer la divergence : c'est ce qui rend sum(D w) = 0.
    enforceWallBC(g);

    ProjectionStats st;
    Field2D div = g.makeCellField();
    computeDivergence(g, div);

    // 2) second membre b = -D w
    Vector b(div.data().size());
    for (size_t k = 0; k < b.size(); ++k) b[k] = -div.data()[k];
    st.divBefore = normInf(b);

    // 3) estimation initiale : pression du pas précédent, ou zéro
    Field2D q = g.makeCellField();
    if (opt.warmStart)
        for (size_t k = 0; k < q.data().size(); ++k) q.data()[k] = scale * g.p.data()[k];

    // 4) A q = b par gradient conjugué sans matrice
    CGOptions cg;
    cg.tol        = opt.tol;
    cg.maxIter    = opt.maxIter > 0 ? opt.maxIter : nx * ny;
    cg.removeMean = true;                       // Neumann pur : noyau = constantes
    cg.history    = residualHistory;
    NegLaplacianMAC A{nx, ny, g.h()};
    CGResult r = conjugateGradient(A, b, q.data(), cg);
    st.iterations = r.iterations;
    st.converged  = r.converged;

    // 5) u = w - G q  (faces intérieures seulement)
    subtractGradient(g, q, 1);

    // 6) pression physique rangée dans la grille
    for (size_t k = 0; k < q.data().size(); ++k) g.p.data()[k] = q.data()[k] / scale;

    st.divAfter = maxAbsDivergence(g);
    return st;
}
