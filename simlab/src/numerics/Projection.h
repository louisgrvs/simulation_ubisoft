#pragma once
#include "grid/MacGrid2D.h"
#include "numerics/ConjugateGradient.h"
#include "grid/CellLabels.h"

// ============================================================================
//  Projection de Helmholtz-Hodge discrète (Stam, éq. 3-5 ; Zhu-Bridson #4.2.3)
// ============================================================================
//  Entrée : w (vitesse sur les faces, éventuellement divergente).
//  On cherche p tel que  u = w - (dt/\rho) G p  vérifie  D u = 0  et  u.n = 0.
//
//      D w - (dt/\rho) D G p = 0   <=>   A q = -D w,   avec q = (dt/ρ) p,  A = -D G
//
//  Puis u = w - G q, et on range p = (\rho/dt) q dans grid.p (pression physique,
//  qui servira de démarrage à chaud au pas suivant et au sable).
//
//  Écart assumé avec Stam : sur sa grille colocalisée, D\cdotG n'est pas résoudre le laplacien,
//  d'où une divergence résiduelle. Ici D\cdot G = -A exactement, donc
//  max|D u| après projection = ||résidu du gradient conjugué||_\infty <= tol.
// ============================================================================

struct ProjectionOptions {
    real density   = 1;
    real tol       = 1e-8;      // divergence maximale tolérée après projection (s^{-1})
    int  maxIter   = 0;         // 0 -> nx*ny (borne théorique de convergence en arithmétique exacte)
    bool warmStart = true;      // part de la pression du pas précédent
};

struct ProjectionStats {
    int  iterations = 0;
    real divBefore  = 0;        // max |D w|
    real divAfter   = 0;        // max |D u|, recalculé sur la grille (pas le résidu de CG)
    bool converged  = false;
};

ProjectionStats project(MacGrid2D& g, real dt, const ProjectionOptions& opt = {},
                        Vector* residualHistory = nullptr);

ProjectionStats project(MacGrid2D& g, const CellLabels& L, real dt,
    const ProjectionOptions& opt = {}, Vector* residualHistory = nullptr);