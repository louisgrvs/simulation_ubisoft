#pragma once
#include <vector>
#include "grid/MacGrid2D.h"

// ============================================================================
//  Opérateurs discrets sur la grille MAC (domaine fermé par 4 parois)
// ============================================================================
//
//  Notations : D = divergence (faces -> cellules), G = gradient (cellules -> faces),
//              A = -D G  (« moins le laplacien », symétrique positif semi-défini).
//
//  Condition de paroi u·n = 0 :
//      u(0, j) = u(nx, j) = 0      (parois gauche et droite)
//      v(i, 0) = v(i, ny) = 0      (parois bas et haut)
//  Ces faces ne sont JAMAIS modifiées par la projection : G est nul sur elles.
//
//  Conséquences exactes (à la précision machine près) :
//   (1) D G = -A  : la projection et le solveur parlent du même opérateur
//                  -> divergence après projection = résidu du gradient conjugué.
//   (2) <D u, p>_cellules = -<u, G p>_faces  : D et -G sont adjoints
//                  (intégration par parties discrète, terme de bord nul).
//   (3) A 1 = 0    : les constantes forment le noyau (Neumann pur),
//                  et sum(D u) = 0 dès que u·n = 0 -> le système est compatible.
// ============================================================================

// Impose u·n = 0 sur les 4 parois (glissement parfait : seule la composante normale est fixée).
void enforceWallBC(MacGrid2D& g);

// div(i,j) = (u(i+1,j) - u(i,j) + v(i,j+1) - v(i,j)) / h      (div : champ de cellules)
void computeDivergence(const MacGrid2D& g, Field2D& div);

// max |div| sur toutes les cellules : la mesure de qualité de la projection.
real maxAbsDivergence(const MacGrid2D& g);

// Gradient d'un champ de cellules q, rangé sur les faces :
//   gx(i,j) = (q(i,j) - q(i-1,j)) / h   pour 1 <= i <= nx-1, et 0 sur les parois
//   gy(i,j) = (q(i,j) - q(i,j-1)) / h   pour 1 <= j <= ny-1, et 0 sur les parois
// gx a le layout de u, gy celui de v.
void computeGradient(const Field2D& q, Field2D& gx, Field2D& gy);

// u <- u - scale * G q  (faces intérieures seulement ; les parois gardent u·n = 0).
void subtractGradient(MacGrid2D& g, const Field2D& q, real scale);

// Laplacien sans sans : y = A x = -Delta_h x, x et y rangés comme un champ de cellules
// (indice k = i + nx*j). Condition de Neumann homogène aux parois : un voisin absent
// est simplement omis (pression fantôme = pression de la cellule, flux nul).
//   (A x)_ij = ( n_ij x_ij - sum_{voisins présents} x_voisin ) / h^2,   n_ij in {2, 3, 4}
// Foncteur pour pouvoir être passé tel quel au gradient conjugué.
struct NegLaplacianMAC {
    int nx, ny;
    real h;
    void operator()(const std::vector<real>& x, std::vector<real>& y) const;
};
