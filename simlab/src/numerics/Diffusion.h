#pragma once
#include "grid/MacGrid2D.h"
#include "numerics/ConjugateGradient.h"

// ============================================================================
//  Diffusion implicite et dissipation (Stam #2.2 étape 3, #2.4, #3.2)
// ============================================================================
//
//  On résout  \partial q/\partial t = \kappa  \Delta q  par Euler implicite, comme Stam :
//
//      (I - \kappa  dt \Delta _h) q_new = q_old
//
//  Le schéma explicite q_new = q_old + \kappa  dt \Delta _h q_old n'est stable que si
//  \kappa  dt / h^2 <= 1/4 ; l'implicite l'est pour tout dt (c'est tout l'esprit de
//  « Stable Fluids »). Chaque mode de Fourier discret de valeur propre \lambda ≥ 0
//  de -\Delta _h est simplement multiplié par 1 / (1 + \kappa  dt \lambda ) \in ]0, 1].
//
//  L'opérateur I - \kappa  dt \Delta _h est symétrique défini positif : on réutilise tel
//  quel le gradient conjugué du jour 4, seul l'opérateur change.
//
//  Conditions aux parois, voisin par voisin (sans matrice, sans cellule fantôme) :
//   - Neumann    : le voisin hors domaine est omis (flux nul à travers la paroi) ;
//   - Dirichlet0 : le voisin hors domaine est une valeur fixée à 0, à distance h.
//
//  Pour la vitesse sur la grille MAC (glissement parfait, « free slip ») :
//   - composante NORMALE (u en x, v en y) : Dirichlet 0, c'est u.n = 0 ;
//   - composante TANGENTIELLE (u en y, v en x) : Neumann, pas de frottement.
//  Pour un scalaire aux centres (fumée) : Neumann partout -> quantité conservée.
// ============================================================================

enum class EdgeCondition { Neumann, Dirichlet0 };

// y = (I - c \Delta _h) x  sur un tableau nx × ny d'inconnues (indice k = i + nx j),
// avec c = \kappa  dt / h² (sans dimension).
//   y_k = x_k + c ( n_k x_k - \sum_{voisins présents} x_voisin )
// où n_k compte les voisins présents, plus les voisins de Dirichlet hors domaine.
struct HelmholtzOperator {
    int  nx = 0, ny = 0;
    real c  = 0;
    EdgeCondition xEdges = EdgeCondition::Neumann;   // bords gauche et droit
    EdgeCondition yEdges = EdgeCondition::Neumann;   // bords bas et haut
    void operator()(const Vector& x, Vector& y) const;
};

struct DiffusionStats {
    int  iterations = 0;     // total des itérations de CG (u + v pour la vitesse)
    real residual   = 0;     // pire résidu final, relatif à max|q_old|
    bool converged  = true;
};

// Scalaire aux centres des cellules, Neumann sur les 4 parois.
// tol est RELATIVE : ||r||_\infty <= tol · max|q_old|.
DiffusionStats diffuseCells(Field2D& q, real kappa, real dt, real tol = 1e-10);

// Vitesse MAC : u et v diffusés séparément (les deux composantes ne se couplent pas).
// Les faces de paroi restent à 0 et ne sont pas des inconnues.
DiffusionStats diffuseVelocity(MacGrid2D& g, real nu, real dt, real tol = 1e-10);

// Dissipation (Stam #2.4, terme -\alpha a), implicite comme dans son Sstep :
//   q <- q / (1 + α dt)
void dissipate(Field2D& q, real alpha, real dt);
