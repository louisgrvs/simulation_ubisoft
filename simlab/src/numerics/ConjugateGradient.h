#pragma once
#include <functional>
#include <vector>
#include "core/Vec2.h"

// ============================================================================
//  Gradient conjugué générique, sans matrice
// ============================================================================
//  Résout A x = b pour A symétrique positive (semi-)définie, donnée uniquement
//  par son action y = A x. Sert pour la pression et pour la diffusion
//  implicite : seul l'opérateur change.
//
//  Cas singulier (Neumann pur, noyau = constantes) : avec removeMean = true,
//  on retire la moyenne de b (compatibilité) et de x à la fin (unicité).
//
//  Critère d'arrêt : ||r||_infty <= tol, en norme infinie car pour la pression,
//  r = -(divergence après projection) exactement : tol est donc directement
//  la divergence maximale tolérée, en s^{-1}
// ============================================================================

using Vector         = std::vector<real>;
using LinearOperator = std::function<void(const Vector& x, Vector& y)>;   // y = A x

struct CGOptions {
    real tol        = 1e-8;
    int  maxIter    = 1000;
    bool removeMean = false;            // true pour un problème de Neumann pur
    Vector* history = nullptr;          // si non nul : ||r||_\infty à chaque itération (0 = initial)
};

struct CGResult {
    int  iterations = 0;
    real residual   = 0;                // ||r||_\infty final (résidu de la récurrence)
    bool converged  = false;
};

// x sert d'estimation initiale (« démarrage à chaud ») et reçoit la solution.
CGResult conjugateGradient(const LinearOperator& A, const Vector& b, Vector& x,
                           const CGOptions& opt = {});

// Petits outils d'algèbre linéaire, utiles aussi dans les tests.
real dotProduct(const Vector& a, const Vector& b);
real normInf(const Vector& a);
void subtractMean(Vector& a);
