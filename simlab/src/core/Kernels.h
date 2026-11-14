#pragma once
// ============================================================================
//  Formule des noyaux d'interpolation 1D 
// ============================================================================
// Avec la méthode FLIP, la première étape est de construire un champs de véocité interpolé,
// comme dans le pseudo code zhu-bridson. Ici, on définit le noyau de ce champs de vélocité interpolé.
// 
//  Un noyau N(r) donne le poids d'un noeud situé à la distance r·h d'un point.
//  En 2D, le poids est le produit tensoriel  w(x) = N(rx) · N(ry).
//
//  La fonction LinearKernel est le « chapeau » N(r) = max(0, 1 - |r|), support [-1, 1].
//    - partition de l'unité : \sum_i N(g - i) = 1 pour tout g réel ;
//    - reproduit exactement les fonctions affines (interpolation bilinéaire) ;
//    - poids >= 0 : P2G est une moyenne pondérée, G2P une combinaison convexe.
//  C'est le noyau de Zhu & Bridson (#4.2.2, « trilinear » en 3D).
//
//  
// ============================================================================
#include <cmath>
#include "core/Vec2.h"

struct LinearKernel {
    static constexpr int support = 2;   // nombre de noeuds par axe touchés par un point

    static real N(real r) {
        r = std::abs(r);
        return r < 1 ? 1 - r : 0;
    }
    // Dérivée de N (discontinue en 0 et en \pm1 : on choisit la valeur à droite en 0).
    static real dN(real r) {
        if (r <= -1 || r >= 1) return 0;
        return r < 0 ? real(1) : real(-1);
    }
};
