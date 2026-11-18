#pragma once
// ============================================================================
//  Transferts particules <-> grille MAC, schéma PIC (Zhu & Bridson #4.2.2, #4.2.4)
// ============================================================================
//
//  P2G (particules -> grille), composante par composante, chacune sur ses faces :
//
//      m_f =  \sum_p w_f(x_p) m_p                       (masse portée par la face f)
//      u_f = \sum_p w_f(x_p) m_p u_p  /  m_f            (moyenne pondérée)
//
//  avec w_f(x) = N((x - x_f)/h) N((y - y_f)/h), N le chapeau de LinearKernel (N(r) = max(0,1-r)).
//  Une face sans aucune particule (m_f = 0) est vide : vitesse 0.
//
//  G2P (grille -> particules), schéma PIC : la vitesse de la particule est
//  remplacée par l'interpolation bilinéaire de la grille, avec les mêmes poids :
//
//      u_p = \sum_f w_f(x_p) u_f
//
//  P2G et G2P utilisent tous deux linearStencil(). Les poids de
//  G2P sont exactement les transposés de ceux de P2G, ce qui donne :
//    - conservation de la quantité de mouvement par chaque transfert ;
//    - E_particules >= E_grille >= E_particules après retour (Cauchy-Schwarz),
//      c'est-à-dire la dissipation de PIC.
//
//  Aux parois, le pochoir d'une particule peut sortir du champ (par exemple
//  j0 = -1 pour u dans la demi-cellule du bas). Convention (CONVENTIONS.md) :
//    - P2G ignore les noeuds hors champ (la division par m_f compense) ;
//    - G2P renormalise par la somme des poids des nœuds présents, ce qui
//      revient à prolonger la composante tangentielle par une constante
//      (glissement parfait).
//  Loin des parois, aucun des deux cas ne se produit.
// ============================================================================
#include "core/Kernels.h"
#include "grid/MacGrid2D.h"
#include "particles/ParticleSet2D.h"

// Pochoir bilinéaire 2 × 2 d'un point x sur un champ de décalage quelconque.
// Le noeud (i0 + a, j0 + b), a, b \in {0, 1}, a le poids wx[a] · wy[b].
// Les indices peuvent sortir du champ : tester avec inField().
struct Stencil2D {
    int  i0 = 0, j0 = 0;
    real wx[2] = { 0, 0 };
    real wy[2] = { 0, 0 };
};

Stencil2D linearStencil(const Field2D& f, Vec2 x); // Coordonnées indice réelles de x dans ce champ (le décalage n'intervient qu'ici).

inline bool inField(const Field2D& f, int i, int j) {
    return i >= 0 && i < f.nx() && j >= 0 && j < f.ny();
}

// Masses de grille produites par P2G. Mêmes dispositions que grid.u et grid.v.
// Elles servent aux diagnostics (énergie, quantité de mouvement de la grille),
// à repérer les faces vides, à FLIP et à l'extrapolation.
struct MacMasses {
    Field2D u, v;
    int emptyU = 0, emptyV = 0;     // nombre de faces sans masse
};

// P2G : écrase grid.u et grid.v (pas grid.p, gardée pour le démarrage à chaud).
// N'impose PAS u\dot n = 0 : c'est le rôle du solveur (enforceWallBC / projection).
void particlesToMac(const ParticleSet2D& particles, MacGrid2D& grid, MacMasses& masses);

// Vitesse de grille interpolée en x, avec les poids de G2P.
Vec2 interpolateMac(const MacGrid2D& grid, Vec2 x);

// G2P, schéma PIC : u_p <- \sum_f w_f(x_p) u_f.
void macToParticlesPIC(const MacGrid2D& grid, ParticleSet2D& particles);

// ----------------------------------------------------------------------------
//  G2P, schéma FLIP et mélange PIC/FLIP (Zhu & Bridson #4.2.4)
// ----------------------------------------------------------------------------
//  PIC remplace la vitesse de la particule par celle de la grille :
//      u_p^PIC  = \sum_f w_f(x_p) u_f^new
//  FLIP n'interpole que la VARIATION de la grille pendant le pas (forces + projection) :
//      u_p^FLIP = u_p + \sum_f w_f(x_p) (u_f^new - u_f^old)
//  Mélange réglable, alpha = flipRatio \in [0, 1] :
//      u_p <- alpha u_p^FLIP + (1 - alpha) u_p^PIC
//  alpha = 0 : PIC exactement (dissipatif, lisse) ; alpha = 1 : FLIP pur (peu
//  dissipatif, mais bruité : rien ne rappelle les écarts entre particules d'une
//  même cellule, que la grille ne voit pas).
//
//  uOld, vOld : vitesse de grille juste après P2G (et l'extrapolation n°1),
//  AVANT les forces. Mêmes dispositions que grid.u et grid.v.
//  Mêmes poids (linearStencil) que P2G et PIC, même traitement des parois.
void macToParticlesFLIP(const MacGrid2D& grid, const Field2D& uOld, const Field2D& vOld,
    ParticleSet2D& particles, real flipRatio);