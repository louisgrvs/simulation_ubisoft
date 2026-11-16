#pragma once
// ============================================================================
//  Extrapolation de vitesse couche par couche (Zhu-Bridson §4.2.3, version simplifiée)
// ============================================================================
//  Zhu-Bridson utilisent le balayage rapide (fast sweeping). On le remplace par :
//
//    pour chaque couche :
//      toute face INVALIDE ayant au moins un voisin VALIDE (4 voisins, même grille)
//      reçoit la MOYENNE de ses voisins valides ;
//      ces faces ne deviennent valides qu'à la FIN de la couche (style Jacobi),
//      pour que le résultat ne dépende pas de l'ordre de parcours.
//
//  Utilisation dans FlipSolver::step :
//    étape 3 (après P2G)        : valide = face ayant reçu un poids > 0
//    étape 8 (après projection) : valide = face fluide (uFaceIsFluid / vFaceIsFluid)
//
//  Nombre de couches : au moins la distance parcourue par une particule en un pas,
//  plus 2 de marge (point milieu RK2 + pochoir de G2P) :
//      layers = max(2, ceil(vmax * dt / h) + 2)
// ============================================================================
#include <cstdint>
#include <vector>
#include "grid/Field2D.h"
#include "grid/MacGrid2D.h"
#include "grid/CellLabels.h"

// Un octet par face, rangé comme le Field2D correspondant (indice i + f.nx() * j).
// 1 = valeur connue, 0 = à extrapoler.
using FaceMask = std::vector<std::uint8_t>;

// Extrapole f sur au plus `layers` couches à partir des faces valides.
// Met à jour `valid`. Renvoie le nombre de couches réellement faites.
int extrapolateLayers(Field2D& f, FaceMask& valid, int layers);

// Applique extrapolateLayers à g.u puis à g.v.
void extrapolateVelocity(MacGrid2D& g, FaceMask& uValid, FaceMask& vValid, int layers);

// Masque de l'étape 3 : valide si le poids de P2G sur la face dépasse eps.
FaceMask maskFromWeights(const Field2D& weights, real eps = 0);

// Masques de l'étape 8 : valide si la face est fluide (une cellule Fluid, aucune Solid).
FaceMask uMaskFromLabels(const CellLabels& L);   // taille (nx+1) × ny
FaceMask vMaskFromLabels(const CellLabels& L);   // taille nx × (ny+1)

// Nombre de couches conseillé pour un pas dt.
int extrapolationLayerCount(real vmax, real dt, real h);