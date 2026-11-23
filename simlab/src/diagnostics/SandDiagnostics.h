#pragma once
// ============================================================================
//  Mesures sur un dépôt de sable (ou une flaque) — jour 13
// ============================================================================
//  Profil de hauteur : dans chaque colonne de cellules, la hauteur de la plus haute
//  particule (0 si la colonne est vide).
//  Angle d'un flanc : régression linéaire h(x) sur les colonnes du flanc dont la
//  hauteur est entre lo·H et hi·H (H : hauteur maximale), en partant du sommet ;
//  on s'arrête à la première colonne plus basse que lo·H.
//  Angle du tas : moyenne des deux flancs (ou le seul flanc mesurable). Mesure locale,
//  bruitée quand le tas ne fait que quelques cellules de haut.
//  Angle équivalent : celui du triangle de même hauteur H et de même base W, où W est la
//  largeur de la zone continue autour du sommet où h ≥ 0,1·H : tan θ_eq = 2 H / W.
//  Plus robuste (deux grandeurs globales), mais il sous-estime la pente d'un tas tronqué.
// ============================================================================
#include <vector>

#include "core/Vec2.h"
#include "particles/ParticleSet2D.h"

struct PileShape {
    real height = 0;          // hauteur maximale du dépôt (m)
    real xMin = 0, xMax = 0;  // extension horizontale des particules (m)
    real spread = 0;          // largeur du dépôt, xMax − xMin (m)
    real peakX = 0;           // abscisse du sommet (m)
    real angleLeftDeg = 0;    // flanc gauche (°), 0 si non mesurable
    real angleRightDeg = 0;   // flanc droit (°)
    real angleDeg = 0;        // moyenne des flancs mesurés
    int  flanks = 0;          // nombre de flancs mesurés (0, 1 ou 2)
    real baseWidth = 0;       // W : largeur continue autour du sommet où h ≥ 0,1·H (m)
    real angleEquivalentDeg = 0;   // atan(2 H / W) (°)
};

// Hauteur de la plus haute particule dans chacune des nx colonnes de largeur h.
std::vector<real> heightProfile(const ParticleSet2D& p, int nx, real h);

// Il faut au moins 3 colonnes dans un flanc pour mesurer sa pente.
// lo = 0,3 : on écarte la fine couche d'étalement au pied du tas, qui aplatirait la pente.
PileShape measurePile(const ParticleSet2D& p, int nx, real h, real lo = real(0.3), real hi = real(0.8));
