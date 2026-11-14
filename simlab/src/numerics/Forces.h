#pragma once
#include "grid/MacGrid2D.h"

// ============================================================================
//  Forces extérieures et sources (Stam #2.2 étape 1 et #2.4)
// ============================================================================
//  Étape « force » de Stam :  w1 = w0 + dt f   (Euler explicite, sans risque :
//  f ne dépend pas de la vitesse). Chaque composante de f est évaluée là où
//  vit la composante de vitesse correspondante : f_x aux faces u, f_y aux faces v.
//
//  Les faces de paroi ne sont pas touchées : u·n = 0 y reste vrai.
//
//  Les « pinceaux » (souris, émetteur) ont un profil gaussien
//      w(x) = exp(-|x - c|^2 / r^2),
//  lisse, pour ne pas injecter de hautes fréquences que la grille ne sait pas porter.
// ============================================================================

real brushWeight(Vec2 x, Vec2 center, real radius);

// u <- u + dt · w(x) · accel  sur les faces intérieures (accel en m/s²).
void addBrushForce(MacGrid2D& g, Vec2 center, real radius, Vec2 accel, real dt);

// q <- q + dt · w(x) · rate  (rate en unités de q par seconde) : terme source S de Stam.
void addBrushSource(Field2D& q, Vec2 center, real radius, real rate, real dt);

// Flottabilité (Fedkiw, Stam, Jensen 2001, éq. 8, avec la fumée comme
// traceur de température) : f_y = \beta s, vers le haut si _\beta > 0.
// s est moyennée sur les deux cellules qui encadrent chaque face v intérieure.
void addBuoyancy(MacGrid2D& g, const Field2D& smoke, real beta, real dt);
