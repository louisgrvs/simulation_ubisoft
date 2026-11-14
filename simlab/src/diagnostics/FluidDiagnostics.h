#pragma once
#include "grid/MacGrid2D.h"

// Grandeurs de validation, toutes en unités physiques.

// Énergie cinétique discrète (masse volumique 1) :
//   E = 1/2 h^2 ( \sum _faces u^2 + \sum _faces v^2 )
// C'est la norme pour laquelle la projection est une projection orthogonale
//  la projection ne peut donc jamais l'augmenter.
real kineticEnergy(const MacGrid2D& g);

// Quantité totale d'un scalaire aux centres : h^2 \sum  q (« masse » de fumée).
real totalAmount(const Field2D& q);

// Norme de la vitesse au centre des cellules (moyenne des deux faces de chaque direction).
void computeSpeed(const MacGrid2D& g, Field2D& out);

// Vorticité \omega  = \partial v/\partial x - \partial u/\partial y (CONVENTIONS.md : positive dans le sens trigonométrique).
// Calculée naturellement aux coins (i h, j h) de la grille MAC, puis moyennée
// sur les 4 coins de chaque cellule. Aux parois (glissement parfait), \omega  = 0.
void computeVorticity(const MacGrid2D& g, Field2D& out);
