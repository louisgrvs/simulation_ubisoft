#pragma once
#include "grid/MacGrid2D.h"
#include "particles/ParticleSet2D.h"
#include "geometry/SolidSdf2D.h"

// ============================================================================
//  Déplacement des particules (Zhu-Bridson #4.2.5)
// ============================================================================
//  RK2 (point milieu) dans le champ de vitesse de la grille, déjà extrapolé.
//  Le pas dt est découpé en sous-pas pour qu'aucune particule ne parcoure plus
//  de cflSub cellules par sous-pas :  dt_s <= cflSub · h / vmax.
//  Après chaque sous-pas, les particules entrées dans un solide en sont repoussées.
// ============================================================================

struct AdvectionStats {
    int  substeps = 0;
    real vmax = 0;   // vitesse maximale sur la grille utilisée pour la CFL
};

AdvectionStats advectParticlesRK2(ParticleSet2D& p, const MacGrid2D& g, real dt,
    real cflSub, const SolidSdf2D& solids, real margin);