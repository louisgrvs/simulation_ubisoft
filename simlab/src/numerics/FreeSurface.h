#pragma once
// ============================================================================
//  Surface libre : étiquettes et conditions aux parois — jour 8
// ============================================================================
#include "grid/CellLabels.h"
#include "grid/MacGrid2D.h"
#include "particles/ParticleSet2D.h"
#include "geometry/SolidSdf2D.h"

// Étiquettes Fluid / Air / Solid à partir des particules et des solides (Zhu-Bridson §4.2.2).
// Les poches d'air isolées d'au plus maxBubbleCells cellules sont comblées (→ Fluid) :
// sans cela, elles implosent sous la pression et le volume d'eau fond (0 : désactivé).
void labelCells(CellLabels& L, const ParticleSet2D& p, real h, const SolidSdf2D* solids,
                int maxBubbleCells = 8);

// Composante normale nulle sur toute face qui touche une cellule Solid (murs compris).
void enforceSolidBC(MacGrid2D& g, const CellLabels& L);

// Une face est fluide si une de ses deux cellules est Fluid et aucune n'est Solid.
// Face u(i, j) : entre les cellules (i−1, j) et (i, j). Face v(i, j) : entre (i, j−1) et (i, j).
bool uFaceIsFluid(const CellLabels& L, int i, int j);
bool vFaceIsFluid(const CellLabels& L, int i, int j);
