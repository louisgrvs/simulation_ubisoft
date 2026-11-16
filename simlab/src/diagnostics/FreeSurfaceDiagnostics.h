#pragma once
// ============================================================================
//  Diagnostics de la surface libre — jour 8
// ============================================================================
//  Les trois vérifications du planning :
//    1. cuve au repos : vitesses parasites très faibles      → maxSpeed, rmsSpeed
//    2. pression proche de ρgh                                → hydroErr
//    3. volume de fluide stable                               → fluidArea, surfaceHeight
//
//  Pression de référence. Le schéma impose p = 0 au CENTRE de la première cellule
//  d'air, pas à la vraie surface. À l'équilibre, la solution exacte du système
//  discret est donc :
//      p_ref(i, j) = ρ g (y_s − y_j),   y_s = centre de la première cellule d'air
//                                             au-dessus de (i, j) dans la colonne i.
//  La surface effective est une demi-cellule au-dessus de la vraie (erreur d'ordre 1,
//  que corrigerait la méthode du fluide fantôme). Comparer à ρg(H − y) donnerait
//  un écart de ρgh/2 qui n'est pas un bogue.
//
//  Ce fichier ne dépend que de la grille, des étiquettes et des particules, pas
//  du solveur : il resservira pour la rupture de barrage (J9) et le sable (J13).
// ============================================================================
#include <ostream>
#include "grid/CellLabels.h"
#include "grid/MacGrid2D.h"
#include "particles/ParticleSet2D.h"

struct FreeSurfaceStats {
    int  fluidCells = 0, airCells = 0, solidCells = 0;
    real fluidArea = 0;        // fluidCells · h² (le « volume » en 2D, m²)
    real surfaceHeight = 0;    // hauteur moyenne de la colonne d'eau (m), sur les colonnes non vides
    real maxSpeed = 0;         // max |u_p| sur les particules (m/s)
    real rmsSpeed = 0;         // sqrt(moyenne de |u_p|²) (m/s)
    real hydroErr = 0;         // max |p − p_ref| / max p_ref, sur les cellules fluides
    real divAfter = 0;         // max |div u| sur les cellules fluides (s⁻¹)
};

// Pression hydrostatique discrète de référence dans la cellule (i, j) (supposée Fluid).
// Renvoie false si la colonne n'a pas de cellule d'air au-dessus (pas de référence).
bool hydrostaticReference(const CellLabels& L, real h, int i, int j, real density, real g, real& pRef);

// Toutes les mesures d'un coup. g = norme de la gravité (m/s²).
FreeSurfaceStats measureFreeSurface(const MacGrid2D& grid, const CellLabels& L, const ParticleSet2D& p,
                                    real density, real g);

// Export CSV : une ligne par instant (en-tête à écrire une fois).
void writeFreeSurfaceCsvHeader(std::ostream& os);
void writeFreeSurfaceCsvRow(std::ostream& os, real t, const FreeSurfaceStats& s);

// Profil de pression dans la colonne i : y, p, p_ref (cellules fluides seulement).
void writePressureProfileCsv(std::ostream& os, const MacGrid2D& grid, const CellLabels& L, int i,
                             real density, real g);
