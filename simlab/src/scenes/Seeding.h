#pragma once
// ============================================================================
//  Initialisation des champs de vélocités (première étapes du pseudo-algo) des particules (Zhu & Bridson #4.2.1)
// ============================================================================
//  L'article met 8 particules par cellule en 3D (2 par axe), perturbées au
//  hasard ; en 2D cela donne 4 particules par cellule.
//
//  Ensemencement stratifié : chaque cellule est découpée en k x k sous-cellules
//  (k = perAxis) et l'on place un point par sous-cellule 
//      x = (i + (a + \xi) / k) h,   a \in {0...k-1},  \xi = 1/2+ jitter·(U - 1/2),  U \in [0, 1[
//  jitter = 0 : réseau régulier (centres des sous-cellules), en gros on fait un tirage uniforme ;
//  jitter = 1 : point uniforme dans toute la sous-cellule, donc la on fait un réseau stratifié.
//  Chaque particule reste donc dans sa sous-cellule : la densité est uniforme
//  à l'échelle h, sans les amas d'un tirage uniforme global.
//
//  Une sous-cellule est retenue si son centre est dans la région : le résultat
//  ne dépend pas du tirage, et le bord de la région est résolu à h/k près.
//
//  Reproductibilité : on tire les nombres directement de std::mt19937 (suite
//  fixée par la norme C++), et non via std::uniform_real_distribution, dont
//  l'algorithme varie d'une bibliothèque standard à l'autre (MSVC != GCC).
//  Même graine => mêmes particules sur toutes les machines.
//
//  Masse : m_p = ρ h^2/ k^2, donc la masse totale vaut \rho × (aire ensemencée).
// ============================================================================
#include <functional>
#include "particles/ParticleSet2D.h"

// Prédicat "ce point est dans la région ".
using Region = std::function<bool(Vec2)>;

namespace region {
Region box(Vec2 lo, Vec2 hi);              // rectangle [lo, hi[
Region disk(Vec2 center, real radius);     // disque ouvert
Region everywhere();
} // namespace region

struct SeedOptions {
    int      perAxis = 2;      // k : k² particules par cellule (2 -> 4, comme l'article)
    real     jitter  = 1;      // 0 : réseau régulier, 1 : uniforme dans la sous-cellule
    unsigned seed    = 1;      // graine de std::mt19937
    real     density = 1;      // ρ (kg/m²)
};

// Ajoute (sans effacer) des particules dans la région, sur une grille nx × ny
// de pas h. Vitesse initiale donnée par velocity(x) (nulle si absente).
// Renvoie le nombre de particules ajoutées.
int seedParticles(ParticleSet2D& particles, int nx, int ny, real h, const Region& inside,
                  const SeedOptions& opt = {},
                  const std::function<Vec2(Vec2)>& velocity = nullptr);
