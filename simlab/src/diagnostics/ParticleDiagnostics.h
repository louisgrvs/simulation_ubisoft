#pragma once
// ============================================================================
//  Grandeurs conservées (ou non) par les transferts — jour 7
// ============================================================================
//  Côté particules :  M = Σ m_p,  P = Σ m_p u_p,  E = ½ Σ m_p |u_p|²,
//                     L = Σ m_p (x_p - c) × u_p        (× : produit vectoriel 2D, scalaire)
//  Côté grille, avec les masses m_f produites par P2G :
//                     P = (Σ_u m_f u_f, Σ_v m_f v_f),  E = ½ Σ m_f u_f²,
//                     L = Σ_u m_f (-(y_f - c_y)) u_f + Σ_v m_f (x_f - c_x) v_f
//  Comparer ces valeurs avant P2G, sur la grille et après G2P localise
//  précisément où se perd l'énergie ou le moment cinétique.
// ============================================================================
#include "numerics/MacTransfer.h"

inline real cross(Vec2 a, Vec2 b) { return a.x * b.y - a.y * b.x; }

real particleMass(const ParticleSet2D& p);
Vec2 particleMomentum(const ParticleSet2D& p);
real particleEnergy(const ParticleSet2D& p);
real particleAngularMomentum(const ParticleSet2D& p, Vec2 center);
real particleMaxSpeed(const ParticleSet2D& p);

Vec2 gridMomentum(const MacGrid2D& g, const MacMasses& m);
real gridEnergy(const MacGrid2D& g, const MacMasses& m);
real gridAngularMomentum(const MacGrid2D& g, const MacMasses& m, Vec2 center);

// Nombre de particules par cellule (champ aux centres), pour l'affichage.
void countParticlesPerCell(const ParticleSet2D& p, Field2D& count);

// ----------------------------------------------------------------------------
//  J9 : énergie potentielle et bruit des particules
// ----------------------------------------------------------------------------
// Énergie potentielle de pesanteur, référence en y = 0 (le fond de la cuve) :
//     E_p = - \sum_p m_p g \cdot x_p        (= \sum_p m_p |g| y_p pour g = (0, -|g|))
// Avec E_c = particleEnergy : E_c + E_p ne peut que décroître dans la rupture de
// barrage (fluide non visqueux, parois fixes). Ce qu'elle perd mesure la dissipation.
real particlePotentialEnergy(const ParticleSet2D& p, Vec2 gravity);

// Bruit : variance des vitesses des particules à l'intérieur de chaque cellule.
//     \bar u_c   = \sum_{p \in c} m_p u_p / M_c
//     \sigma_c^2 = \sum_{p \in c} m_p |u_p - \bar u_c|^2 / M_c          (m²/s²)
// C'est l'énergie que la grille ne voit pas : P2G ne garde (à peu près) que \bar u_c.
// PIC l'efface à chaque pas, FLIP la conserve et la laisse s'accumuler.
// Un champ lisse donne déjà un plancher d'ordre h² |\nabla u|² / 12 (gradient
// à l'intérieur de la cellule) : on compare les schémas entre eux, pas à zéro.
// Les cellules de moins de minCount particules sont ignorées (variance 0).
struct NoiseStats {
    real meanVariance = 0;   // \sum_c M_c \sigma_c^2 / \sum_c M_c : moyenne pondérée par la masse (m²/s²)
    real relative = 0;       // meanVariance / <|u|²>, sans dimension (0 si fluide au repos)
    int  cells = 0;          // cellules prises en compte
};
NoiseStats cellVelocityVariance(const ParticleSet2D& p, Field2D& variance, int minCount = 2);

// ----------------------------------------------------------------------------
//  J11 : forme d'un fluide libre (disque en rotation, simulation complète)
// ----------------------------------------------------------------------------
// Centre de masse :      x_G = \sum_p m_p x_p / \sum_p m_p
// Rayon de giration :    R_g = sqrt( \sum_p m_p |x_p - c|^2 / \sum_p m_p )
// Pour un disque plein uniforme de rayon R, R_g = R / sqrt(2) (≈ 0,212 pour R = 0,3).
// R_g qui grandit : le disque s'étale ; qui diminue : les particules s'agglutinent.
// Le moment cinétique d'un fluide libre se mesure autour de x_G : autour d'un
// point fixe, une dérive du disque ajouterait un terme x_G × P qui n'est pas de la rotation.
Vec2 particleCenterOfMass(const ParticleSet2D& p);   // (0, 0) si aucune particule
real particleRadiusOfGyration(const ParticleSet2D& p, Vec2 center);
