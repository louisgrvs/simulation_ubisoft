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
