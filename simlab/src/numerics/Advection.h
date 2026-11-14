#pragma once
#include "grid/MacGrid2D.h"

// Advection semi-lagrangienne (Stam 1999, #2.2 « méthode des caractéristiques »).
// On résout  dq/dt + u·grad(q) = 0  : q est constant le long des trajectoires.
// Pour chaque échantillon x, on remonte la trajectoire sur dt pour trouver le
// point de départ x_back, puis  q_new(x) = q_old(x_back).

enum class Integrator {
    Euler,  // x_back = x - dt u(x)                         ordre 1 (spirale vers l'extérieur en rotation)
    RK2     // x_mid  = x - dt/2 u(x) ; x_back = x - dt u(x_mid)   ordre 2 (point milieu)
};

// Remonte la caractéristique passant par x sur un temps dt.
// Le point obtenu est borné au domaine [0, nx h] x [0, ny h].
Vec2 traceBack(const MacGrid2D& grid, Vec2 x, real dt, Integrator scheme = Integrator::RK2);

// qOut(x_ij) = q( traceBack(x_ij) ) pour chaque échantillon (i, j) de q.
// Fonctionne pour tout Field2D (centres, faces…) : on utilise q.position(i, j).
// q et qOut doivent avoir la même disposition et être DEUX tableaux distincts
// (sinon on lirait des valeurs déjà mises à jour).
void advectScalar(const MacGrid2D& grid, const Field2D& q, Field2D& qOut, real dt,
                  Integrator scheme = Integrator::RK2);

// Pas de temps tel que  maxSpeed * dt / h = cfl.
// Le semi-lagrangien est stable pour tout CFL ; le CFL règle ici la précision.
real cflTimeStep(const MacGrid2D& grid, real cfl);
