#include "numerics/MacTransfer.h"
#include <cmath>

Stencil2D linearStencil(const Field2D& f, Vec2 x) {
    // Coordonnées indice réelles de x dans ce champ (le décalage n'intervient qu'ici).
    const real gx = x.x / f.h() - f.ox();
    const real gy = x.y / f.h() - f.oy();
    Stencil2D s;
    s.i0 = int(std::floor(gx));
    s.j0 = int(std::floor(gy));
    const real fx = gx - s.i0, fy = gy - s.j0;     // dans [0, 1[
    // N(fx) pour le nœ=oeud i0, N(1 - fx) pour le noeud i0 + 1 (chapeau linéaire).
    s.wx[0] = 1 - fx;
    s.wx[1] = fx;
    s.wy[0] = 1 - fy;
    s.wy[1] = fy;
    return s;
}

namespace {

// Disperse une composante (comp = 0 : x, 1 : y) des vitesses sur un champ de faces.
// C'est le code d'interpolation linéaire de u pour P2G : u_f = \sum_p u_p w_f(x_p)m_p/m_f. 
void scatterComponent(const ParticleSet2D& p, int comp, Field2D& vel, Field2D& mass, int& empty) {
    vel.fill(0);      // accumule d'abord le numérateur \sum w m u
    mass.fill(0);     // puis le dénominateur \sum w m
    for (std::size_t k = 0; k < p.size(); ++k) {
        const Stencil2D s = linearStencil(vel, p.x[k]);
        const real uk = comp == 0 ? p.v[k].x : p.v[k].y;
        for (int b = 0; b < 2; ++b)
            for (int a = 0; a < 2; ++a) {
                const int i = s.i0 + a, j = s.j0 + b;
                if (!inField(vel, i, j)) continue;
                const real wm = s.wx[a] * s.wy[b] * p.m[k];
                mass(i, j) += wm;
                vel(i, j)  += wm * uk;
            }
    }
    empty = 0;
    for (std::size_t n = 0; n < vel.data().size(); ++n) {
        const real m = mass.data()[n];
        if (m > 0) {
            vel.data()[n] /= m;
        } else {
            vel.data()[n] = 0;
            ++empty;
        }
    }
}

// Interpole un champ de faces en x avec le pochoir linéaire.
real gatherComponent(const Field2D& f, Vec2 x) {
    const Stencil2D s = linearStencil(f, x);
    real sum = 0, wsum = 0;
    bool complete = true;
    for (int b = 0; b < 2; ++b)
        for (int a = 0; a < 2; ++a) {
            const int i = s.i0 + a, j = s.j0 + b;
            if (!inField(f, i, j)) { complete = false; continue; }
            const real w = s.wx[a] * s.wy[b];
            sum  += w * f(i, j);
            wsum += w;
        }
    // Pochoir complet : poids exactement transposés de P2G, sans division.
    // Pochoir tronqué (près d'une paroi) : on renormalise.
    if (complete) return sum;
    return wsum > 0 ? sum / wsum : real(0);
}

} // namespace

//le code pour P2G :
void particlesToMac(const ParticleSet2D& particles, MacGrid2D& grid, MacMasses& masses) {
    if (!masses.u.sameLayout(grid.u)) masses.u = grid.u;
    if (!masses.v.sameLayout(grid.v)) masses.v = grid.v;
    scatterComponent(particles, 0, grid.u, masses.u, masses.emptyU); //Donc là on traite u
    scatterComponent(particles, 1, grid.v, masses.v, masses.emptyV); //Là on traite v
}

Vec2 interpolateMac(const MacGrid2D& grid, Vec2 x) {
    return {gatherComponent(grid.u, x), gatherComponent(grid.v, x)};
}

//là on acutalise la particle selon G2P.
void macToParticlesPIC(const MacGrid2D& grid, ParticleSet2D& particles) {
    for (std::size_t k = 0; k < particles.size(); ++k) particles.v[k] = interpolateMac(grid, particles.x[k]);
}
