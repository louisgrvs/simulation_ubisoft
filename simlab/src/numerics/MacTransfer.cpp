#include "numerics/MacTransfer.h"
#include <cassert>
#include <cmath>
#include "core/Vec2.h"

Stencil2D linearStencil(const Field2D& f, Vec2 x) {
    // Coordonnées indice réelles de x dans ce champ (le décalage n'intervient qu'ici).
    const real gx = x.x / f.h() - f.ox();
    const real gy = x.y / f.h() - f.oy();
    Stencil2D s;
    s.i0 = int(std::floor(gx));
    s.j0 = int(std::floor(gy));
    const real fx = gx - s.i0, fy = gy - s.j0;     // dans [0, 1[
    // N(fx) pour le noeud i0, N(1 - fx) pour le noeud i0 + 1 (chapeau linéaire).
    s.wx[0] = 1 - fx;
    s.wx[1] = fx;
    s.wy[0] = 1 - fy;
    s.wy[1] = fy;
    // Dérivées par rapport à la position : d fx / dx = 1/h, donc -1/h et +1/h.
    // Le gradient du noeud (a,b) vaut (dwx[a]·wy[b], wx[a]·dwy[b]).
    s.dwx[0] = -1 / f.h();
    s.dwx[1] = 1 / f.h();
    s.dwy[0] = -1 / f.h();
    s.dwy[1] = 1 / f.h();
    return s;
}

namespace {

    // Disperse une composante (comp = 0 : x, 1 : y) des vitesses sur un champ de faces.
    // PIC  : u_f = \sum_p w_f(x_p) m_p u_p / m_f
    // APIC : u_f = \sum_p w_f(x_p) m_p (u_p + c_p · (x_f - x_p)) / m_f   (affine = true)
    void scatterComponent(const ParticleSet2D& p, int comp, bool affine,
        Field2D& vel, Field2D& mass, int& empty) {
        vel.fill(0);      // accumule d'abord le numérateur \sum w m (u + c·(x_f - x_p))
        mass.fill(0);     // puis le dénominateur \sum w m
        for (std::size_t k = 0; k < p.size(); ++k) {
            const Stencil2D s = linearStencil(vel, p.x[k]);
            const real uk = comp == 0 ? p.v[k].x : p.v[k].y;
            const Vec2 cap = affine ? (comp == 0 ? p.cx[k] : p.cy[k]) : Vec2{ 0, 0 };
            const Vec2 xp = p.x[k];
            for (int b = 0; b < 2; ++b)
                for (int a = 0; a < 2; ++a) {
                    const int i = s.i0 + a, j = s.j0 + b;
                    if (!inField(vel, i, j)) continue;
                    const real wm = s.wx[a] * s.wy[b] * p.m[k];
                    const Vec2 xf = vel.position(i, j);   // position de LA face (décalage du champ)
                    mass(i, j) += wm;
                    vel(i, j) += wm * (uk + dot(cap, xf - xp));
                }
        }
        empty = 0;
        for (std::size_t n = 0; n < vel.data().size(); ++n) {
            const real m = mass.data()[n];
            if (m > 0) {
                vel.data()[n] /= m;
            }
            else {
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
                sum += w * f(i, j);
                wsum += w;
            }
        // Pochoir complet : poids exactement transposés de P2G, sans division.
        // Pochoir tronqué (près d'une paroi) : on renormalise.
        if (complete) return sum;
        return wsum > 0 ? sum / wsum : real(0);
    }

    // Gradient en x du champ bilinéaire interpolé : \sum_f \nabla w_f(x) f_f.
    // C'est la partie affine c_p d'APIC (Jiang et al. #5.3, #6).
    Vec2 gatherGradient(const Field2D& f, Vec2 x) {
        const Stencil2D s = linearStencil(f, x);
        const bool xCut = s.i0 < 0 || s.i0 + 1 >= f.nx();   // une colonne du pochoir manque
        const bool yCut = s.j0 < 0 || s.j0 + 1 >= f.ny();   // une ligne du pochoir manque
        Vec2 grad = { 0, 0 };
        real wxPresent = 0, wyPresent = 0;   // somme des wx (colonnes) et wy (lignes) présents
        for (int a = 0; a < 2; ++a)
            if (s.i0 + a >= 0 && s.i0 + a < f.nx()) wxPresent += s.wx[a];
        for (int b = 0; b < 2; ++b)
            if (s.j0 + b >= 0 && s.j0 + b < f.ny()) wyPresent += s.wy[b];

        for (int b = 0; b < 2; ++b)
            for (int a = 0; a < 2; ++a) {
                const int i = s.i0 + a, j = s.j0 + b;
                if (!inField(f, i, j)) continue;
                const Vec2 nablaw = { s.dwx[a] * s.wy[b], s.wx[a] * s.dwy[b] };   // règle du produit
                grad += f(i, j) * nablaw;
            }
        // Pochoir complet : rien à corriger (wxPresent = wyPresent = 1).
        // Pochoir tronqué : la vitesse est prolongée par une constante selon l'axe
        // tronqué, donc la dérivée selon cet axe vaut 0 ; selon l'autre axe, on
        // divise par le poids de la ligne (ou colonne) présente, ce qui donne la
        // différence finie (f(i0+1) - f(i0)) / h sur cette ligne.
        if (yCut) {
            grad.y = 0;
            grad.x = wyPresent > 0 ? grad.x / wyPresent : real(0);
        }
        if (xCut) {
            grad.x = 0;
            grad.y = wxPresent > 0 ? grad.y / wxPresent : real(0);
        }
        return grad;
    }

} // namespace

// P2G : PIC/FLIP (affine = false) ou APIC (affine = true).
void particlesToMac(const ParticleSet2D& particles, MacGrid2D& grid, MacMasses& masses, bool affine) {
    if (!masses.u.sameLayout(grid.u)) masses.u = grid.u;
    if (!masses.v.sameLayout(grid.v)) masses.v = grid.v;
    scatterComponent(particles, 0, affine, grid.u, masses.u, masses.emptyU); // composante u, avec cx
    scatterComponent(particles, 1, affine, grid.v, masses.v, masses.emptyV); // composante v, avec cy
}

Vec2 interpolateMac(const MacGrid2D& grid, Vec2 x) {
    return { gatherComponent(grid.u, x), gatherComponent(grid.v, x) }; // v_ap^{n+1} = \sum_i w_aip^n \tilde{v}_ai^{n+1}
}

// G2P PIC : la vitesse de la particule est remplacée par celle de la grille.
void macToParticlesPIC(const MacGrid2D& grid, ParticleSet2D& particles) {
    for (std::size_t k = 0; k < particles.size(); ++k) particles.v[k] = interpolateMac(grid, particles.x[k]);
}

// G2P FLIP / mélange PIC-FLIP (#4.2.4). Le flipratio c'est la constante alpha du mélange des deux.
void macToParticlesFLIP(const MacGrid2D& grid, const Field2D& uOld, const Field2D& vOld,
    ParticleSet2D& particles, real flipRatio) {
    assert(uOld.sameLayout(grid.u) && vOld.sameLayout(grid.v));
    // Variation de la grille pendant le pas, face par face : Delta u = u_new - u_old.
    Field2D du = grid.u, dv = grid.v;
    for (std::size_t n = 0; n < du.data().size(); ++n) du.data()[n] -= uOld.data()[n];
    for (std::size_t n = 0; n < dv.data().size(); ++n) dv.data()[n] -= vOld.data()[n];

    const real a = flipRatio;
    for (std::size_t k = 0; k < particles.size(); ++k) {
        const Vec2 x = particles.x[k];
        const Vec2 pic = interpolateMac(grid, x);                               // u_p^PIC
        const Vec2 flip = particles.v[k] + Vec2{ gatherComponent(du, x), gatherComponent(dv, x) };
        particles.v[k] = a * flip + (1 - a) * pic;
    }
}

// G2P APIC (Jiang et al. 2015, #6) : vitesse comme PIC, plus c_a,p = \sum_i \nabla w_aip \tilde{v}_ai.
void macToParticlesAPIC(const MacGrid2D& grid, ParticleSet2D& particles) {
    for (std::size_t k = 0; k < particles.size(); ++k) {
        const Vec2 x = particles.x[k];
        particles.v[k] = interpolateMac(grid, x);      // même que pour PIC
        particles.cx[k] = gatherGradient(grid.u, x);   // (du/dx, du/dy)
        particles.cy[k] = gatherGradient(grid.v, x);   // (dv/dx, dv/dy)
    }
}