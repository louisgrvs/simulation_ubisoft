#pragma once
#include <algorithm>
#include <cmath>
#include "grid/MacGrid2D.h"

// Conditions initiales du jour 5.
namespace smoke_scene {

// Tourbillon cellulaire « en boîte », à divergence DISCRÈTE nulle et u·n = 0 :
//   ψ(x, y) = A (Ly/π) sin(πx/Lx) sin(πy/Ly),   u = ∂ψ/∂y,  v = -∂ψ/∂x,
// posé aux coins des cellules comme au jour 4 (max |u| ≈ A).
//
// Pourquoi ce champ : ω = -Δψ est proportionnel à ψ, donc (u·∇)ω = 0. C'est une
// solution STATIONNAIRE d'Euler. Avec ν = 0, la vraie énergie est constante :
// toute perte mesurée est de la dissipation NUMÉRIQUE. Avec ν > 0 (glissement
// parfait), c'est un mode propre de Stokes : E(t) = E(0) exp(-2 ν k² t),
// k² = π²/Lx² + π²/Ly².
inline void setCellularVortex(MacGrid2D& g, real amplitude) {
    const real h = g.h(), Lx = g.size().x, Ly = g.size().y;
    auto psi = [&](int i, int j) {
        return amplitude * (Ly / kPi) * std::sin(kPi * i * h / Lx) * std::sin(kPi * j * h / Ly);
    };
    for (int j = 0; j < g.u.ny(); ++j)
        for (int i = 0; i < g.u.nx(); ++i) g.u(i, j) = (psi(i, j + 1) - psi(i, j)) / h;
    for (int j = 0; j < g.v.ny(); ++j)
        for (int i = 0; i < g.v.nx(); ++i) g.v(i, j) = -(psi(i + 1, j) - psi(i, j)) / h;
}

inline real cellularVortexK2(const MacGrid2D& g) {
    const real Lx = g.size().x, Ly = g.size().y;
    return kPi * kPi / (Lx * Lx) + kPi * kPi / (Ly * Ly);
}

// Disque de fumée lissé sur une épaisseur ~h (évite un bord en marche d'escalier).
inline void addSmokeDisk(Field2D& s, Vec2 center, real radius, real value) {
    for (int j = 0; j < s.ny(); ++j)
        for (int i = 0; i < s.nx(); ++i) {
            const real d = norm(s.position(i, j) - center);
            const real t = std::clamp((radius - d) / s.h() + real(0.5), real(0), real(1));
            s(i, j) += value * t;
        }
}

} // namespace smoke_scene
