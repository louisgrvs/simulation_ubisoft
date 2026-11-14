#pragma once
#include <algorithm>
#include <cmath>
#include "grid/MacGrid2D.h"

namespace smoke_scene {

// Tourbillon cellulaire « en boîte », à divergence discrète nulle et u.n = 0 :
//   \psi(x, y) = A (Ly/\pi ) sin(\pi x/Lx) sin(\pi y/Ly),   u = \partial \psi/\partial y,  v = -\partial \psi/\partial x,
// posé aux coins des cellules comme au jour 4 (max |u| \approx A).
//
// Pourquoi ce champ : \omega = -\Delta\psi est proportionnel à \psi, donc (u.\nabla )\omega = 0. C'est une
// solution stationnaire d'Euler. Avec v = 0, la vraie énergie est constante :
// toute perte mesurée est de la dissipation numérique. Avec v > 0 (glissement
// parfait), c'est un mode propre de Stokes : E(t) = E(0) exp(-2 v k^2 t),
// k^2 = \pi ^2/Lx^2 + \pi ^2/Ly^2.
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
