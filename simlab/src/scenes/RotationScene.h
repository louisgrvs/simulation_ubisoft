#pragma once
#include <cmath>
#include "grid/MacGrid2D.h"

// Scène test : rotation solide imposée autour du centre du domaine,
// et disque fendu de Zalesak (1979) comme scalaire transporté.
// Le disque est un champ 0/1 à bords nets : idéal pour voir la diffusion
// numérique (bords qui s'adoucissent) et la qualité des trajectoires (fente).
namespace rotation_scene {

// u(x) = omega * (-(y - cy), x - cx) : divergence nulle, un tour en 2 pi / omega.
// C'est un champ AFFINE : l'interpolation bilinéaire le reproduit exactement.
inline void setRotation(MacGrid2D& g, real omega) {
    const Vec2 c = 0.5 * g.size();
    g.setVelocity([=](Vec2 x) { return Vec2{-omega * (x.y - c.y), omega * (x.x - c.x)}; });
}

// Disque de rayon 0.15 L centré en (0.5 L, 0.75 L), fendu d'une encoche
// de largeur 0.05 L montant depuis le bas du disque jusqu'à y = 0.85 L.
// angle != 0 : solution EXACTE après une rotation d'angle `angle` autour du centre
// (on évalue le disque initial au point d'où vient x).
inline void setSlottedDisk(Field2D& q, real L, real angle = 0) {
    const Vec2 c{0.5 * L, 0.5 * L};
    for (int j = 0; j < q.ny(); ++j)
        for (int i = 0; i < q.nx(); ++i) {
            Vec2 x = (c + rotate(q.position(i, j) - c, -angle)) / L;
            Vec2 d = x - Vec2{0.5, 0.75};
            bool inDisk = norm(d) <= 0.15;
            bool inSlot = std::abs(d.x) <= 0.025 && x.y <= 0.85;
            q(i, j) = (inDisk && !inSlot) ? 1.0 : 0.0;
        }
}

// Masse (intégrale de q) : le semi-lagrangien NE la conserve PAS exactement.
inline real mass(const Field2D& q) { return q.sum() * q.h() * q.h(); }

inline Vec2 centroid(const Field2D& q) {
    Vec2 c{0, 0};
    real m = 0;
    for (int j = 0; j < q.ny(); ++j)
        for (int i = 0; i < q.nx(); ++i) {
            c += q(i, j) * q.position(i, j);
            m += q(i, j);
        }
    return c / m;
}

// Erreur L1 entre deux champs de même disposition : sum |a - b| h^2.
inline real l1Error(const Field2D& a, const Field2D& b) {
    real e = 0;
    for (std::size_t k = 0; k < a.data().size(); ++k) e += std::abs(a.data()[k] - b.data()[k]);
    return e * a.h() * a.h();
}

} // namespace rotation_scene
