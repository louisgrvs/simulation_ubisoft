#include "geometry/SolidSdf2D.h"
#include <algorithm>
#include <cmath>

namespace {
    // Distance signée à un rectangle : < 0 dedans, > 0 dehors.
    real boxDistance(Vec2 x, Vec2 min, Vec2 max) {
        const Vec2 c = 0.5 * (min + max), half = 0.5 * (max - min);
        const real qx = std::abs(x.x - c.x) - half.x;
        const real qy = std::abs(x.y - c.y) - half.y;
        const real outside = std::hypot(std::max(qx, real(0)), std::max(qy, real(0)));
        const real inside = std::min(std::max(qx, qy), real(0));
        return outside + inside;
    }
} // namespace

SolidSdf2D::SolidSdf2D(Vec2 boxSize) : size_(boxSize) {}

void SolidSdf2D::addBox(Vec2 min, Vec2 max) { obstacles_.push_back({ min, max }); }

real SolidSdf2D::phi(Vec2 x) const {
    // Murs de la cuve : φ > 0 à l'intérieur.
    real d = std::min({ x.x, size_.x - x.x, x.y, size_.y - x.y });
    // Obstacles : le fluide est à l'extérieur, donc φ = distance signée au rectangle.
    for (const auto& b : obstacles_) d = std::min(d, boxDistance(x, b.min, b.max));
    return d;
}

Vec2 SolidSdf2D::normal(Vec2 x) const {
    // Pas des différences centrées : bien plus petit que la marge de repoussement
    // (≈ 10⁻³ h), sinon près d'un coin les deux points tombent de part et d'autre d'une arête.
    const real e = real(1e-7) * std::min(size_.x, size_.y);
    const Vec2 n{ phi({ x.x + e, x.y }) - phi({ x.x - e, x.y }),
                 phi({ x.x, x.y + e }) - phi({ x.x, x.y - e }) };
    const real len = norm(n);
    return len > 0 ? n / len : Vec2{ 0, 0 };
}

void pushOutOfSolids(ParticleSet2D& p, const SolidSdf2D& s, real margin) {
    for (std::size_t k = 0; k < p.size(); ++k) {
        // Quelques itérations : dans un coin rentrant (sol + paroi de l'obstacle), sortir
        // d'un côté peut faire entrer de l'autre. En pratique, 1 ou 2 suffisent.
        for (int it = 0; it < 4; ++it) {
            const real d = s.phi(p.x[k]);
            if (d >= margin) break;
            const Vec2 n = s.normal(p.x[k]);
            p.x[k] += (margin - d) * n;              // retour à la distance « marge »
            const real vn = dot(p.v[k], n);
            if (vn < 0) p.v[k] = p.v[k] - vn * n;    // on retire la vitesse qui rentre dans le mur
        }
    }
}