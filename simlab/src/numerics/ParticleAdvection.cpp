#include "numerics/ParticleAdvection.h"
#include <algorithm>
#include <cmath>

namespace {

    // Vitesse interpolée en x à partir des faces u et v.
    // -> Reprends EXACTEMENT ce qu'utilise le advectParticles du jour 7 (par ex. Field2D::sample).
    Vec2 velocityAt(const MacGrid2D& g, Vec2 x) {
        return { g.u.sample(x), g.v.sample(x) };
    }

    // Plus grande composante de vitesse sur toute la grille.
    real maxGridSpeed(const MacGrid2D& g) {
        real m = 0;
        for (real u : g.u.data()) m = std::max(m, std::abs(u));
        for (real v : g.v.data()) m = std::max(m, std::abs(v));
        return m;
    }

} // namespace

AdvectionStats advectParticlesRK2(ParticleSet2D& p, const MacGrid2D& g, real dt,
    real cflSub, const SolidSdf2D& solids, real margin) {
    AdvectionStats st;
    st.vmax = maxGridSpeed(g);

    // Le champ de la grille ne change pas pendant le pas : le nombre de sous-pas se calcule une fois.
    const real maxDisp = st.vmax * dt;                    // distance maximale parcourue en dt
    const int  n = std::max(1, int(std::ceil(maxDisp / (cflSub * g.h()))));
    const real ds = dt / n;
    st.substeps = n;

    for (int s = 0; s < n; ++s) {
        for (std::size_t k = 0; k < p.x.size(); ++k) {
            const Vec2 x0 = p.x[k];
            const Vec2 xm = x0 + real(0.5) * ds * velocityAt(g, x0);   // point milieu
            p.x[k] = x0 + ds * velocityAt(g, xm);                       // pas complet
        }
        pushOutOfSolids(p, solids, margin);               // après CHAQUE sous-pas
    }
    return st;
}