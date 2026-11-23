#include "diagnostics/SandDiagnostics.h"

#include <algorithm>
#include <cmath>

std::vector<real> heightProfile(const ParticleSet2D& p, int nx, real h) {
    std::vector<real> top(std::size_t(nx), 0);
    for (const Vec2& x : p.x) {
        const int i = std::clamp(int(std::floor(x.x / h)), 0, nx - 1);
        top[i] = std::max(top[i], x.y);
    }
    return top;
}

namespace {

    // Pente du flanc qui part du sommet dans la direction step (−1 : gauche, +1 : droite).
    // Renvoie false si le flanc a moins de 3 colonnes utilisables.
    bool flankSlope(const std::vector<real>& top, int peak, int step, real h, real H, real lo, real hi,
                    real& slope) {
        real sx = 0, sy = 0, sxx = 0, sxy = 0;
        int m = 0;
        for (int i = peak; i >= 0 && i < int(top.size()); i += step) {
            if (top[i] < lo * H) break;        // pied du flanc
            if (top[i] > hi * H) continue;     // sommet arrondi : on l'ignore
            const real x = (i + real(0.5)) * h, y = top[i];
            sx += x; sy += y; sxx += x * x; sxy += x * y; ++m;
        }
        if (m < 3) return false;
        const real den = m * sxx - sx * sx;
        if (!(den > 0)) return false;
        slope = (m * sxy - sx * sy) / den;
        return true;
    }

} // namespace

PileShape measurePile(const ParticleSet2D& p, int nx, real h, real lo, real hi) {
    PileShape r;
    if (p.x.empty()) return r;
    r.xMin = p.x[0].x;
    r.xMax = p.x[0].x;
    for (const Vec2& x : p.x) {
        r.xMin = std::min(r.xMin, x.x);
        r.xMax = std::max(r.xMax, x.x);
    }
    r.spread = r.xMax - r.xMin;

    const std::vector<real> top = heightProfile(p, nx, h);
    const int peak = int(std::max_element(top.begin(), top.end()) - top.begin());
    r.height = top[peak];
    r.peakX = (peak + real(0.5)) * h;

    real slope = 0, sum = 0;
    if (flankSlope(top, peak, -1, h, r.height, lo, hi, slope)) {
        r.angleLeftDeg = std::atan(std::abs(slope)) * 180 / kPi;
        sum += r.angleLeftDeg; ++r.flanks;
    }
    if (flankSlope(top, peak, +1, h, r.height, lo, hi, slope)) {
        r.angleRightDeg = std::atan(std::abs(slope)) * 180 / kPi;
        sum += r.angleRightDeg; ++r.flanks;
    }
    r.angleDeg = r.flanks ? sum / r.flanks : 0;

    // Base continue autour du sommet, au-dessus de 0,1·H.
    int a = peak, b = peak;
    while (a > 0 && top[a - 1] >= real(0.1) * r.height) --a;
    while (b + 1 < int(top.size()) && top[b + 1] >= real(0.1) * r.height) ++b;
    r.baseWidth = (b - a + 1) * h;
    r.angleEquivalentDeg = r.height > 0 ? std::atan(2 * r.height / r.baseWidth) * 180 / kPi : 0;
    return r;
}
