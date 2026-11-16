#include "numerics/FreeSurface.h"
#include <cmath>
#include <cstdint>
#include <vector>

void labelCells(CellLabels& L, const ParticleSet2D& p, real h, const SolidSdf2D* solids,
                int maxBubbleCells) {
    const int nx = L.nx(), ny = L.ny();

    // Passe 1 : une seule boucle sur les particules, chacune coche sa case.
    std::vector<std::uint8_t> hasParticle(std::size_t(nx) * ny, 0);
    for (const Vec2& x : p.x) {
        const int i = int(std::floor(x.x / h));
        const int j = int(std::floor(x.y / h));
        if (i >= 0 && i < nx && j >= 0 && j < ny) hasParticle[i + nx * j] = 1;
    }

    // Passe 2 : boucle sur les cases.
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const Vec2 c{(i + real(0.5)) * h, (j + real(0.5)) * h};   // centre de la case
            if (solids && solids->phi(c) < 0)  L.set(i, j, CellType::Solid);
            else if (hasParticle[i + nx * j])  L.set(i, j, CellType::Fluid);
            else                               L.set(i, j, CellType::Air);
        }

    // Passe 3 : les petites bulles isolées deviennent Fluid.
    // Une cellule vide entourée d'eau aurait p = 0 : sous la pression hydrostatique,
    // le fluide voisin s'y précipiterait en un pas (vitesse ~ g·profondeur·dt/h),
    // entasserait les particules, viderait d'autres cellules, et le volume fondrait.
    // On regroupe les cellules d'air en composantes connexes (4 voisins) ; une
    // composante d'au plus maxBubbleCells cellules est une bulle → Fluid.
    if (maxBubbleCells > 0) {
        std::vector<std::uint8_t> seen(std::size_t(nx) * ny, 0);
        std::vector<int> component, stack;
        const int di[4] = { -1, 1, 0, 0 }, dj[4] = { 0, 0, -1, 1 };
        for (int j0 = 0; j0 < ny; ++j0)
            for (int i0 = 0; i0 < nx; ++i0) {
                const int k0 = i0 + nx * j0;
                if (seen[k0] || !L.isAir(i0, j0)) continue;
                component.clear();
                stack.assign(1, k0);
                seen[k0] = 1;
                while (!stack.empty()) {
                    const int k = stack.back();
                    stack.pop_back();
                    component.push_back(k);
                    const int i = k % nx, j = k / nx;
                    for (int d = 0; d < 4; ++d) {
                        const int ii = i + di[d], jj = j + dj[d];
                        if (!L.isAir(ii, jj)) continue;          // hors domaine = Solid
                        const int kk = ii + nx * jj;
                        if (!seen[kk]) { seen[kk] = 1; stack.push_back(kk); }
                    }
                }
                if (int(component.size()) <= maxBubbleCells)
                    for (int k : component) L.set(k % nx, k / nx, CellType::Fluid);
            }
    }
}

void enforceSolidBC(MacGrid2D& g, const CellLabels& L) {
    const int nx = g.nx(), ny = g.ny();
    // Faces u : (nx+1) × ny. Sur une grille MAC, la normale d'une face u est ±x.
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i <= nx; ++i)
            if (L.isSolid(i - 1, j) || L.isSolid(i, j)) g.u(i, j) = 0;
    // Faces v : nx × (ny+1). Normale ±y.
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i < nx; ++i)
            if (L.isSolid(i, j - 1) || L.isSolid(i, j)) g.v(i, j) = 0;
}

bool uFaceIsFluid(const CellLabels& L, int i, int j) {
    const CellType a = L(i - 1, j), b = L(i, j);
    if (a == CellType::Solid || b == CellType::Solid) return false;
    return a == CellType::Fluid || b == CellType::Fluid;
}

bool vFaceIsFluid(const CellLabels& L, int i, int j) {
    const CellType a = L(i, j - 1), b = L(i, j);
    if (a == CellType::Solid || b == CellType::Solid) return false;
    return a == CellType::Fluid || b == CellType::Fluid;
}
