#include "diagnostics/FreeSurfaceDiagnostics.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "numerics/MacOperators.h"   // maxAbsDivergence(g, L)

bool hydrostaticReference(const CellLabels& L, real h, int i, int j, real density, real g, real& pRef) {
    // On remonte la colonne jusqu'à la première cellule qui n'est pas Fluid.
    int jj = j;
    while (jj < L.ny() && L.isFluid(i, jj)) ++jj;
    if (jj >= L.ny() || !L.isAir(i, jj)) return false;   // plafond ou solide : pas de p = 0 au-dessus
    pRef = density * g * real(jj - j) * h;                 // ρ g (y_s − y_j), y en centres de cellules
    return true;
}

FreeSurfaceStats measureFreeSurface(const MacGrid2D& grid, const CellLabels& L, const ParticleSet2D& p,
                                    real density, real g) {
    FreeSurfaceStats s;
    const int nx = L.nx(), ny = L.ny();
    const real h = grid.h();

    // 1. Étiquettes et volume.
    s.fluidCells = L.count(CellType::Fluid);
    s.airCells   = L.count(CellType::Air);
    s.solidCells = L.count(CellType::Solid);
    s.fluidArea  = s.fluidCells * h * h;

    // Hauteur de la colonne d'eau : nombre de cellules fluides contiguës depuis le fond.
    real sumHeight = 0;
    int columns = 0;
    for (int i = 0; i < nx; ++i) {
        int j = 0;
        while (j < ny && L.isFluid(i, j)) ++j;
        if (j > 0) { sumHeight += j * h; ++columns; }
    }
    s.surfaceHeight = columns > 0 ? sumHeight / columns : 0;

    // 2. Vitesses des particules.
    real sum2 = 0;
    for (const Vec2& v : p.v) {
        const real v2 = dot(v, v);
        s.maxSpeed = std::max(s.maxSpeed, std::sqrt(v2));
        sum2 += v2;
    }
    s.rmsSpeed = p.size() > 0 ? std::sqrt(sum2 / real(p.size())) : 0;

    // 3. Écart à la pression hydrostatique, rapporté à la plus grande pression de référence.
    real maxErr = 0, maxRef = 0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            if (!L.isFluid(i, j)) continue;
            real pRef;
            if (!hydrostaticReference(L, h, i, j, density, g, pRef)) continue;
            maxErr = std::max(maxErr, std::abs(grid.p(i, j) - pRef));
            maxRef = std::max(maxRef, pRef);
        }
    s.hydroErr = maxRef > 0 ? maxErr / maxRef : 0;

    s.divAfter = maxAbsDivergence(grid, L);
    return s;
}

void writeFreeSurfaceCsvHeader(std::ostream& os) {
    os << std::setprecision(12);   // 6 chiffres par défaut : trop peu pour des écarts de 1e-12
    os << "t,fluid_cells,fluid_area,surface_height,max_speed,rms_speed,hydro_err,div_after\n";
}

void writeFreeSurfaceCsvRow(std::ostream& os, real t, const FreeSurfaceStats& s) {
    os << t << ',' << s.fluidCells << ',' << s.fluidArea << ',' << s.surfaceHeight << ','
       << s.maxSpeed << ',' << s.rmsSpeed << ',' << s.hydroErr << ',' << s.divAfter << '\n';
}

void writePressureProfileCsv(std::ostream& os, const MacGrid2D& grid, const CellLabels& L, int i,
                             real density, real g) {
    os << std::setprecision(15) << "y,p,p_ref\n";
    const real h = grid.h();
    for (int j = 0; j < L.ny(); ++j) {
        if (!L.isFluid(i, j)) continue;
        real pRef = 0;
        hydrostaticReference(L, h, i, j, density, g, pRef);
        os << (j + real(0.5)) * h << ',' << grid.p(i, j) << ',' << pRef << '\n';
    }
}
