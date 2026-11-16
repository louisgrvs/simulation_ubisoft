// Extrapole f sur au plus `layers` couches à partir des faces valides.
// Met à jour `valid` (les faces remplies deviennent valides).
// Renvoie le nombre de couches réellement faites (moins si plus rien à remplir).
// But : prolonger la vitesse du fluide vers l'air
#include "numerics/Extrapolation.h"
#include <algorithm>
#include <cmath>
#include "numerics/FreeSurface.h"   // uFaceIsFluid, vFaceIsFluid

int extrapolateLayers(Field2D& f, FaceMask& valid, int layers) {
    const int nx = f.nx(), ny = f.ny();   // dimensions du champ de FACES (nx+1 pour u, etc.)

    // Les 4 voisins : gauche, droite, bas, haut.
    const int di[4] = { -1, 1, 0, 0 };
    const int dj[4] = { 0, 0, -1, 1 };

    for (int layer = 1; layer <= layers; ++layer) {        // pour couche = 1 à layers
        FaceMask next = valid;                              // nouveau = copie de valid
        int filled = 0;                                     // rempli = 0

        for (int j = 0; j < ny; ++j)                        // pour chaque face (i, j)
            for (int i = 0; i < nx; ++i) {
                const int k = i + nx * j;
                if (valid[k]) continue;                     // déjà valide : on passe

                real sum = 0;                               // somme = 0 ; n = 0
                int n = 0;
                for (int d = 0; d < 4; ++d) {               // pour chacun des 4 voisins
                    const int ii = i + di[d], jj = j + dj[d];
                    if (ii < 0 || ii >= nx || jj < 0 || jj >= ny) continue;   // hors du champ
                    if (valid[ii + nx * jj]) {              // voisin valide (ancien masque)
                        sum += f(ii, jj);
                        ++n;
                    }
                }

                if (n > 0) {                                // si n > 0
                    f(i, j) = sum / n;                      //   moyenne des voisins valides
                    next[k] = 1;                            //   valide à la fin de la couche
                    ++filled;
                }
            }

        valid = next;                                       // valid = nouveau
        if (filled == 0) return layer - 1;                  // plus rien à remplir
    }
    return layers;
}

void extrapolateVelocity(MacGrid2D& g, FaceMask& uValid, FaceMask& vValid, int layers) {
    extrapolateLayers(g.u, uValid, layers);
    extrapolateLayers(g.v, vValid, layers);
}

FaceMask maskFromWeights(const Field2D& weights, real eps) {
    const auto& w = weights.data();
    FaceMask m(w.size(), 0);
    for (std::size_t k = 0; k < w.size(); ++k) m[k] = (w[k] > eps) ? 1 : 0;
    return m;
}

FaceMask uMaskFromLabels(const CellLabels& L) {
    const int nx = L.nx(), ny = L.ny();
    FaceMask m(std::size_t(nx + 1) * ny, 0);               // faces u : (nx+1) × ny
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i <= nx; ++i)
            if (uFaceIsFluid(L, i, j)) m[i + (nx + 1) * j] = 1;
    return m;
}

FaceMask vMaskFromLabels(const CellLabels& L) {
    const int nx = L.nx(), ny = L.ny();
    FaceMask m(std::size_t(nx) * (ny + 1), 0);             // faces v : nx × (ny+1)
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i < nx; ++i)
            if (vFaceIsFluid(L, i, j)) m[i + nx * j] = 1;
    return m;
}

int extrapolationLayerCount(real vmax, real dt, real h) {
    return std::max(2, int(std::ceil(vmax * dt / h)) + 2);
}