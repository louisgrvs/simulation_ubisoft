#include "numerics/SandFriction.h"
// ============================================================================
//  Sable : frottement aux parois (éq. 4) et contrainte de frottement (éq. 2, 3) — J13
// ============================================================================
#include <algorithm>
#include <cmath>
#include <vector>
#include "grid/Field2D.h"
#include "numerics/FreeSurface.h"   // uFaceIsFluid, vFaceIsFluid

namespace {

    // Énergie cinétique de la grille 1/2 \sum |u|² h² (faces u puis faces v).
    // Sert au contrôle « la friction ne crée pas d'énergie ».
    // Les deux grilles n'ont PAS la même taille : u est (nx+1) × ny, v est nx × (ny+1),
    // d'où deux boucles séparées.
    real faceEnergy(const MacGrid2D& g) {
        real res = 0;
        for (int j = 0; j < g.ny(); ++j)
            for (int i = 0; i <= g.nx(); ++i) res += g.u(i, j) * g.u(i, j);
        for (int j = 0; j <= g.ny(); ++j)
            for (int i = 0; i < g.nx(); ++i) res += g.v(i, j) * g.v(i, j);
        return real(0.5) * res * g.h() * g.h();   // real(0.5) : 1/2 entre entiers vaut 0
    }

    // Composante \sigma_xy au coin (i, j) (coin bas-gauche de la cellule (i, j)) :
    //   SOMME des 4 cellules voisines / 4, une cellule sans contrainte comptant 0
    //   (pas la moyenne des seules cellules présentes !). Coin sur le bord du domaine : 0.
    // Pourquoi : c'est l'opération ADJOINTE de « moyenne des 4 coins → centre » faite
    // dans strainRateAtCell. Avec elle, \sum u·\Delta u = −(dt/\rho) \sum \sigma:D h² ≤ 0
    // exactement, et le test d'énergie passe.
    // Coins : (nx+1) × (ny+1) ; les 4 cellules autour du coin (i, j) sont
    // (i−1, j−1), (i, j−1), (i−1, j), (i, j).
    real stressXYAtNode(const Field2D& sxy, int i, int j) {
        if (i <= 0 || j <= 0 || i >= sxy.nx() || j >= sxy.ny()) return 0;   // coin sur le bord
        return real(0.25) * (sxy(i - 1, j - 1) + sxy(i, j - 1) + sxy(i - 1, j) + sxy(i, j));
    }

    // \partial u/\partial y + \partial v/\partial x au coin (i, j) ; 0 sur le bord du domaine
    // (même convention que stressXYAtNode, pour que les deux opérations restent adjointes).
    real shearAtNode(const MacGrid2D& g, int i, int j) {
        if (i <= 0 || j <= 0 || i >= g.nx() || j >= g.ny()) return 0;
        const real h = g.h();
        const real dudy = (g.u(i, j) - g.u(i, j - 1)) / h;   // u(i, j) et u(i, j−1) encadrent le coin
        const real dvdx = (g.v(i, j) - g.v(i - 1, j)) / h;   // v(i, j) et v(i−1, j) aussi
        return dudy + dvdx;
    }

} // namespace

// ---------------------------------------------------------------------------
//  Briques
// ---------------------------------------------------------------------------
real frobeniusNorm(const Sym2& S) {
    return std::sqrt(S.xx * S.xx + 2 * S.xy * S.xy + S.yy * S.yy);   // la racine !
}

Sym2 deviator(const Sym2& S) {
    // Déviateur = on retire la partie sphérique (tr S / 2) I : ce n'est pas la
    // symétrisation (Sym2 est déjà symétrique). En 2D, tr/2 et non tr/3.
    const real m = real(0.5) * (S.xx + S.yy);
    return { S.xx - m, S.xy, S.yy - m };
}

Sym2 strainRateAtCell(const MacGrid2D& g, int i, int j) {
    const real h = g.h();
    Sym2 D;
    D.xx = (g.u(i + 1, j) - g.u(i, j)) / h;
    D.yy = (g.v(i, j + 1) - g.v(i, j)) / h;
    // Les 4 coins de la cellule (i, j) : (i, j), (i+1, j), (i, j+1), (i+1, j+1).
    real s = 0;
    for (int b = j; b <= j + 1; ++b)
        for (int a = i; a <= i + 1; ++a) s += shearAtNode(g, a, b);
    D.xy = real(0.5) * (s / 4);   // 1/2 (moyenne de \partial u/\partial y + \partial v/\partial x)
    return D;
}

real frictionCoefficient2D(real frictionAngleDeg) {
    return std::sqrt(real(2)) * std::sin(frictionAngleDeg * kPi / 180);
}

Sym2 frictionalStressAtCell(const Sym2& D, real p, real k, real rigidScale, real minStrainRate,
    bool* rigid) {
    if (rigid) *rigid = false;
    const Sym2 Dd = deviator(D);
    const real n = frobeniusNorm(Dd);
    if (n < minStrainRate) return {};                       // pas de division par 0
    // Éq. 3 : sous le seuil, la contrainte qui arrête la cellule en un pas.
    if (rigidScale > 0 && rigidScale * n <= k * p) {
        if (rigid) *rigid = true;
        return { rigidScale * Dd.xx, rigidScale * Dd.xy, rigidScale * Dd.yy };
    }
    // Éq. 2 : sur le seuil, \sigma_f = k p D̂ (= 2 \eta dev D avec \eta = k p / (2 n)).
    const real c = k * p / n;
    return { c * Dd.xx, c * Dd.xy, c * Dd.yy };               // Sym2 n'a pas d'opérateur *
}

real wallFrictionFactor(real uT, real uN, real mu) {
    if (uN >= 0 || uT == 0) return 1;                       // s'éloigne du mur, ou rien à freiner
    return std::max(real(0), 1 - mu * std::abs(uN) / std::abs(uT));
}

// ---------------------------------------------------------------------------
//  Étape 5 : groupes rigides
// ---------------------------------------------------------------------------
int labelRigidGroups(const std::vector<char>& mask, int nx, int ny, std::vector<int>& group) {
    group.assign(std::size_t(nx) * ny, -1);
    int count = 0;
    std::vector<int> stack;
    for (int start = 0; start < nx * ny; ++start) {
        if (!mask[start] || group[start] >= 0) continue;
        // Parcours en profondeur avec une pile (comme les poches d'air de labelCells, J8).
        group[start] = count;
        stack.assign(1, start);
        while (!stack.empty()) {
            const int c = stack.back();
            stack.pop_back();
            const int i = c % nx, j = c / nx;
            const int ni[4] = { i + 1, i - 1, i, i }, nj[4] = { j, j, j + 1, j - 1 };
            for (int k = 0; k < 4; ++k) {
                if (ni[k] < 0 || nj[k] < 0 || ni[k] >= nx || nj[k] >= ny) continue;
                const int q = ni[k] + nx * nj[k];
                if (mask[q] && group[q] < 0) { group[q] = count; stack.push_back(q); }
            }
        }
        ++count;
    }
    return count;
}

namespace {

    // Groupe d'une cellule (−1 hors du domaine ou hors groupe).
    int groupOf(const std::vector<int>& group, int nx, int ny, int i, int j) {
        return (i < 0 || j < 0 || i >= nx || j >= ny) ? -1 : group[i + nx * j];
    }
    // Groupe d'une face : celui d'une de ses deux cellules (deux cellules rigides voisines
    // sont toujours dans le même groupe, donc il n'y a jamais de conflit).
    int uFaceGroup(const std::vector<int>& group, int nx, int ny, int i, int j) {
        const int a = groupOf(group, nx, ny, i - 1, j);
        return a >= 0 ? a : groupOf(group, nx, ny, i, j);
    }
    int vFaceGroup(const std::vector<int>& group, int nx, int ny, int i, int j) {
        const int a = groupOf(group, nx, ny, i, j - 1);
        return a >= 0 ? a : groupOf(group, nx, ny, i, j);
    }

} // namespace

std::vector<RigidMotion2D> projectRigidGroups(MacGrid2D& g, const CellLabels& L,
    const std::vector<int>& group, int groupCount,
    int* facesWritten) {
    const int nx = g.nx(), ny = g.ny();
    const real h = g.h();
    std::vector<RigidMotion2D> motion(groupCount);
    if (facesWritten) *facesWritten = 0;
    if (groupCount == 0) return motion;

    // Centre c de chaque groupe : centre des cellules.
    std::vector<int> cells(groupCount, 0);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int q = group[i + nx * j];
            if (q < 0) continue;
            motion[q].c += Vec2{ (i + real(0.5)) * h, (j + real(0.5)) * h };
            ++cells[q];
        }
    for (int q = 0; q < groupCount; ++q) motion[q].c = (real(1) / cells[q]) * motion[q].c;

    // Sommes des moindres carrés. Avec x̂ = x − c_x, ŷ = y − c_y :
    //   faces u : u ≈ U.x − ω ŷ      faces v : v ≈ U.y + ω x̂
    struct Sums { real nu = 0, su = 0, sy = 0, syy = 0, syu = 0, nv = 0, sv = 0, sx = 0, sxx = 0, sxv = 0; };
    std::vector<Sums> S(groupCount);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            if (!uFaceIsFluid(L, i, j)) continue;
            const int q = uFaceGroup(group, nx, ny, i, j);
            if (q < 0) continue;
            const real y = (j + real(0.5)) * h - motion[q].c.y, u = g.u(i, j);
            S[q].nu += 1; S[q].su += u; S[q].sy += y; S[q].syy += y * y; S[q].syu += y * u;
        }
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i < nx; ++i) {
            if (!vFaceIsFluid(L, i, j)) continue;
            const int q = vFaceGroup(group, nx, ny, i, j);
            if (q < 0) continue;
            const real x = (i + real(0.5)) * h - motion[q].c.x, v = g.v(i, j);
            S[q].nv += 1; S[q].sv += v; S[q].sx += x; S[q].sxx += x * x; S[q].sxv += x * v;
        }

    // Équations normales (3 × 3), résolues par élimination de U :
    //   U.x = (Σu + ω Σŷ) / n_u ,   U.y = (Σv − ω Σx̂) / n_v ,
    //   ω [Σŷ² + Σx̂² − (Σŷ)²/n_u − (Σx̂)²/n_v] = Σx̂ v − Σŷ u + Σŷ Σu / n_u − Σx̂ Σv / n_v.
    // Le crochet est une inertie (≥ 0). S'il est nul (groupe d'une cellule, par exemple),
    // la rotation n'est pas mesurable : ω = 0.
    for (int q = 0; q < groupCount; ++q) {
        const Sums& t = S[q];
        real inertia = t.syy + t.sxx, rhs = t.sxv - t.syu;
        if (t.nu > 0) { inertia -= t.sy * t.sy / t.nu; rhs += t.sy * t.su / t.nu; }
        if (t.nv > 0) { inertia -= t.sx * t.sx / t.nv; rhs -= t.sx * t.sv / t.nv; }
        RigidMotion2D& m = motion[q];
        m.omega = inertia > real(1e-12) * (t.syy + t.sxx + h * h) ? rhs / inertia : 0;
        m.U.x = t.nu > 0 ? (t.su + m.omega * t.sy) / t.nu : 0;
        m.U.y = t.nv > 0 ? (t.sv - m.omega * t.sx) / t.nv : 0;
    }

    // Réécriture des faces du groupe.
    int written = 0;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            if (!uFaceIsFluid(L, i, j)) continue;
            const int q = uFaceGroup(group, nx, ny, i, j);
            if (q < 0) continue;
            const RigidMotion2D& m = motion[q];
            g.u(i, j) = m.U.x - m.omega * ((j + real(0.5)) * h - m.c.y);
            ++written;
        }
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i < nx; ++i) {
            if (!vFaceIsFluid(L, i, j)) continue;
            const int q = vFaceGroup(group, nx, ny, i, j);
            if (q < 0) continue;
            const RigidMotion2D& m = motion[q];
            g.v(i, j) = m.U.y + m.omega * ((i + real(0.5)) * h - m.c.x);
            ++written;
        }
    if (facesWritten) *facesWritten = written;
    return motion;
}

// ---------------------------------------------------------------------------
//  Éq. 4 : frottement aux parois
// ---------------------------------------------------------------------------
WallFrictionStats applyWallFriction(MacGrid2D& g, const CellLabels& L, real mu) {
    WallFrictionStats st;
    const int nx = g.nx(), ny = g.ny();
    // On lit les vitesses d'AVANT (copies) : le résultat ne dépend pas de l'ordre de parcours.
    const Field2D u0 = g.u, v0 = g.v;
    // Une face touche un solide si une de ses deux cellules est Solid (hors domaine = Solid).
    auto vTouchesSolid = [&](int i, int j) { return L.isSolid(i, j - 1) || L.isSolid(i, j); };
    auto uTouchesSolid = [&](int i, int j) { return L.isSolid(i - 1, j) || L.isSolid(i, j); };

    // Applique le facteur de l'éq. 4 à t, calculé avec la tangentielle d'avant tOld.
    auto brake = [&](real& t, real tOld, real uN) {
        const real f = wallFrictionFactor(tOld, uN, mu);
        if (f < 1) {
            t *= f;
            ++st.frictionFaces;
            if (f == 0) ++st.stoppedFaces;
        }
        };
    // Moyenne des normales (copie) des deux faces voisines qui touchent un solide.
    // Renvoie false si aucune ne touche de solide.
    auto meanNormalV = [&](int i0, int i1, int j, real& uN) {
        real s = 0; int c = 0;
        for (int a : { i0, i1 })
            if (a >= 0 && a < nx && vTouchesSolid(a, j)) { s += v0(a, j); ++c; }
        if (c) uN = s / c;
        return c > 0;
        };
    auto meanNormalU = [&](int i, int j0, int j1, real& uN) {
        real s = 0; int c = 0;
        for (int b : { j0, j1 })
            if (b >= 0 && b < ny && uTouchesSolid(i, b)) { s += u0(i, b); ++c; }
        if (c) uN = s / c;
        return c > 0;
        };

    // Faces u : tangentielles aux sols et plafonds.
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            if (!uFaceIsFluid(L, i, j)) continue;
            real uN = 0;
            if ((L.isSolid(i - 1, j - 1) || L.isSolid(i, j - 1)) && meanNormalV(i - 1, i, j, uN))
                brake(g.u(i, j), u0(i, j), uN);           // sol : n = +y
            if ((L.isSolid(i - 1, j + 1) || L.isSolid(i, j + 1)) && meanNormalV(i - 1, i, j + 1, uN))
                brake(g.u(i, j), u0(i, j), -uN);          // plafond : n = −y
        }
    // Faces v : tangentielles aux murs gauche et droit.
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i < nx; ++i) {
            if (!vFaceIsFluid(L, i, j)) continue;
            real uN = 0;
            if ((L.isSolid(i - 1, j - 1) || L.isSolid(i - 1, j)) && meanNormalU(i, j - 1, j, uN))
                brake(g.v(i, j), v0(i, j), uN);           // mur de gauche : n = +x
            if ((L.isSolid(i + 1, j - 1) || L.isSolid(i + 1, j)) && meanNormalU(i + 1, j - 1, j, uN))
                brake(g.v(i, j), v0(i, j), -uN);          // mur de droite : n = −x
        }
    return st;
}

// ---------------------------------------------------------------------------
//  Éq. 2 et 3 : contrainte de frottement
// ---------------------------------------------------------------------------
FrictionalStressStats applyFrictionalStress(MacGrid2D& g, const CellLabels& L, real dt,
    const FrictionalStressOptions& opt) {
    FrictionalStressStats st;
    st.energyBefore = faceEnergy(g);
    const int nx = g.nx(), ny = g.ny();
    const real h = g.h(), rho = opt.density;
    const real k = frictionCoefficient2D(opt.frictionAngleDeg);
    const real rigidScale = opt.rigidStress ? rho * h * h / dt : 0;

    // 1. Contrainte dans chaque cellule Fluid (0 dans l'air et les solides), et masque
    //    des cellules rigides (éq. 3).
    Field2D sxx = g.makeCellField(), syy = g.makeCellField(), sxy = g.makeCellField();
    std::vector<char> rigidMask(std::size_t(nx) * ny, 0);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            if (!L.isFluid(i, j)) continue;
            const real p = std::max(g.p(i, j), real(0));    // pas de traction dans le sable
            bool rigid = false;
            const Sym2 s = frictionalStressAtCell(strainRateAtCell(g, i, j), p, k, rigidScale,
                opt.minStrainRate, &rigid);
            sxx(i, j) = s.xx; syy(i, j) = s.yy; sxy(i, j) = s.xy;
            const real n = frobeniusNorm(s);
            if (n > 0) ++st.stressCells;
            if (rigid) { ++st.rigidCells; rigidMask[i + nx * j] = 1; }
            st.maxStress = std::max(st.maxStress, n);
        }

    // 2. Étape 5 : groupes de cellules rigides (si demandé).
    std::vector<int> group;
    const bool groups = opt.rigidGroups && opt.rigidStress;
    st.rigidGroups = groups ? labelRigidGroups(rigidMask, nx, ny, group) : 0;
    if (!groups) group.assign(std::size_t(nx) * ny, -1);

    // 3. \Delta u = (dt/\rho) \nabla·\sigma sur les faces fluides HORS des groupes rigides.
    //    Toutes les contraintes sont calculées avant la moindre écriture (étape 1), on peut
    //    donc écrire directement ; les faces des groupes sont réécrites à l'étape 4.
    const real c = dt / (rho * h);
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i <= nx; ++i)
            if (uFaceIsFluid(L, i, j) && uFaceGroup(group, nx, ny, i, j) < 0)   // exclut aussi les murs
                g.u(i, j) += c * ((sxx(i, j) - sxx(i - 1, j))
                    + stressXYAtNode(sxy, i, j + 1) - stressXYAtNode(sxy, i, j));
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i < nx; ++i)
            if (vFaceIsFluid(L, i, j) && vFaceGroup(group, nx, ny, i, j) < 0)   // exclut aussi les murs
                g.v(i, j) += c * ((stressXYAtNode(sxy, i + 1, j) - stressXYAtNode(sxy, i, j))
                    + syy(i, j) - syy(i, j - 1));

    // 4. Étape 5 : chaque groupe rigide prend le mouvement rigide le plus proche.
    if (st.rigidGroups > 0) projectRigidGroups(g, L, group, st.rigidGroups, &st.rigidFaces);

    st.energyAfter = faceEnergy(g);
    return st;
}