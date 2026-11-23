#pragma once
// ============================================================================
//  Sable : frottement aux parois et contrainte de frottement — jour 13
//  (Zhu & Bridson 2005, §3.1 à 3.3, éq. 1, 2 et 4)
// ============================================================================
//
//  Place dans FlipSolver::step (matériau Sand seulement) :
//
//    5.   gravité + souris
//    5b.  applyWallFriction        (éq. 4)   AVANT la projection : on a besoin de
//                                            la vitesse normale avant qu'elle soit annulée
//    6-7. projection                         (son enforceSolidBC annule u·n)
//    7b.  applyFrictionalStress    (éq. 2)   APRÈS la projection : utilise grid.p
//    8-9. extrapolation n°2 + parois
//
//  Constantes en 2D (déformation plane, sable sans cohésion) :
//    contrainte moyenne   σ_m = tr(σ) / 2                       (et non tr/3)
//    seuil (éq. 1 en 2D)  |dev σ|_F ≤ k p,     k = √2 sin φ     (et non √3 en 3D)
//    éq. 2 en 2D          σ_f = k p D̂,         D̂ = dev D / |dev D|_F
//    vérification         cisaillement simple ⇒ |σ_xy| = sin φ · p
//
//  Convention de signe (CONVENTIONS.md : σ = −p I) : σ_f est DISSIPATIVE, c'est une
//  viscosité effective  η = k p / (2 |D|)  ⇒  Δu = (Δt/ρ) ∇·σ_f fait baisser l'énergie.
//
//  Cellules sous le seuil (éq. 3, §3.2) : σ_rigid = (ρ h² / Δt) dev D, la contrainte qui
//  arrête le glissement de la cellule en un pas. Si |σ_rigid|_F ≤ k p, la cellule est
//  « rigide » et l'on stocke σ_rigid au lieu de σ_f (étape 4 de l'algorithme §3.2).
//  Sans cela, σ_f (de norme k p quelle que soit |D|) dépasse et INVERSE les petites
//  vitesses : le schéma explicite crée de l'énergie. C'est l'équivalent d'une borne
//  η ≤ ρ h² / (2 Δt) sur la viscosité effective.
//
//  Groupes rigides (étape 5 de §3.2, avancée du J14 au J13) : les cellules rigides voisines
//  (4-voisinage) forment des groupes ; les faces de chaque groupe sont remplacées par le
//  mouvement rigide U + ω × (x − c) le plus proche (moindres carrés sur les faces), et ne
//  reçoivent pas la mise à jour (Δt/ρ) ∇·σ. C'est ce qui arrête le sable : sans les groupes,
//  σ_rigid seule laisse le tas fluer comme un fluide très visqueux, et il s'aplatit.
//  La projection moindres carrés conserve la quantité de mouvement et le moment cinétique
//  du groupe (masses égales sur les faces), et ne peut pas créer d'énergie.
// ============================================================================
#include "grid/MacGrid2D.h"
#include "grid/CellLabels.h"

#include <vector>

// ---------------------------------------------------------------------------
//  Frottement aux parois (éq. 4, d'après Bridson et al. 2002)
// ---------------------------------------------------------------------------
//      u_T ← max(0, 1 − μ |u·n| / |u_T|) u_T      si u·n < 0 (poussé vers le mur)
//
//  n : normale sortant du solide, vers le fluide. Hors du domaine = Solid,
//  donc les 4 murs et les obstacles sont traités de la même façon.
struct WallFrictionStats {
    int  frictionFaces = 0;     // faces tangentielles effectivement freinées
    int  stoppedFaces = 0;     // faces ramenées exactement à 0 (μ|u_N| ≥ |u_T|)
};

WallFrictionStats applyWallFriction(MacGrid2D& g, const CellLabels& L, real mu);

// ---------------------------------------------------------------------------
//  Contrainte de frottement (éq. 2), appliquée après la projection
// ---------------------------------------------------------------------------
struct FrictionalStressOptions {
    real frictionAngleDeg = 30;           // φ, en degrés
    real density = 1;            // ρ (le même que dans la projection)
    bool rigidStress = true;         // éq. 3 : σ_rigid dans les cellules sous le seuil
    bool rigidGroups = true;         // étape 5 : groupes rigides projetés (exige rigidStress)
    real minStrainRate = real(1e-8);   // |D| en dessous : σ = 0 (pas de division par 0)
};

struct FrictionalStressStats {
    int  stressCells = 0;     // cellules Fluid où σ_f ≠ 0
    int  rigidCells = 0;     // parmi elles, celles sous le seuil (σ = σ_rigid, éq. 3)
    int  rigidGroups = 0;     // nombre de groupes rigides (étape 5)
    int  rigidFaces = 0;     // faces réécrites par un mouvement rigide
    real maxStress = 0;     // max |σ_f|_F
    real energyBefore = 0;     // ½ Σ |u|² h² sur les faces (contrôle : doit baisser)
    real energyAfter = 0;
};

// Lit grid.p (pression physique, rangée par project()), modifie grid.u et grid.v.
FrictionalStressStats applyFrictionalStress(MacGrid2D& g, const CellLabels& L, real dt,
    const FrictionalStressOptions& opt = {});

// ---------------------------------------------------------------------------
//  Briques exposées pour les tests (test_sand.cpp)
// ---------------------------------------------------------------------------
// Tenseur symétrique 2×2 : [[xx, xy], [xy, yy]].
struct Sym2 {
    real xx = 0, xy = 0, yy = 0;
};

// |S|_F = √(xx² + 2 xy² + yy²)
real frobeniusNorm(const Sym2& S);

// Partie déviatorique en 2D : S − (tr S / 2) I
Sym2 deviator(const Sym2& S);

// Taux de déformation D = (∇u + ∇uᵀ)/2 au centre de la cellule (i, j).
//   ∂u/∂x, ∂v/∂y : directement depuis les faces de la cellule ;
//   ∂u/∂y, ∂v/∂x : calculés aux 4 coins de la cellule, puis moyennés.
//   Coin sur le bord du domaine : ∂u/∂y = ∂v/∂x = 0 (et σ_xy = 0 au même coin, voir le .cpp).
Sym2 strainRateAtCell(const MacGrid2D& g, int i, int j);

// k = √2 sin φ (φ en degrés) : coefficient du seuil en 2D.
real frictionCoefficient2D(real frictionAngleDeg);

// Contrainte d'une cellule (éq. 2 et 3 en 2D) :
//   σ_rigid = rigidScale · dev D           (rigidScale = ρ h² / Δt ; ≤ 0 : désactivé)
//   si |σ_rigid|_F ≤ k p : σ = σ_rigid, *rigid = true          (sous le seuil, éq. 3)
//   sinon               : σ = σ_f = k p dev D / |dev D|_F        (sur le seuil, éq. 2)
//   p : pression de la cellule, déjà bornée à p ≥ 0
//   |dev D|_F < minStrainRate  ⇒  σ = 0
Sym2 frictionalStressAtCell(const Sym2& D, real p, real k, real rigidScale, real minStrainRate,
    bool* rigid = nullptr);

// Étape 5 : numérote les composantes connexes (4-voisinage) d'un masque de cellules.
//   mask[i + nx j] ≠ 0 : cellule rigide.  group[i + nx j] : numéro du groupe, −1 sinon.
//   Renvoie le nombre de groupes.
int labelRigidGroups(const std::vector<char>& mask, int nx, int ny, std::vector<int>& group);

// Mouvement rigide 2D : u(x) = U + ω × (x − c), soit (U.x − ω (y − c_y), U.y + ω (x − c_x)).
struct RigidMotion2D {
    Vec2 U{ 0, 0 };
    real omega = 0;
    Vec2 c{ 0, 0 };
};

// Étape 5 : pour chaque groupe, ajuste le mouvement rigide aux faces FLUIDES du groupe
// (une face appartient au groupe si une de ses deux cellules en fait partie), au sens des
// moindres carrés, puis réécrit ces faces. Renvoie le mouvement de chaque groupe.
std::vector<RigidMotion2D> projectRigidGroups(MacGrid2D& g, const CellLabels& L,
    const std::vector<int>& group, int groupCount,
    int* facesWritten = nullptr);

// Facteur de l'éq. 4 : max(0, 1 − μ |uN| / |uT|) si uN < 0, sinon 1.
// (uN mesuré avec la normale sortant du solide ; protège le cas uT = 0.)
real wallFrictionFactor(real uT, real uN, real mu);