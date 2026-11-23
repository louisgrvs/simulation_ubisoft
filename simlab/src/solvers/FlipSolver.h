#pragma once
// ============================================================================
//  Solveur hybride particules / grille MAC (Zhu & Bridson #4.2)
//  jour 8 : surface libre ; jour 9 : FLIP et rupture de barrage ; jour 10 : APIC
//  jour 11 : disque en rotation libre (pas complet), diagnostics de forme
//  jour 13 : sable (Zhu & Bridson §3) : frottement aux parois, contrainte de frottement
// ============================================================================
//
//  Un pas :
//
//    1. étiquettes Fluid / Air / Solid à partir des particules        (#4.2.2)
//    2. P2G (en APIC : avec la partie affine c_p · (x_f - x_p), J10)   (#4.2.2)
//    3. extrapolation n°1 (faces sans poids) + parois
//    4. uOld_, vOld_ = vitesse avant les forces (pour FLIP)
//    5. gravité + souris                                              (#4.2.3)
//    5b. [Sable, J13] frottement aux parois (éq. 4), AVANT que la projection annule u·n
//    6-7. projection avec p = 0 dans l'air                            (#4.2.3)
//    7b. [Sable, J13] contrainte de frottement (éq. 2, 3) et groupes rigides (étape 5)
//    7c. [Sable, J13] deuxième projection : le frottement n'est pas à divergence nulle
//    8-9. extrapolation n°2 (depuis les faces fluides) + parois
//   10. G2P : PIC, mélange alpha FLIP + (1 - alpha) PIC (J9),          (#4.2.4)
//       ou APIC : vitesse comme PIC + gradient c_p (J10)            (Jiang et al. #6)
//   11. déplacement RK2 en sous-pas CFL + repoussement hors des parois (#4.2.5)
//
//  La scène RotatingDisk ne fait QUE les transferts (étapes 2 et 10).
//  RotatingDiskFree (J11) part du même état mais fait le pas complet, sans gravité :
//  le disque tourne librement, avec une surface libre tout autour.
//  La gravité n'agit que dans TankAtRest, DamBreak et ColumnCollapse (J13) ; les scènes du J7 restent inchangées.
// ============================================================================
#include <string>
#include <vector>
#include "grid/MacGrid2D.h"
#include "grid/CellLabels.h"          // J8
#include "geometry/SolidSdf2D.h"       // J8
#include "numerics/MacTransfer.h"
#include "particles/ParticleSet2D.h"
#include "solvers/ISolver2D.h"
#include "diagnostics/ParticleDiagnostics.h"
#include "numerics/SandFriction.h"     // J13

enum class TransferScheme { PIC, FLIP, APIC };   // J9 : FLIP (mélange réglé par flipRatio) ; J10 : APIC

// J13 : Water = le solveur des jours 8 à 11, inchangé ; Sand = étapes 5b et 7b en plus.
enum class Material { Water, Sand };

enum class ParticleScene {
    BoxVortex,      // boîte pleine, tourbillon cellulaire initial
    BoxAtRest,      // boîte pleine au repos : à remuer avec la souris
    RotatingDisk,   // disque en rotation rigide, TRANSFERTS SEULS (mesure de PIC ; exact en APIC)
    TankAtRest,     // cuve remplie jusqu'à tankHeight, avec gravité
    DamBreak,       // J9 : colonne d'eau damWidth × damHeight contre la paroi gauche, lâchée à t = 0
    RotatingDiskFree, // J11 : même disque, pas complet sans gravité (les particules bougent)
    ColumnCollapse    // J13 : colonne columnWidth × columnHeight posée au CENTRE du sol, lâchée à t = 0
};

struct FlipParams {
    // Grille : domaine [0, 1] × [0, ny/nx], h = 1/nx.
    int nx = 64, ny = 64;
    TransferScheme scheme = TransferScheme::PIC;

    // Particules (Zhu & Bridson #4.2.1).
    int      particlesPerAxis = 2;     // 2 -> 4 par cellule
    real     jitter = 1;               // 0 : réseau régulier, 1 : stratifié aléatoire
    unsigned seed = 1;
    real     density = 1;

    // Pas de temps. Les particules sont déplacées en sous-pas de CFL <= substepCfl
    // dans la vitesse de grille du pas (#4.2.5).
    real dt = real(1) / 60;
    bool autoDt = false;
    real cflTarget = 2;
    real substepCfl = 1;

    // J9 : le paramètre alpha pour mélanger FLIP et PIC (le schéma est choisi par `scheme`, plus haut).
    // Ignoré en APIC (J10), qui ne se mélange pas avec FLIP.
    real flipRatio = real(0.95);   // alpha, part de FLIP (0 : PIC, 1 : FLIP pur)

    // Scènes.
    real vortexAmplitude = 1;          // max |u| du tourbillon (m/s)
    real diskRadius = real(0.3);       // fraction de Lx (APIC #8 : rayon 0,3 dans [0,1]^2)
    real diskAngularVelocity = 2;      // \omega (rad/s)
    real tankHeight = real(0.5);       // hauteur d'eau, en fraction de Ly
    //      (0.5 avec ny = 64 : 32 cellules, surface alignée)
    real damWidth = real(0.25);       // J9 : largeur de la colonne d'eau, en fraction de Lx
    real damHeight = real(0.5);        // J9 : hauteur de la colonne, en fraction de Ly (rapport 2, Martin & Moyce)

    // Sable (Zhu & Bridson §3).                                              // J13
    Material material = Material::Water;
    real frictionAngleDeg = 30;        // φ (degrés), éq. 1 et 2
    real wallFriction = real(0.6);     // μ sable-paroi, éq. 4 (fig. 2 de l'article)
    bool rigidStress = true;           // éq. 3 : σ_rigid dans les cellules sous le seuil (SandFriction.h)
    bool rigidGroups = true;           // étape 5 : groupes rigides projetés sur un mouvement rigide
    bool reprojectAfterFriction = true;  // 2e projection après le frottement (volume conservé)
    real columnWidth = real(0.2);      // largeur de la colonne, en fraction de Lx
    real columnHeight = real(0.4);     // hauteur de la colonne, en fraction de Ly (rapport 2)

    // Surface libre.                                                       // J8
    Vec2 gravity = { 0, real(-9.81) };
    real pushMargin = real(1e-3);      // distance minimale aux parois, en fraction de h (celle du J7)
    bool tankObstacle = false;         // marche posée au fond de la cuve (facultatif)
    int  maxBubbleCells = 8;           // poches d'air isolées comblées (0 : jamais, pour voir l'effet)

    // Souris et projection.
    real brushRadius = real(0.05);
    real brushGain = 25;
    real projectionTol = 1e-6;
};

struct FlipBrush {
    bool active = false;
    Vec2 pos;
    Vec2 velocity;       // vitesse du curseur (m/s)
};

struct FlipStepStats {
    real dt = 0;
    real cfl = 0;                  // max|u_p| dt / h
    int  substeps = 0;
    int  pressureIterations = 0;
    real divBefore = 0, divAfter = 0;
    int  emptyFaces = 0;           // faces sans particule après P2G
    int  fluidCells = 0;           // J8
    int  extrapolationLayers = 0;  // J8
    real energyParticlesBefore = 0;  // avant P2G
    real energyGrid = 0;             // juste après P2G (avant forces et projection)
    real energyParticlesAfter = 0;   // après G2P
    real angularMomentumGrid = 0;    // J11 : L de la grille juste après P2G, autour du centre de masse
    int  frictionFaces = 0;          // J13 : faces freinées par l'éq. 4
    int  stressCells = 0;            // J13 : cellules où σ_f ≠ 0
    int  rigidCells = 0;             // J13 : parmi elles, celles sous le seuil (σ_rigid, éq. 3)
    int  rigidGroups = 0;            // J13 : groupes rigides (étape 5)
    real energyBeforeFriction = 0;   // J13 : énergie de la grille avant / après l'éq. 2
    real energyAfterFriction = 0;
    double msLabels = 0, msP2G = 0, msForces = 0, msProjection = 0,
        msExtrapolation = 0, msG2P = 0, msAdvection = 0, msFriction = 0;   // msFriction : J13
    double msTotal() const {
        return msLabels + msP2G + msForces + msProjection + msExtrapolation + msG2P + msAdvection
            + msFriction;
    }
};

class FlipSolver final : public ISolver2D {
public:
    explicit FlipSolver(const FlipParams& params = {}, ParticleScene scene = ParticleScene::BoxVortex);

    // --- ISolver2D -----------------------------------------------------------
    std::string name() const override {
        switch (params_.scheme) {
        case TransferScheme::PIC:  return "PIC (particules + grille MAC, surface libre)";
        case TransferScheme::FLIP: return "FLIP/PIC (particules + grille MAC, surface libre)";
        case TransferScheme::APIC: return "APIC (particules + grille MAC, surface libre)";
        }
        return "";
    }
    void reset() override;
    void step(real dt) override;
    real stableDt() const override;
    real time() const override { return time_; }
    std::vector<Diagnostic> diagnostics() const override;
    std::vector<FieldView> fields() const override;

    // --- Pilotage propre à ce solveur ----------------------------------------
    FlipParams& params() { return params_; }
    const FlipParams& params() const { return params_; }
    ParticleScene scene() const { return scene_; }
    void setScene(ParticleScene s) { scene_ = s; reset(); }
    void setBrush(const FlipBrush& b) { brush_ = b; }
    bool transfersOnly() const { return scene_ == ParticleScene::RotatingDisk; } //Nous dis si pas complet ou pas
    bool usesGravity()   const {
        return scene_ == ParticleScene::TankAtRest || scene_ == ParticleScene::DamBreak
            || scene_ == ParticleScene::ColumnCollapse;
    } //Nous dis si gravité ou pas (J13 : + colonne)
    bool isSand() const { return params_.material == Material::Sand; }                        // J13
    // J11 : de l'air autour du fluide (p = 0 sur sa surface). Les boîtes pleines n'en ont pas,
    // et le disque en transferts seuls ne fait pas de projection.
    bool hasFreeSurface() const { return usesGravity() || scene_ == ParticleScene::RotatingDiskFree; }
    Vec2 activeGravity() const { return usesGravity() ? params_.gravity : Vec2{ 0, 0 }; }
    real potentialEnergy() const;                                                         //E_p des particules
    const NoiseStats& noise() const { return noiseStats_; }                               //bruit

    real chosenDt() const;
    Vec2 domainCenter() const { return 0.5 * grid_.size(); }
    Vec2 centerOfMass() const { return particleCenterOfMass(particles_); }                // J11
    real initialRadiusOfGyration() const { return gyration0_; }                           // J11 : R_g à t = 0

    const ParticleSet2D& particles() const { return particles_; }
    ParticleSet2D& particles() { return particles_; }   // pour les tests
    const MacGrid2D& grid() const { return grid_; }
    const MacMasses& masses() const { return masses_; }
    const CellLabels& labels() const { return labels_; }          // J8
    const SolidSdf2D& solids() const { return solids_; }          // J8
    const FlipStepStats& lastStats() const { return stats_; }
    long                 stepCount() const { return steps_; }

private:
    void updateDerivedFields();

    FlipParams params_;
    ParticleScene scene_;
    MacGrid2D grid_;
    MacMasses masses_;
    Field2D uOld_, vOld_;                    // vitesse de grille avant les forces (pour FLIP, jour 9)
    ParticleSet2D particles_;
    CellLabels labels_;
    SolidSdf2D solids_{ Vec2{ 1, 1 } };          //redimensionné dans reset()
    Field2D speed_, vorticity_, divergence_, count_;
    Field2D labelsField_;                    // étiquettes converties pour l'affichage
    Field2D noise_;                          // variance des vitesses dans chaque cellule
    NoiseStats noiseStats_;                  // sa moyenne pondérée par la masse
    FlipBrush brush_;
    FlipStepStats stats_;
    real time_ = 0;
    long steps_ = 0;
    real gyration0_ = 0;                     // J11 : rayon de giration à la réinitialisation
};