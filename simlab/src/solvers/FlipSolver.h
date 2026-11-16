#pragma once
// ============================================================================
//  Solveur hybride particules / grille MAC (Zhu & Bridson #4.2)
//  jour 8 : surface libre ; jour 9 : FLIP et rupture de barrage
// ============================================================================
//
//  Un pas :
//
//    1. étiquettes Fluid / Air / Solid à partir des particules        (#4.2.2)
//    2. P2G                                                           (#4.2.2)
//    3. extrapolation n°1 (faces sans poids) + parois
//    4. uOld_, vOld_ = vitesse avant les forces (pour FLIP)
//    5. gravité + souris                                              (#4.2.3)
//    6-7. projection avec p = 0 dans l'air                            (#4.2.3)
//    8-9. extrapolation n°2 (depuis les faces fluides) + parois
//   10. G2P : PIC, ou mélange alpha FLIP + (1 - alpha) PIC (J9)        (#4.2.4)
//   11. déplacement RK2 en sous-pas CFL + repoussement hors des parois (#4.2.5)
//
//  La scène RotatingDisk ne fait QUE les transferts (étapes 2 et 10).
//  La gravité n'agit que dans TankAtRest et DamBreak (les scènes du J7 restent inchangées).
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

enum class TransferScheme { PIC, FLIP };   // J9 : FLIP (mélange réglé par flipRatio) ; APIC : jour 10

enum class ParticleScene {
    BoxVortex,      // boîte pleine, tourbillon cellulaire initial
    BoxAtRest,      // boîte pleine au repos : à remuer avec la souris
    RotatingDisk,   // disque en rotation rigide, TRANSFERTS SEULS (mesure de PIC)
    TankAtRest,     // cuve remplie jusqu'à tankHeight, avec gravité
    DamBreak        // J9 : colonne d'eau damWidth × damHeight contre la paroi gauche, lâchée à t = 0
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
    real flipRatio = real(0.95);   // alpha, part de FLIP (0 : PIC, 1 : FLIP pur)

    // Scènes.
    real vortexAmplitude = 1;          // max |u| du tourbillon (m/s)
    real diskRadius = real(0.3);       // fraction de Lx (APIC #8 : rayon 0,3 dans [0,1]^2)
    real diskAngularVelocity = 2;      // \omega (rad/s)
    real tankHeight = real(0.5);       // hauteur d'eau, en fraction de Ly
    //      (0.5 avec ny = 64 : 32 cellules, surface alignée)
    real damWidth = real(0.25);       // J9 : largeur de la colonne d'eau, en fraction de Lx
    real damHeight = real(0.5);        // J9 : hauteur de la colonne, en fraction de Ly (rapport 2, Martin & Moyce)

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
    double msLabels = 0, msP2G = 0, msForces = 0, msProjection = 0,
        msExtrapolation = 0, msG2P = 0, msAdvection = 0;
    double msTotal() const {
        return msLabels + msP2G + msForces + msProjection + msExtrapolation + msG2P + msAdvection;
    }
};

class FlipSolver final : public ISolver2D {
public:
    explicit FlipSolver(const FlipParams& params = {}, ParticleScene scene = ParticleScene::BoxVortex);

    // --- ISolver2D -----------------------------------------------------------
    std::string name() const override {   // J9 : le nom suit le schéma choisi
        return params_.scheme == TransferScheme::PIC ? "PIC (particules + grille MAC, surface libre)"
            : "FLIP/PIC (particules + grille MAC, surface libre)";
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
    bool transfersOnly() const { return scene_ == ParticleScene::RotatingDisk; }
    bool usesGravity()   const { return scene_ == ParticleScene::TankAtRest || scene_ == ParticleScene::DamBreak; }
    Vec2 activeGravity() const { return usesGravity() ? params_.gravity : Vec2{ 0, 0 }; }   // J9
    real potentialEnergy() const;                                                         // J9 : E_p des particules
    const NoiseStats& noise() const { return noiseStats_; }                               // J9 : bruit

    real chosenDt() const;
    Vec2 domainCenter() const { return 0.5 * grid_.size(); }

    const ParticleSet2D& particles() const { return particles_; }
    ParticleSet2D& particles() { return particles_; }   // pour les tests
    const MacGrid2D& grid() const { return grid_; }
    const MacMasses& masses() const { return masses_; }
    const CellLabels& labels() const { return labels_; }          // J8
    const SolidSdf2D& solids() const { return solids_; }          // J8
    const FlipStepStats& lastStats() const { return stats_; }
    long                 stepCount() const { return steps_; }

private:
    void updateDerivedFields();        // J8 : advectParticles a disparu (→ ParticleAdvection)

    FlipParams params_;
    ParticleScene scene_;
    MacGrid2D grid_;
    MacMasses masses_;
    Field2D uOld_, vOld_;                    // vitesse de grille avant les forces (pour FLIP, jour 9)
    ParticleSet2D particles_;
    CellLabels labels_;                      // J8
    SolidSdf2D solids_{ Vec2{ 1, 1 } };          // J8 : redimensionné dans reset()
    Field2D speed_, vorticity_, divergence_, count_;
    Field2D labelsField_;                    // J8 : étiquettes converties pour l'affichage
    Field2D noise_;                          // J9 : variance des vitesses dans chaque cellule
    NoiseStats noiseStats_;                  // J9 : sa moyenne pondérée par la masse
    FlipBrush brush_;
    FlipStepStats stats_;
    real time_ = 0;
    long steps_ = 0;
};