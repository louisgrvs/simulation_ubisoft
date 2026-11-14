#pragma once
// ============================================================================
//  Solveur hybride particules / grille MAC (Zhu - Bridson #4.2)
// ============================================================================
//
//  Un pas, en schéma PIC, dans une boîte entierement remplie de fluide :
//
//    particules ──P2G──> grille ──force──> ──projection──> grille ──G2P (PIC)──> particules
//                                                                  │
//    particules <──────────── déplacement RK2 dans la vitesse de grille, en sous-pas 
//
//  L'advection, qui dissipait tant chez Stam, n'est plus qu'un déplacement de
//  points : aucune interpolation de la quantité transportée. La dissipation
//  restante vient des transferts eux-mêmes (G2P∘P2G est une moyenne), et c'est
//  ce que mesure la scène " disque en rotation ".
//
//  Limites assumées  :
//   - toute la boîte est du fluide : la projection du jour 4 s'applique telle
//     quelle. Étiquettes fluide/air, p = 0 dans l'air, gravité;
//   - seul PIC est disponible. FLIP aura besoin de la vitesse de grille
//     juste après P2G : elle est déjà gardée dans uOld_ / vOld_.
//
//  La scène RotatingDisk ne fait que les transferts (P2G puis G2P), sans force,
//  sans projection et sans déplacer les particules : on isole ainsi la
//  dissipation propre à PIC.
// ============================================================================
#include <string>
#include <vector>
#include "grid/MacGrid2D.h"
#include "numerics/MacTransfer.h"
#include "particles/ParticleSet2D.h"
#include "solvers/ISolver2D.h"

enum class TransferScheme { PIC };   // FLIP : jour 9, APIC : jour 10

enum class ParticleScene {
    BoxVortex,      // boîte pleine, tourbillon cellulaire initial
    BoxAtRest,      // boîte pleine au repos : à remuer avec la souris
    RotatingDisk    // disque en rotation rigide, transferts seul (mesure de PIC)
};

struct FlipParams {
    // Grille : domaine [0, 1] × [0, ny/nx], h = 1/nx.
    int nx = 64, ny = 64;
    TransferScheme scheme = TransferScheme::PIC;

    // Particules (Zhu - Bridson #4.2.1).
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

    // Scènes.
    real vortexAmplitude = 1;          // max |u| du tourbillon (m/s)
    real diskRadius = real(0.3);       // fraction de Lx (APIC §8 : rayon 0,3 dans [0,1]^2)
    real diskAngularVelocity = 2;      // \omega (rad/s)

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
    real energyParticlesBefore = 0;  // avant P2G
    real energyGrid = 0;             // juste après P2G (avant forces et projection)
    real energyParticlesAfter = 0;   // après G2P
    double msP2G = 0, msForces = 0, msProjection = 0, msG2P = 0, msAdvection = 0;
    double msTotal() const { return msP2G + msForces + msProjection + msG2P + msAdvection; }
};

class FlipSolver final : public ISolver2D {
public:
    explicit FlipSolver(const FlipParams& params = {}, ParticleScene scene = ParticleScene::BoxVortex);

    // --- ISolver2D -----------------------------------------------------------
    std::string name() const override { return "PIC (particules + grille MAC)"; }
    void reset() override;
    void step(real dt) override;
    real stableDt() const override;
    real time() const override { return time_; }
    std::vector<Diagnostic> diagnostics() const override;
    std::vector<FieldView> fields() const override;

    // --- Pilotage propre à ce solveur ----------------------------------------
    FlipParams&       params()       { return params_; }
    const FlipParams& params() const { return params_; }
    ParticleScene scene() const { return scene_; }
    void setScene(ParticleScene s) { scene_ = s; reset(); }
    void setBrush(const FlipBrush& b) { brush_ = b; }
    bool transfersOnly() const { return scene_ == ParticleScene::RotatingDisk; }

    real chosenDt() const;
    Vec2 domainCenter() const { return 0.5 * grid_.size(); }

    const ParticleSet2D& particles() const { return particles_; }
    ParticleSet2D&       particles()       { return particles_; }   // pour les tests
    const MacGrid2D&     grid() const { return grid_; }
    const MacMasses&     masses() const { return masses_; }
    const FlipStepStats& lastStats() const { return stats_; }
    long                 stepCount() const { return steps_; }

private:
    void advectParticles(real dt);
    void updateDerivedFields();

    FlipParams params_;
    ParticleScene scene_;
    MacGrid2D grid_;
    MacMasses masses_;
    Field2D uOld_, vOld_;                    // vitesse de grille après P2G (pour FLIP, jour 9)
    ParticleSet2D particles_;
    Field2D speed_, vorticity_, divergence_, count_;
    FlipBrush brush_;
    FlipStepStats stats_;
    real time_ = 0;
    long steps_ = 0;
};
