#pragma once
// ============================================================================
//  Fumée 2D « Stable Fluids » (Stam 1999) sur grille MAC
// ============================================================================
//
//  Un pas de temps, dans l'ordre de Stam (#2.2 et pseudo-code #3.2) :
//
//   Vstep  w0 ──force──> w1 ──auto-advection──> w2 ──diffusion ν──> w3 ──projection──> w4
//          (souris,       (semi-lagrangien RK2,  (implicite, CG)    (CG)
//           flottabilité)  jour 3)
//
//   Sstep  s0 ──source──> s1 ──advection par w4──> s2 ──diffusion κ──> s3 ──dissipation alpha──> s4
//
//  Chaque sous-étape est inconditionnellement stable : le pas de temps ne règle
//  que la précision, jamais la stabilité. C'est ce que vérifie le test
//  « rien n'explose à grand pas de temps ».
//
//  Écart assumé avec Stam : la grille MAC (vitesses sur les faces) au lieu de
//  sa grille colocalisée ; voir le jour 4.
// ============================================================================
#include <string>
#include <vector>
#include "grid/MacGrid2D.h"
#include "numerics/Advection.h"
#include "solvers/ISolver2D.h"

enum class SmokeScene {
    Plume,    // boîte vide + émetteur en bas : le panache classique
    Bubble,   // une bulle de fumée qui monte par flottabilité
    Vortex    // tourbillon cellulaire stationnaire : mesure de la dissipation
};

struct StableFluidsParams {
    // Grille : domaine [0, 1] × [0, ny/nx], cellules carrées h = 1/nx.
    int nx = 64, ny = 64;

    // Pas de temps (s). Avec autoDt, on prend min(dt, cflTarget · h / max|u|).
    real dt        = real(1) / 60;
    bool autoDt    = false;
    real cflTarget = 2;

    // Coefficients physiques (Stam §2.4), en unités SI du domaine.
    real viscosity   = 0;     // nu   (m²/s)   diffusion de la vitesse
    real diffusivity = 0;     // kappa  (m²/s)   diffusion de la fumée
    real dissipation = 0;     // alpha   (1/s)    disparition de la fumée
    real buoyancy    = 4;     // beta   (m/s² par unité de fumée), flottabilité

    // Émetteur fixe (scène Plume).
    bool emitterOn     = true;
    Vec2 emitterPos    = {0.5, 0.1};
    real emitterRadius = 0.05;
    real emitterRate   = 3;   // unités de fumée par seconde au centre

    // Pinceau de la souris.
    real brushRadius    = 0.05;
    real brushGain      = 25;   // f = gain · w(x) · V_souris   (1/s)
    real brushSmokeRate = 6;

    // Projection et advection.
    real projectionTol = 1e-6;  // max |div| toléré après projection (1/s)
    Integrator integrator = Integrator::RK2;
};

// Ce que la souris applique au prochain pas (coordonnées physiques).
struct Brush {
    bool active = false;
    Vec2 pos;
    Vec2 velocity;          // vitesse du curseur (m/s)
    bool push   = true;     // force proportionnelle à la vitesse du curseur
    bool inject = true;     // source de fumée
};

// Statistiques du dernier pas, pour le panneau et les mesures.
struct StepStats {
    real dt = 0;
    real cfl = 0;                  // max|u| dt / h au début du pas
    int  pressureIterations = 0;
    real divBefore = 0;            // max |div w3| avant projection
    real divAfter  = 0;            // max |div w4| après projection
    int  viscosityIterations = 0;
    int  smokeDiffusionIterations = 0;
    double msForces = 0, msAdvection = 0, msDiffusion = 0, msProjection = 0, msSmoke = 0;
    double msTotal() const { return msForces + msAdvection + msDiffusion + msProjection + msSmoke; }
};

class StableFluidsSolver final : public ISolver2D {
public:
    explicit StableFluidsSolver(const StableFluidsParams& params = {},
                                SmokeScene scene = SmokeScene::Plume);

    // --- ISolver2D -----------------------------------------------------------
    std::string name() const override { return "Stable Fluids (MAC)"; }
    void reset() override;
    void step(real dt) override;
    real stableDt() const override;
    real time() const override { return time_; }
    std::vector<Diagnostic> diagnostics() const override;
    std::vector<FieldView> fields() const override;

    // --- Pilotage propre à ce solveur ----------------------------------------
    StableFluidsParams&       params()       { return params_; }
    const StableFluidsParams& params() const { return params_; }
    SmokeScene scene() const { return scene_; }
    void setScene(SmokeScene s) { scene_ = s; reset(); }
    void setBrush(const Brush& b) { brush_ = b; }

    // Pas de temps qui sera réellement utilisé (dt fixe, ou CFL si autoDt).
    real chosenDt() const;

    const MacGrid2D& grid()  const { return grid_; }
    MacGrid2D&       grid()        { return grid_; }   // pour les tests
    const Field2D&   smoke() const { return smoke_; }
    Field2D&         smoke()       { return smoke_; }
    const StepStats& lastStats() const { return stats_; }
    long             stepCount() const { return steps_; }

    real kineticEnergy() const;
    real smokeAmount() const;

private:
    // Vstep
    void addForces(real dt);
    void advectVelocity(real dt);
    void diffuseVelocity(real dt);
    void projectVelocity(real dt);
    // Sstep
    void updateSmoke(real dt);
    // Champs affichés seulement (vitesse, vorticité, divergence)
    void updateDerivedFields();

    StableFluidsParams params_;
    SmokeScene scene_;
    MacGrid2D grid_;
    Field2D smoke_, smokeTmp_, uTmp_, vTmp_;
    Field2D speed_, vorticity_, divergence_;
    Brush brush_;
    StepStats stats_;
    real time_ = 0;
    long steps_ = 0;
};
