#include "solvers/StableFluidsSolver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <utility>

#include "diagnostics/FluidDiagnostics.h"
#include "numerics/Diffusion.h"
#include "numerics/Forces.h"
#include "numerics/MacOperators.h"
#include "numerics/Projection.h"
#include "scenes/SmokeScenes.h"

namespace {

using Clock = std::chrono::steady_clock;
double msSince(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

} // namespace

StableFluidsSolver::StableFluidsSolver(const StableFluidsParams& params, SmokeScene scene)
    : params_(params), scene_(scene), grid_(params.nx, params.ny, real(1) / params.nx) {
    reset();
}

void StableFluidsSolver::reset() {
    const int nx = std::max(params_.nx, 4), ny = std::max(params_.ny, 4);
    params_.nx = nx;
    params_.ny = ny;
    grid_ = MacGrid2D(nx, ny, real(1) / nx);   // reconstruite : la taille a pu changer

    smoke_      = grid_.makeCellField();
    smokeTmp_   = grid_.makeCellField();
    speed_      = grid_.makeCellField();
    vorticity_  = grid_.makeCellField();
    divergence_ = grid_.makeCellField();
    uTmp_ = grid_.u;
    vTmp_ = grid_.v;

    const Vec2 L = grid_.size();
    switch (scene_) {
    case SmokeScene::Plume:
        break;                                         // tout part de l'émetteur
    case SmokeScene::Bubble:
        smoke_scene::addSmokeDisk(smoke_, {0.5 * L.x, 0.25 * L.y}, 0.12 * L.x, 1);
        break;
    case SmokeScene::Vortex:
        smoke_scene::setCellularVortex(grid_, 1);
        // Une bande de fumée pour voir le tourbillon l'enrouler.
        smoke_scene::addSmokeDisk(smoke_, {0.5 * L.x, 0.75 * L.y}, 0.1 * L.x, 1);
        break;
    }
    brush_ = {};
    stats_ = {};
    time_  = 0;
    steps_ = 0;
    updateDerivedFields();
}

real StableFluidsSolver::stableDt() const {
    // Le schéma est stable pour tout dt : ce pas vise seulement une précision
    // raisonnable (CFL cible). Plafonné pour le cas d'un fluide au repos.
    const real umax = grid_.maxSpeed();
    const real dtCfl = umax > 0 ? params_.cflTarget * grid_.h() / umax : real(0.1);
    return std::min(dtCfl, real(0.1));
}

real StableFluidsSolver::chosenDt() const {
    return params_.autoDt ? std::min(params_.dt, stableDt()) : params_.dt;
}

void StableFluidsSolver::step(real dt) {
    if (!(dt > 0)) return;
    stats_ = {};
    stats_.dt  = dt;
    stats_.cfl = grid_.maxSpeed() * dt / grid_.h();

    // --- Vstep ----------------------------------------------------------------
    auto t0 = Clock::now();
    addForces(dt);
    stats_.msForces = msSince(t0);

    t0 = Clock::now();
    advectVelocity(dt);
    stats_.msAdvection = msSince(t0);

    t0 = Clock::now();
    diffuseVelocity(dt);
    stats_.msDiffusion = msSince(t0);

    t0 = Clock::now();
    projectVelocity(dt);
    stats_.msProjection = msSince(t0);

    // --- Sstep ----------------------------------------------------------------
    t0 = Clock::now();
    updateSmoke(dt);
    stats_.msSmoke = msSince(t0);

    time_ += dt;
    ++steps_;
    updateDerivedFields();
}

// Étape 1 : w1 = w0 + dt f.
void StableFluidsSolver::addForces(real dt) {
    if (params_.buoyancy != 0) addBuoyancy(grid_, smoke_, params_.buoyancy, dt);

    if (brush_.active && brush_.push) {
        // f = gain · w(x) · V_souris : sur un pas, le fluide sous le pinceau
        // reçoit dt · gain de la vitesse du curseur (une fraction si dt·gain < 1).
        addBrushForce(grid_, brush_.pos, params_.brushRadius,
                      params_.brushGain * brush_.velocity, dt);
    }
    enforceWallBC(grid_);
}

// Étape 2 : auto-advection. On remonte les caractéristiques dans le champ w1
// lui-même (Stam : « the velocity field is moved along itself »). advectScalar
// lit w1 dans grid_ et écrit dans des tampons : u et v sont advectés par le
// MÊME w1, ce qui impose de ne pas écraser u avant d'avoir advecté v.
void StableFluidsSolver::advectVelocity(real dt) {
    advectScalar(grid_, grid_.u, uTmp_, dt, params_.integrator);
    advectScalar(grid_, grid_.v, vTmp_, dt, params_.integrator);
    std::swap(grid_.u, uTmp_);
    std::swap(grid_.v, vTmp_);
    enforceWallBC(grid_);
}

// Étape 3 : (I - ν dt Δ) w3 = w2.
void StableFluidsSolver::diffuseVelocity(real dt) {
    DiffusionStats st = ::diffuseVelocity(grid_, params_.viscosity, dt);
    stats_.viscosityIterations = st.iterations;
}

// Étape 4 : w4 = P(w3), jour 4.
void StableFluidsSolver::projectVelocity(real dt) {
    ProjectionOptions opt;
    opt.tol = params_.projectionTol;
    ProjectionStats st = project(grid_, dt, opt);
    stats_.pressureIterations = st.iterations;
    stats_.divBefore = st.divBefore;
    stats_.divAfter  = st.divAfter;
}

// Sstep (Stam §2.4) : source, transport par w4, diffusion κ, dissipation α.
void StableFluidsSolver::updateSmoke(real dt) {
    if (scene_ == SmokeScene::Plume && params_.emitterOn)
        addBrushSource(smoke_, params_.emitterPos, params_.emitterRadius, params_.emitterRate, dt);
    if (brush_.active && brush_.inject)
        addBrushSource(smoke_, brush_.pos, params_.brushRadius, params_.brushSmokeRate, dt);

    advectScalar(grid_, smoke_, smokeTmp_, dt, params_.integrator);
    std::swap(smoke_, smokeTmp_);

    stats_.smokeDiffusionIterations = diffuseCells(smoke_, params_.diffusivity, dt).iterations;
    dissipate(smoke_, params_.dissipation, dt);
}

void StableFluidsSolver::updateDerivedFields() {
    computeSpeed(grid_, speed_);
    computeVorticity(grid_, vorticity_);
    computeDivergence(grid_, divergence_);
}

real StableFluidsSolver::kineticEnergy() const { return ::kineticEnergy(grid_); }
real StableFluidsSolver::smokeAmount() const { return totalAmount(smoke_); }

std::vector<Diagnostic> StableFluidsSolver::diagnostics() const {
    return {
        {"temps (s)",              time_},
        {"pas",                    double(steps_)},
        {"dt (s)",                 stats_.dt},
        {"CFL",                    stats_.cfl},
        {"energie cinetique",      kineticEnergy()},
        {"masse de fumee",         smokeAmount()},
        {"max|div| avant proj.",   stats_.divBefore},
        {"max|div| apres proj.",   stats_.divAfter},
        {"iter. CG pression",      double(stats_.pressureIterations)},
        {"iter. CG viscosite",     double(stats_.viscosityIterations)},
        {"iter. CG diffusion",     double(stats_.smokeDiffusionIterations)},
        {"ms par pas",             stats_.msTotal()},
    };
}

std::vector<FieldView> StableFluidsSolver::fields() const {
    return {
        {"Fumee",              &smoke_,      false},
        {"Norme de la vitesse", &speed_,     false},
        {"Vorticite",          &vorticity_,  true},
        {"Pression",           &grid_.p,     true},
        {"Divergence",         &divergence_, true},
    };
}
