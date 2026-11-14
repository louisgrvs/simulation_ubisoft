#include "solvers/FlipSolver.h"

#include <algorithm>
#include <chrono>
#include <cmath>

#include "diagnostics/FluidDiagnostics.h"
#include "diagnostics/ParticleDiagnostics.h"
#include "numerics/Forces.h"
#include "numerics/MacOperators.h"
#include "numerics/Projection.h"
#include "scenes/Seeding.h"

namespace {

using Clock = std::chrono::steady_clock;
double msSince(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

} // namespace

FlipSolver::FlipSolver(const FlipParams& params, ParticleScene scene)
    : params_(params), scene_(scene), grid_(params.nx, params.ny, real(1) / params.nx) {
    reset();
}

void FlipSolver::reset() {
    params_.nx = std::max(params_.nx, 4);
    params_.ny = std::max(params_.ny, 4);
    params_.particlesPerAxis = std::clamp(params_.particlesPerAxis, 1, 4);
    params_.jitter = std::clamp(params_.jitter, real(0), real(1));
    const int nx = params_.nx, ny = params_.ny;
    grid_ = MacGrid2D(nx, ny, real(1) / nx);
    masses_ = MacMasses{grid_.u, grid_.v};
    uOld_ = grid_.u;
    vOld_ = grid_.v;
    speed_      = grid_.makeCellField();
    vorticity_  = grid_.makeCellField();
    divergence_ = grid_.makeCellField();
    count_      = grid_.makeCellField();

    SeedOptions opt;
    opt.perAxis = params_.particlesPerAxis;
    opt.jitter  = params_.jitter;
    opt.seed    = params_.seed;
    opt.density = params_.density;

    const real h = grid_.h();
    const Vec2 L = grid_.size();
    particles_.clear();
    particles_.reserve(std::size_t(nx) * ny * opt.perAxis * opt.perAxis);

    switch (scene_) {
    case ParticleScene::BoxVortex: {
        // Tourbillon cellulaire continu, \psi  = A (Ly/\pi ) sin(\pi x/Lx) sin(\pi y/Ly) :
        // u = \partial \psi /\partial y, v = -\partial \psi /\partial x (même champ qu'au jour 5, ici évalué aux particules).
        const real A = params_.vortexAmplitude;
        auto vel = [=](Vec2 x) {
            const real ax = kPi * x.x / L.x, ay = kPi * x.y / L.y;
            return Vec2{A * std::sin(ax) * std::cos(ay), -A * (L.y / L.x) * std::cos(ax) * std::sin(ay)};
        };
        seedParticles(particles_, nx, ny, h, region::everywhere(), opt, vel);
        break;
    }
    case ParticleScene::BoxAtRest: //champs qui bouge pas qu'on remu avec la souris
        seedParticles(particles_, nx, ny, h, region::everywhere(), opt);
        break;

    case ParticleScene::RotatingDisk: {
        const Vec2 c = domainCenter();
        const real w = params_.diskAngularVelocity;
        auto vel = [=](Vec2 x) { const Vec2 r = x - c; return Vec2{-w * r.y, w * r.x}; };   // ω × r
        seedParticles(particles_, nx, ny, h, region::disk(c, params_.diskRadius * L.x), opt, vel);
        break;
    }
    }
    brush_ = {};
    stats_ = {};
    time_ = 0;
    steps_ = 0;
    updateDerivedFields();
}

real FlipSolver::stableDt() const {
    const real vmax = particleMaxSpeed(particles_);
    const real dtCfl = vmax > 0 ? params_.cflTarget * grid_.h() / vmax : real(0.1);
    return std::min(dtCfl, real(0.1));
}

real FlipSolver::chosenDt() const {
    return params_.autoDt ? std::min(params_.dt, stableDt()) : params_.dt;
}

void FlipSolver::step(real dt) {
    //On suit la résolution comme dans le zhu-brigston, décrit en 4.2 :PIC methods
    if (!(dt > 0)) return;
    stats_ = {};
    stats_.dt = dt;
    stats_.cfl = particleMaxSpeed(particles_) * dt / grid_.h();
    stats_.energyParticlesBefore = particleEnergy(particles_);

    // 1. P2G (#4.2.2 zhu brigston : transfering to the grid)
    auto t0 = Clock::now();
    particlesToMac(particles_, grid_, masses_);
    stats_.energyGrid = gridEnergy(grid_, masses_);
    stats_.emptyFaces = masses_.emptyU + masses_.emptyV;
    stats_.msP2G = msSince(t0);

    if (!transfersOnly()) {
        // 2. Forces sur la grille. u·n = 0 d'abord : P2G a pu écrire sur les parois.
        t0 = Clock::now();
        enforceWallBC(grid_);
        uOld_ = grid_.u;     // FLIP : Delta u = u_nouveau - u_ancien
        vOld_ = grid_.v;
        if (brush_.active)
            addBrushForce(grid_, brush_.pos, params_.brushRadius, params_.brushGain * brush_.velocity, dt);
        stats_.msForces = msSince(t0);

        // 3. Projection : toute la boîte est fluide.
        t0 = Clock::now();
        ProjectionOptions opt;
        opt.tol = params_.projectionTol;
        opt.density = params_.density;
        const ProjectionStats ps = project(grid_, dt, opt);
        stats_.pressureIterations = ps.iterations;
        stats_.divBefore = ps.divBefore;
        stats_.divAfter = ps.divAfter;
        stats_.msProjection = msSince(t0);
    }

    // 4. G2P, schéma PIC (#4.2.4)
    t0 = Clock::now();
    macToParticlesPIC(grid_, particles_);
    stats_.energyParticlesAfter = particleEnergy(particles_);
    stats_.msG2P = msSince(t0);

    // 5. Déplacement des particules (#4.2.5)
    if (!transfersOnly()) {
        t0 = Clock::now();
        advectParticles(dt);
        stats_.msAdvection = msSince(t0);
    }

    time_ += dt;
    ++steps_;
    updateDerivedFields();
}

// RK2 (point milieu) dans la vitesse de grille, en sous-pas de CFL <= substepCfl.
// Le champ de grille est figé pendant le pas, comme dans l'article.
void FlipSolver::advectParticles(real dt) {
    const real h = grid_.h();
    const real vmax = grid_.maxSpeed();
    const int n = std::max(1, int(std::ceil(vmax * dt / (params_.substepCfl * h))));
    const real sub = dt / n;
    stats_.substeps = n;

    // Les particules restent dans le domaine, à une petite marge des parois.
    const Vec2 L = grid_.size();
    const real margin = real(1e-3) * h;
    for (std::size_t k = 0; k < particles_.size(); ++k) {
        Vec2 x = particles_.x[k];
        for (int s = 0; s < n; ++s) {
            const Vec2 xMid = x + (real(0.5) * sub) * interpolateMac(grid_, x);
            x = x + sub * interpolateMac(grid_, xMid);
            x.x = std::clamp(x.x, margin, L.x - margin);
            x.y = std::clamp(x.y, margin, L.y - margin);
        }
        particles_.x[k] = x;
    }
}

void FlipSolver::updateDerivedFields() {
    computeSpeed(grid_, speed_);
    computeVorticity(grid_, vorticity_);
    computeDivergence(grid_, divergence_);
    countParticlesPerCell(particles_, count_);
}

std::vector<Diagnostic> FlipSolver::diagnostics() const {
    const Vec2 P = particleMomentum(particles_);
    return {
        {"temps (s)",              time_},
        {"pas",                    double(steps_)},
        {"particules",             double(particles_.size())},
        {"dt (s)",                 stats_.dt},
        {"CFL",                    stats_.cfl},
        {"sous-pas",               double(stats_.substeps)},
        {"E particules",           particleEnergy(particles_)},
        {"E grille (apres P2G)",   stats_.energyGrid},
        {"E apres G2P/avant P2G", stats_.energyParticlesBefore > 0
                                       ? stats_.energyParticlesAfter / stats_.energyParticlesBefore : 1.0},
        {"L particules",           particleAngularMomentum(particles_, domainCenter())},
        {"|P| particules",         norm(P)},
        {"faces vides",            double(stats_.emptyFaces)},
        {"max|div| apres proj.",   stats_.divAfter},
        {"iter. CG pression",      double(stats_.pressureIterations)},
        {"ms par pas",             stats_.msTotal()},
    };
}

std::vector<FieldView> FlipSolver::fields() const {
    return {
        {"Norme de la vitesse",    &speed_,      false},
        {"Vorticite",              &vorticity_,  true},
        {"Particules par cellule", &count_,      false},
        {"Pression",               &grid_.p,     true},
        {"Divergence",             &divergence_, true},
    };
}
