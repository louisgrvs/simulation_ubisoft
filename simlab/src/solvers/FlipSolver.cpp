#include "solvers/FlipSolver.h"

#include <algorithm>
#include <chrono>
#include <cmath>

#include "diagnostics/FluidDiagnostics.h"
#include "diagnostics/FreeSurfaceDiagnostics.h"
#include "diagnostics/ParticleDiagnostics.h"
#include "numerics/Extrapolation.h"
#include "numerics/Forces.h"
#include "numerics/FreeSurface.h"
#include "numerics/MacOperators.h"
#include "numerics/ParticleAdvection.h"
#include "numerics/Projection.h"
#include "numerics/SandFriction.h"
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
    params_.tankHeight = std::clamp(params_.tankHeight, real(0), real(1));
    params_.flipRatio = std::clamp(params_.flipRatio, real(0), real(1));   // J9
    params_.damWidth = std::clamp(params_.damWidth, real(0), real(1));   // J9
    params_.damHeight = std::clamp(params_.damHeight, real(0), real(1));   // J9
    // J13 : \phi dans [0°, 90°], \mu >= 0, colonne dans le domaine.
    params_.frictionAngleDeg = std::clamp(params_.frictionAngleDeg, real(0), real(90));
    params_.wallFriction = std::max(params_.wallFriction, real(0));      // pas de borne haute
    params_.columnWidth = std::clamp(params_.columnWidth, real(0), real(1));    // fractions de Lx et Ly
    params_.columnHeight = std::clamp(params_.columnHeight, real(0), real(1));

    const int nx = params_.nx, ny = params_.ny;
    grid_ = MacGrid2D(nx, ny, real(1) / nx);
    masses_ = MacMasses{ grid_.u, grid_.v };
    uOld_ = grid_.u;
    vOld_ = grid_.v;
    speed_ = grid_.makeCellField();
    vorticity_ = grid_.makeCellField();
    divergence_ = grid_.makeCellField();
    count_ = grid_.makeCellField();
    noise_ = grid_.makeCellField();       // J9

    // surface libre.
    labels_ = CellLabels(nx, ny);
    solids_ = SolidSdf2D(grid_.size());
    if (scene_ == ParticleScene::TankAtRest && params_.tankObstacle) {
        const Vec2 Ls = grid_.size();
        solids_.addBox(Vec2{ real(0.4) * Ls.x, 0 }, Vec2{ real(0.6) * Ls.x, real(0.25) * Ls.y });   // marche au fond
    }
    labelsField_ = grid_.makeCellField();

    SeedOptions opt;
    opt.perAxis = params_.particlesPerAxis;
    opt.jitter = params_.jitter;
    opt.seed = params_.seed;
    opt.density = params_.density;

    const real h = grid_.h();
    const Vec2 L = grid_.size();
    particles_.clear();
    particles_.reserve(std::size_t(nx) * ny * opt.perAxis * opt.perAxis);

    switch (scene_) {
    case ParticleScene::BoxVortex: {
        // Tourbillon cellulaire continu, \psi = A (Ly/\pi) sin(\pi x/Lx) sin(\pi y/Ly) :
        // u = \partial \psi/\partial y, v = -\partial \psi/\partial x (même champ qu'au jour 5, ici évalué aux particules).
        const real A = params_.vortexAmplitude;
        auto vel = [=](Vec2 x) {
            const real ax = kPi * x.x / L.x, ay = kPi * x.y / L.y;
            return Vec2{ A * std::sin(ax) * std::cos(ay), -A * (L.y / L.x) * std::cos(ax) * std::sin(ay) };
            };
        seedParticles(particles_, nx, ny, h, region::everywhere(), opt, vel);
        break;
    }
    case ParticleScene::BoxAtRest:   // boîte pleine au repos, qu'on remue avec la souris
        seedParticles(particles_, nx, ny, h, region::everywhere(), opt);
        break;
    case ParticleScene::RotatingDiskFree:
    case ParticleScene::RotatingDisk: {
        const Vec2 c = domainCenter();
        const real w = params_.diskAngularVelocity;
        auto vel = [=](Vec2 x) { const Vec2 r = x - c; return Vec2{ -w * r.y, w * r.x }; };   // ω × r
        seedParticles(particles_, nx, ny, h, region::disk(c, params_.diskRadius * L.x), opt, vel);
        // en APIC, la partie affine exacte de la rotation (C = grad u, lignes cx et cy).
        // après l'ensemencement : avant, particles_ est vide et la boucle ne fait rien.
        // Sinon c reste à 0 (add() le met à 0), pour que PIC et FLIP ne soient pas modifiés.
        if (params_.scheme == TransferScheme::APIC)
            for (std::size_t k = 0; k < particles_.size(); ++k) {
                particles_.cx[k] = { 0, -w };   // (du/dx, du/dy) pour u = -w (y - c_y)
                particles_.cy[k] = { w, 0 };    // (dv/dx, dv/dy) pour v =  w (x - c_x)
            }
        break;
    }

    case ParticleScene::TankAtRest:  //  cuve remplie jusqu'à tankHeight·Ly, au repos
    {
        // L'eau remplit la cuve jusqu'à tankHeight·Ly, sauf dans l'obstacle éventuel.
        const Region water = region::box(Vec2{ 0, 0 }, Vec2{ L.x, params_.tankHeight * L.y });
        const SolidSdf2D sdf = solids_;
        seedParticles(particles_, nx, ny, h, [=](Vec2 x) { return water(x) && sdf.phi(x) > 0; }, opt);
        break;
    }
    case ParticleScene::DamBreak:    // colonne d'eau au repos contre la paroi gauche
    {
        const Region dam = region::box(Vec2{ 0, 0 }, Vec2{ params_.damWidth * L.x, params_.damHeight * L.y });
        const SolidSdf2D sdf = solids_;
        seedParticles(particles_, nx, ny, h, [=](Vec2 x) { return dam(x) && sdf.phi(x) > 0; }, opt);
        break;
    }
    case ParticleScene::ColumnCollapse:   // J13 : colonne au repos, posée au centre du sol
    {
        //  - région : region::box(lo, hi) avec lo = (Lx/2 − largeur/2, 0), hi = (Lx/2 + largeur/2, hauteur),
        //    largeur = columnWidth·Lx, hauteur = columnHeight·Ly
        //  - même ensemencement que DamBreak (avec le test sdf.phi(x) > 0, au cas où)
        //  - vitesse initiale nulle
        //  Même scène pour l'eau et le sable : seul params_.material change.
        const real w = params_.columnWidth * L.x;     // largeur en mètres (fraction × Lx)
        const real H = params_.columnHeight * L.y;    // hauteur en mètres
        const Vec2 lo{ real(0.5) * L.x - real(0.5) * w, 0 };
        const Vec2 hi{ real(0.5) * L.x + real(0.5) * w, H };
        const Region col = region::box(lo, hi);
        const SolidSdf2D sdf = solids_;
        seedParticles(particles_, nx, ny, h, [=](Vec2 x) { return col(x) && sdf.phi(x) > 0; }, opt);
        break;
    }
    }   // fin du switch (scene_)
    brush_ = {};
    stats_ = {};
    time_ = 0;
    steps_ = 0;
    gyration0_ = particleRadiusOfGyration(particles_, centerOfMass());   // J11
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

// énergie potentielle des particules, référence au fond (nulle sans gravité).
real FlipSolver::potentialEnergy() const {
    return particlePotentialEnergy(particles_, activeGravity());
}

void FlipSolver::step(real dt) {
    // Boucle de Zhu-Bridson #4.2, avec surface libre.
    if (!(dt > 0)) return;
    stats_ = {};
    stats_.dt = dt;
    stats_.cfl = particleMaxSpeed(particles_) * dt / grid_.h();
    stats_.energyParticlesBefore = particleEnergy(particles_);

    const Vec2 gravity = activeGravity();   // nulle hors des scènes à gravité
    // Vitesse maximale attendue pendant le pas (les particules + ce qu'ajoute la gravité).
    const real vmaxEst = particleMaxSpeed(particles_) + norm(gravity) * dt;
    const int layers = extrapolationLayerCount(vmaxEst, dt, grid_.h());
    stats_.extrapolationLayers = layers;

    // 1. Étiquettes Fluid / Air / Solid (#4.2.2)
    auto t0 = Clock::now();
    labelCells(labels_, particles_, grid_.h(), &solids_, params_.maxBubbleCells);
    stats_.fluidCells = labels_.count(CellType::Fluid);
    stats_.msLabels = msSince(t0);

    // 2. P2G (#4.2.2) ; en APIC, avec la partie affine c_p · (x_f - x_p) (J10)
    t0 = Clock::now();
    particlesToMac(particles_, grid_, masses_, params_.scheme == TransferScheme::APIC);   // terme affine en APIC
    stats_.energyGrid = gridEnergy(grid_, masses_);
    // J11 : la mesure honnête de L, la même pour les trois schémas. En APIC, une partie
    // du moment cinétique est portée par cx, cy : la somme des x_p × m_p v_p ne la voit pas,
    // la grille si (P2G l'y a déposée).
    stats_.angularMomentumGrid = gridAngularMomentum(grid_, masses_, centerOfMass());
    stats_.emptyFaces = masses_.emptyU + masses_.emptyV;
    stats_.msP2G = msSince(t0);
    uOld_ = grid_.u;   // en transferts seuls, la grille ne change pas (Delta u = 0)
    vOld_ = grid_.v;

    if (!transfersOnly()) {
        // 3. Extrapolation n°1 : remplit les faces qui n'ont reçu aucun poids.
        t0 = Clock::now();
        FaceMask uValid = maskFromWeights(masses_.u);
        FaceMask vValid = maskFromWeights(masses_.v);
        extrapolateVelocity(grid_, uValid, vValid, layers);
        enforceSolidBC(grid_, labels_);   // P2G et l'extrapolation ont pu écrire sur les parois

        // 4. Vitesse avant les forces, pour FLIP : après l'extrapolation, avant la gravité.
        uOld_ = grid_.u;
        vOld_ = grid_.v;

        // 5. Forces (#4.2.3)
        addGravity(grid_, gravity, dt);
        if (brush_.active)
            addBrushForce(grid_, brush_.pos, params_.brushRadius, params_.brushGain * brush_.velocity, dt);
        stats_.msForces = msSince(t0);

        // 5b. [Sable] frottement aux parois (éq. 4) : lit u·n AVANT que project() l'annule.  // J13
        if (isSand()) {
            t0 = Clock::now();
            // \mu : coefficient de frottement sable-paroi (éq. 4), 0,6 dans la fig. 2 de l'article.
            const WallFrictionStats st = applyWallFriction(grid_, labels_, params_.wallFriction);
            stats_.frictionFaces = st.frictionFaces;
            stats_.msFriction += msSince(t0);
        }

        // 6-7. Projection avec p = 0 dans l'air (enforceSolidBC est appelé au début de project).
        t0 = Clock::now();
        ProjectionOptions opt;
        opt.tol = params_.projectionTol;
        opt.density = params_.density;
        const ProjectionStats ps = project(grid_, labels_, dt, opt);
        stats_.pressureIterations = ps.iterations;
        stats_.divBefore = ps.divBefore;
        stats_.divAfter = ps.divAfter;
        stats_.msProjection = msSince(t0);

        // 7b. [Sable] contrainte de frottement (éq. 2), avec grid_.p de la projection.  // J13
        if (isSand()) {
            t0 = Clock::now();
            FrictionalStressOptions fo;   // pas "opt" : ce nom est déjà pris par la projection
            fo.density = params_.density;
            fo.frictionAngleDeg = params_.frictionAngleDeg;
            fo.rigidStress = params_.rigidStress;
            fo.rigidGroups = params_.rigidGroups;
            const FrictionalStressStats stressStat = applyFrictionalStress(grid_, labels_, dt, fo);
            stats_.stressCells = stressStat.stressCells;
            stats_.rigidCells = stressStat.rigidCells;
            stats_.rigidGroups = stressStat.rigidGroups;
            stats_.energyBeforeFriction = stressStat.energyBefore;
            stats_.energyAfterFriction = stressStat.energyAfter;
            // 7c. Deuxième projection : la contrainte et les groupes rigides ne sont pas à
            //     divergence nulle. Sans elle, le sable perd du volume à chaque pas (les
            //     particules se tassent : jusqu'à −65 % en 2 s pour φ = 45°).
            if (params_.reprojectAfterFriction) {
                const ProjectionStats ps2 = project(grid_, labels_, dt, opt);
                stats_.pressureIterations += ps2.iterations;
                stats_.divAfter = ps2.divAfter;
            }
            stats_.msFriction += msSince(t0);
        }

        // 8-9. Extrapolation n°2 depuis les faces fluides, puis parois à nouveau.
        t0 = Clock::now();
        uValid = uMaskFromLabels(labels_);
        vValid = vMaskFromLabels(labels_);
        extrapolateVelocity(grid_, uValid, vValid, layers);
        enforceSolidBC(grid_, labels_);
        stats_.msExtrapolation = msSince(t0);
    }

    // 10.G2P (#4.2.4) : PIC, mélange FLIP/PIC ou APIC 
    t0 = Clock::now();
    if (params_.scheme == TransferScheme::APIC) {
        macToParticlesAPIC(grid_, particles_);
    }
    else {
        if (params_.scheme == TransferScheme::FLIP)
            macToParticlesFLIP(grid_, uOld_, vOld_, particles_, params_.flipRatio);
        else
            macToParticlesPIC(grid_, particles_);
        // PIC et FLIP ne calculent pas c. On le remet à 0 pour qu'un retour
        // à chaud vers APIC (interface) ne reparte pas de gradients périmés.
        std::fill(particles_.cx.begin(), particles_.cx.end(), Vec2{ 0, 0 });
        std::fill(particles_.cy.begin(), particles_.cy.end(), Vec2{ 0, 0 });
    }
    stats_.energyParticlesAfter = particleEnergy(particles_);
    stats_.msG2P = msSince(t0);

    // 11. Déplacement RK2 en sous-pas + repoussement hors des parois (#4.2.5)
    if (!transfersOnly()) {
        t0 = Clock::now();
        const AdvectionStats adv = advectParticlesRK2(particles_, grid_, dt, params_.substepCfl,
            solids_, params_.pushMargin * grid_.h());
        stats_.substeps = adv.substeps;
        stats_.msAdvection = msSince(t0);
    }

    time_ += dt;
    ++steps_;
    updateDerivedFields();
}

void FlipSolver::updateDerivedFields() {
    computeSpeed(grid_, speed_);
    computeVorticity(grid_, vorticity_);
    computeDivergence(grid_, divergence_);
    countParticlesPerCell(particles_, count_);
    labels_.toField(labelsField_);   // 0 = air, 1 = fluide, 2 = solide
    noiseStats_ = cellVelocityVariance(particles_, noise_);
}

std::vector<Diagnostic> FlipSolver::diagnostics() const {
    const Vec2 P = particleMomentum(particles_);
    const Vec2 G = centerOfMass();                                                  // J11
    const real Rg = particleRadiusOfGyration(particles_, G);                        // J11
    const real h = grid_.h();
    const real Ec = particleEnergy(particles_), Ep = potentialEnergy();
    return {
        { "temps (s)",              time_ },
        { "pas",                    double(steps_) },
        { "particules",             double(particles_.size()) },
        { "dt (s)",                 stats_.dt },
        { "CFL",                    stats_.cfl },
        { "sous-pas",               double(stats_.substeps) },
        { "E particules",           Ec },
        { "E potentielle",          Ep },                                           // J9
        { "E totale (c + p)",       Ec + Ep },                                      // J9
        { "bruit <sigma2> (m2/s2)", noiseStats_.meanVariance },                     // J9
        { "bruit relatif",          noiseStats_.relative },                         // J9
        { "E grille (apres P2G)",   stats_.energyGrid },
        { "E apres G2P/avant P2G", stats_.energyParticlesBefore > 0
                                       ? stats_.energyParticlesAfter / stats_.energyParticlesBefore : 1.0 },
        { "L particules (/ x_G)",   particleAngularMomentum(particles_, G) },     // J11 : autour du centre de masse
        { "L grille (apres P2G)",   stats_.angularMomentumGrid },                   // J11
        { "centre de masse x",      G.x },                                          // J11
        { "centre de masse y",      G.y },                                          // J11
        { "rayon de giration R_g",  Rg },                                           // J11 : R/sqrt(2) pour un disque
        { "R_g / R_g(0)",           gyration0_ > 0 ? Rg / gyration0_ : 1.0 },       // J11
        { "|P| particules",         norm(P) },
        { "max|u| particules",      particleMaxSpeed(particles_) },                 // J8 : vitesses parasites
        { "cellules fluides",       double(stats_.fluidCells) },                    // J8
        { "aire fluide (m2)",       double(stats_.fluidCells) * h * h },            // J8 : volume en 2D
        { "couches extrapolees",    double(stats_.extrapolationLayers) },           // J8
        { "ecart a rho g h",        usesGravity()                                  // J8
                                       ? measureFreeSurface(grid_, labels_, particles_, params_.density,
                                                            norm(params_.gravity)).hydroErr
                                       : 0.0 },
        { "faces vides",            double(stats_.emptyFaces) },
        { "max|div| apres proj.",   stats_.divAfter },
        { "iter. CG pression",      double(stats_.pressureIterations) },
        { "faces freinees (eq. 4)", double(stats_.frictionFaces) },                // J13
        { "cellules sigma_f",       double(stats_.stressCells) },                  // J13
        { "cellules rigides",       double(stats_.rigidCells) },                   // J13
        { "groupes rigides",        double(stats_.rigidGroups) },                  // J13
        { "E apres/avant eq. 2",    stats_.energyBeforeFriction > 0                // J13 : doit rester <= 1
                                       ? stats_.energyAfterFriction / stats_.energyBeforeFriction : 1.0 },
        { "ms par pas",             stats_.msTotal() },
    };
}

std::vector<FieldView> FlipSolver::fields() const {
    return {
        { "Norme de la vitesse",    &speed_,       false },
        { "Vorticite",              &vorticity_,   true },
        { "Particules par cellule", &count_,       false },
        { "Pression",               &grid_.p,      true },
        { "Divergence",             &divergence_,  true },
        { "Etiquettes",             &labelsField_, false },   // J8
        { "Bruit (variance/cellule)", &noise_,     false },   // J9
    };
}