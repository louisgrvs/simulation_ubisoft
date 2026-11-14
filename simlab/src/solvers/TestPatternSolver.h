#pragma once
// Faux solveur de J2 : des champs analytiques pour valider l'orientation,
// les valeurs affichées et la boucle de l'application.
#include "solvers/ISolver2D.h"

namespace sim {

class TestPatternSolver final : public ISolver2D {
public:
    // Par défaut : domaine 2 m × 1 m, non carré exprès pour détecter
    // une inversion de nx et ny.
    explicit TestPatternSolver(int nx = 64, int ny = 32, Real h = Real(1) / Real(32));

    std::string name() const override { return "Motif de test"; }
    void reset() override;
    void step(Real dt) override;
    Real stableDt() const override { return Real(1) / Real(60); }
    Real time() const override { return t_; }
    std::vector<Diagnostic> diagnostics() const override;
    std::vector<FieldView> fields() const override;

private:
    void updateWave();

    Real t_ = 0;
    Field2D<Real> x_;       // f = x       : sombre à gauche, clair à droite
    Field2D<Real> y_;       // f = y       : sombre en bas, clair en haut
    Field2D<Real> marker_;  // disque en haut à gauche
    Field2D<Real> wave_;    // sin(2πx/Lx)·sin(πy/Ly)·cos(t), signé
};

} // namespace sim
