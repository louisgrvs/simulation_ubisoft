#include "solvers/TestPatternSolver.h"

#include <algorithm>
#include <cmath>

namespace sim {

TestPatternSolver::TestPatternSolver(int nx, int ny, Real h)
    : x_(nx, ny, h, Stagger::Center),
      y_(nx, ny, h, Stagger::Center),
      marker_(nx, ny, h, Stagger::Center),
      wave_(nx, ny, h, Stagger::Center) {
    reset();
}

void TestPatternSolver::reset() {
    t_ = 0;
    const Real Lx = Real(x_.nx()) * x_.h();
    const Real Ly = Real(x_.ny()) * x_.h();

    x_.setFrom([](Vec2 p) { return p.x; });
    y_.setFrom([](Vec2 p) { return p.y; });

    const Vec2 c(Real(0.25) * Lx, Real(0.75) * Ly);
    const Real r = Real(0.15) * Ly;
    marker_.setFrom([&](Vec2 p) { return glm::length(p - c) < r ? Real(1) : Real(0); });

    updateWave();
}

void TestPatternSolver::step(Real dt) {
    t_ += dt;
    updateWave();
}

void TestPatternSolver::updateWave() {
    const Real Lx = Real(wave_.nx()) * wave_.h();
    const Real Ly = Real(wave_.ny()) * wave_.h();
    const Real k = std::cos(t_);
    wave_.setFrom([&](Vec2 p) {
        return std::sin(2 * kPi * p.x / Lx) * std::sin(kPi * p.y / Ly) * k;
    });
}

std::vector<Diagnostic> TestPatternSolver::diagnostics() const {
    double integral = 0.0;  // somme en double (CONVENTIONS.md, section 8)
    for (Real v : wave_.data()) integral += double(v);
    integral *= double(wave_.h()) * double(wave_.h());

    const auto [lo, hi] = std::minmax_element(wave_.data().begin(), wave_.data().end());
    return {
        {"t (s)", double(t_)},
        {"min onde", double(*lo)},
        {"max onde", double(*hi)},
        {"integrale onde", integral},  // doit rester ~0 par antisymétrie en x
    };
}

std::vector<FieldView> TestPatternSolver::fields() const {
    return {
        {"x", &x_, false},
        {"y", &y_, false},
        {"marqueur", &marker_, false},
        {"onde", &wave_, true},
    };
}

} // namespace sim
