// Démo sans OpenGL (jour 4) : projection d'un champ aléatoire pour plusieurs tailles de grille.
// Écrit projection_residuals.csv (historique du résidu de CG) et projection_summary.csv.
// Usage : projection_demo [dossier_de_sortie]
#include <chrono>
#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include "grid/MacGrid2D.h"
#include "numerics/MacOperators.h"
#include "numerics/Projection.h"

int main(int argc, char** argv) {
    const std::string out = argc > 1 ? argv[1] : ".";
    std::ofstream res(out + "/projection_residuals.csv"), sum(out + "/projection_summary.csv");
    res << "n,iteration,residual\n";
    sum << "n,iterations,div_before,div_after,time_ms\n";

    std::printf("%6s %10s %12s %12s %10s\n", "n", "iter. CG", "max|div| av.", "max|div| ap.", "temps ms");
    for (int n : {16, 32, 64, 128, 256}) {
        MacGrid2D g(n, n, 1.0 / n);
        std::mt19937 rng(42);
        std::uniform_real_distribution<real> U(-1, 1);
        for (real& a : g.u.data()) a = U(rng);
        for (real& a : g.v.data()) a = U(rng);

        ProjectionOptions opt;
        opt.tol = 1e-8;
        Vector history;
        auto t0 = std::chrono::steady_clock::now();
        ProjectionStats st = project(g, 0.01, opt, &history);
        auto t1 = std::chrono::steady_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        for (size_t k = 0; k < history.size(); ++k) res << n << ',' << k << ',' << history[k] << '\n';
        sum << n << ',' << st.iterations << ',' << st.divBefore << ',' << st.divAfter << ',' << ms << '\n';
        std::printf("%6d %10d %12.3e %12.3e %10.1f\n", n, st.iterations, st.divBefore, st.divAfter, ms);
    }
    return 0;
}
