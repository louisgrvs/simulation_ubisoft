// Démo console (sans OpenGL) du jour 3 :
// disque de Zalesak transporté par une rotation solide, un tour complet,
// pour Euler et RK2, à plusieurs nombres de CFL.
//
// Sorties dans <dossier> :
//   <run>.csv             t, max, min, masse, écart du barycentre (en h), erreur L1 à la solution exacte
//   <run>_<k>.pgm         instantanés à 0, 1/4, 1/2 et 1 tour (images en niveaux de gris)
//   summary.csv           une ligne par run
//
// Usage : advection_demo [n=128] [dossier=out_advection]

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "grid/MacGrid2D.h"
#include "numerics/Advection.h"
#include "scenes/RotationScene.h"

namespace fs = std::filesystem;
using namespace rotation_scene;

// PGM binaire (P5), y vers le haut : la ligne 0 de l'image est j = ny-1.
static void writePGM(const fs::path& file, const Field2D& q) {
    std::ofstream out(file, std::ios::binary);
    out << "P5\n" << q.nx() << " " << q.ny() << "\n255\n";
    for (int j = q.ny() - 1; j >= 0; --j)
        for (int i = 0; i < q.nx(); ++i) {
            real a = std::clamp(q(i, j), real(0), real(1));
            out.put(char(static_cast<unsigned char>(a * 255.0 + 0.5)));
        }
}

int main(int argc, char** argv) {
    const int n = argc > 1 ? std::stoi(argv[1]) : 128;
    const fs::path dir = argc > 2 ? argv[2] : "out_advection";
    fs::create_directories(dir);

    const real L = 1.0, h = L / n;
    const real omega = 2 * kPi, T = 2 * kPi / omega;  // un tour complet

    MacGrid2D grid(n, n, h);
    setRotation(grid, omega);

    Field2D q0 = grid.makeCellField();
    setSlottedDisk(q0, L);
    const real m0 = mass(q0);
    const Vec2 c0 = centroid(q0);

    std::printf("Grille %dx%d, h = %.4g, |div u|max = %.2e\n", n, n, h, grid.maxAbsDivergence());
    std::printf("%-12s %6s %7s %9s %9s %9s %9s %10s %8s\n",
                "run", "pas", "CFL", "erreurL1", "masse%", "max", "min", "barycentre", "ms");

    std::ofstream summary(dir / "summary.csv");
    summary << "run,integrator,cfl,steps,l1,mass_drift,max,min,centroid_err_h,ms\n";

    struct Run { Integrator scheme; const char* name; real cfl; };
    const std::vector<Run> runs = {
        {Integrator::Euler, "euler", 0.5}, {Integrator::Euler, "euler", 2.0}, {Integrator::Euler, "euler", 8.0},
        {Integrator::RK2,   "rk2",   0.5}, {Integrator::RK2,   "rk2",   2.0}, {Integrator::RK2,   "rk2",   8.0},
    };

    for (const Run& r : runs) {
        // Nombre de pas multiple de 4 pour tomber pile sur 1/4, 1/2 et 1 tour.
        int steps = int(std::ceil(T / cflTimeStep(grid, r.cfl)));
        steps = 4 * ((steps + 3) / 4);
        const real dt = T / steps;
        const real cfl = grid.maxSpeed() * dt / h;

        char tag[64];
        std::snprintf(tag, sizeof tag, "%s_cfl%g", r.name, r.cfl);
        std::ofstream csv(dir / (std::string(tag) + ".csv"));
        csv << "t,max,min,mass_drift,centroid_err_h,l1\n";

        Field2D q = q0, tmp = grid.makeCellField(), exact = grid.makeCellField();
        const Vec2 center = 0.5 * grid.size();
        auto log = [&](real t) {
            setSlottedDisk(exact, L, omega * t);                    // solution exacte à l'instant t
            Vec2 cExact = center + rotate(c0 - center, omega * t);  // barycentre exact
            csv << t << "," << q.maxValue() << "," << q.minValue() << "," << (mass(q) - m0) / m0 << ","
                << norm(centroid(q) - cExact) / h << "," << l1Error(q, exact) << "\n";
        };
        writePGM(dir / (std::string(tag) + "_0.pgm"), q);
        log(0);

        double ms = 0;  // temps de l'advection seule (sans les diagnostics)
        for (int s = 1; s <= steps; ++s) {
            auto t0 = std::chrono::steady_clock::now();
            advectScalar(grid, q, tmp, dt, r.scheme);
            ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            std::swap(q, tmp);
            log(s * dt);
            if (s % (steps / 4) == 0 && s / (steps / 4) != 3)
                writePGM(dir / (std::string(tag) + "_" + std::to_string(s / (steps / 4)) + ".pgm"), q);
        }

        real l1 = l1Error(q, q0), drift = (mass(q) - m0) / m0, cerr = norm(centroid(q) - c0) / h;
        std::printf("%-12s %6d %7.2f %9.4f %+8.2f%% %9.4f %9.4f %9.2fh %8.0f\n",
                    tag, steps, cfl, l1, 100 * drift, q.maxValue(), q.minValue(), cerr, ms);
        summary << tag << "," << r.name << "," << cfl << "," << steps << "," << l1 << "," << drift
                << "," << q.maxValue() << "," << q.minValue() << "," << cerr << "," << ms << "\n";
    }
    std::printf("Sorties écrites dans %s\n", dir.string().c_str());
    return 0;
}
