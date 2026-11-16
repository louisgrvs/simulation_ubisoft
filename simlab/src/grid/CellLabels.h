#pragma once
#include <cstdint>
#include <vector>
#include "grid/Field2D.h"

enum class CellType : std::uint8_t { Air = 0, Fluid = 1, Solid = 2 };

class CellLabels {
public:
    CellLabels(int nx = 0, int ny = 0) : nx_(nx), ny_(ny), t_(std::size_t(nx)* ny, CellType::Air) {} //au debut on dit toute les cellules sont de l'air

    int nx() const { return nx_; }
    int ny() const { return ny_; }

    // Hors du domaine -> Solid : les murs de la cuve deviennent des conditions de Neumann.
    CellType operator()(int i, int j) const {
        if (i < 0 || i >= nx_ || j < 0 || j >= ny_) return CellType::Solid;
        return t_[i + nx_ * j];                      // indice = i + nx·j (CONVENTIONS.md)
    }
    void set(int i, int j, CellType t) { t_[i + nx_ * j] = t; }
    void fill(CellType t) { std::fill(t_.begin(), t_.end(), t); }

    bool isFluid(int i, int j) const { return (*this)(i, j) == CellType::Fluid; }
    bool isAir(int i, int j) const { return (*this)(i, j) == CellType::Air; }
    bool isSolid(int i, int j) const { return (*this)(i, j) == CellType::Solid; }

    int count(CellType t) const { return int(std::count(t_.begin(), t_.end(), t)); }

    // Pour l'affichage : f doit être un champ aux centres (grid.makeCellField()).
    void toField(Field2D& f) const {
        for (int j = 0; j < ny_; ++j)
            for (int i = 0; i < nx_; ++i) f(i, j) = real(int((*this)(i, j)));
    }

private:
    int nx_, ny_;
    std::vector<CellType> t_;
};