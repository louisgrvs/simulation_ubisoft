#pragma once
#include <cassert>
#include <cstddef>
#include <vector>
#include "core/Vec2.h"

// Tableau 2D d'échantillons réguliers, de pas h, DÉCALÉ de (ox, oy) cellules.
//
//   L'échantillon (i, j) se trouve au point physique  ((i + ox) h, (j + oy) h).
//
//   (ox, oy) = (0.5, 0.5) : centres des cellules     (p, fumée, température…)
//   (ox, oy) = (0.0, 0.5) : faces verticales          (u)
//   (ox, oy) = (0.5, 0.0) : faces horizontales        (v)
//   (ox, oy) = (0.0, 0.0) : nœuds                     (MPM, plus tard)
//
// Toute la logique « où vit cette inconnue » est concentrée ici : MacGrid2D
// n'est qu'un assemblage de trois Field2D avec des décalages différents :
// - une grille pour la pression (le centre des cellules)
// - une grille pour la vélocité horizontale (faces verticales)
// - une grille pour la vélocité verticale (faces horizontales)
class Field2D {
public:
    Field2D() = default;
    Field2D(int nx, int ny, real h, real ox, real oy, real init = 0);

    int  nx() const { return nx_; }
    int  ny() const { return ny_; }
    real h()  const { return h_; }
    real ox() const { return ox_; }
    real oy() const { return oy_; }

    real&       operator()(int i, int j)       { return data_[idx(i, j)]; }
    const real& operator()(int i, int j) const { return data_[idx(i, j)]; }

    // Position physique de l'échantillon (i, j).
    Vec2 position(int i, int j) const { return {(i + ox_) * h_, (j + oy_) * h_}; }

    // Interpolation bilinéaire en un point physique quelconque.
    // Hors du rectangle des échantillons : on borne (clamp) -> valeur du bord.
    // Le résultat est une combinaison convexe de 4 voisins (poids >= 0, somme 1).
    real sample(Vec2 x) const;

    void fill(real value);
    bool sameLayout(const Field2D& o) const;

    real minValue() const;
    real maxValue() const;
    real sum() const;

    std::vector<real>&       data()       { return data_; }
    const std::vector<real>& data() const { return data_; }

private:
    std::size_t idx(int i, int j) const {
        assert(i >= 0 && i < nx_ && j >= 0 && j < ny_);
        return std::size_t(i) + std::size_t(nx_) * std::size_t(j);  // i rapide : ligne par ligne
    }

    int  nx_ = 0, ny_ = 0;
    real h_ = 1, ox_ = 0, oy_ = 0;
    std::vector<real> data_;
};
