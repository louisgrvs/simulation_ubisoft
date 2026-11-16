#pragma once
#include <vector>
#include "core/Vec2.h"
#include "particles/ParticleSet2D.h"

class SolidSdf2D {                    // \phi > 0 dans le domaine fluide, \phi < 0 dans le solide
public:
    explicit SolidSdf2D(Vec2 boxSize);
    void addBox(Vec2 min, Vec2 max);  // obstacle facultatif dans la cuve
    real phi(Vec2 x) const;           // min(x, Lx − x, y, Ly − y, obstacles…)
    Vec2 normal(Vec2 x) const;        // \nabla \phi, par différences centrées
private:
    struct Box { Vec2 min, max; };   // un rectangle : coin bas-gauche, coin haut-droit
    Vec2 size_;                       // dimensions de la cuve (Lx, Ly)
    std::vector<Box> obstacles_;      // liste des obstacles ajoutés par addBox
};
void pushOutOfSolids(ParticleSet2D& p, const SolidSdf2D& s, real margin);