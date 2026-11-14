#pragma once
// ============================================================================
//  Particules 2D
// ============================================================================
//  Un tableau par attribut (" structure de tableaux ") : lisible, efficace en
//  cache quand une boucle ne lit qu'un attribut, et prêt à être copié tel quel
//  dans un tampon GPU.
//
//  position x_p, vitesse u_p, masse m_p.
//  À venir : partie affine c_p , F_E et J_P .
//  Tout nouvel attribut doit être ajouté dans add(), clear() et reserve().
// ============================================================================
#include <cstddef>
#include <vector>
#include "core/Vec2.h"

class ParticleSet2D {
public:
    std::vector<Vec2> x;   // positions (m)
    std::vector<Vec2> v;   // vitesses (m/s)
    std::vector<real> m;   // masses (kg par mètre d'épaisseur en 2D)

    std::size_t size() const { return x.size(); }
    bool empty() const { return x.empty(); }

    void clear() {
        x.clear();
        v.clear();
        m.clear();
    }
    void reserve(std::size_t n) {
        x.reserve(n);
        v.reserve(n);
        m.reserve(n);
    }
    void add(Vec2 position, Vec2 velocity, real mass) {
        x.push_back(position);
        v.push_back(velocity);
        m.push_back(mass);
    }
};
