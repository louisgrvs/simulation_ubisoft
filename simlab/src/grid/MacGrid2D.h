#pragma once
#include <functional>
#include "grid/Field2D.h"

// Grille MAC (Marker-And-Cell, Harlow & Welch 1965) décalée : nx x ny cellules de pas h.
// Domaine physique : [0, nx h] x [0, ny h]. Un beau schéma : 
//
//              v(i,j+1)
//          +------o------+
//          |             |
//   u(i,j) o    p(i,j)   o u(i+1,j)
//          |             |
//          +------o------+
//               v(i,j)
//
//  p : centres des cellules            nx     x ny       en ((i+1/2h, (j+1/2h)
//  u : milieux des faces verticales    (nx+1) x ny       en ( i    h, (j+1/2h)
//  v : milieux des faces horizontales  nx     x (ny+1)   en ((i+1/2h,  j    h)
//
// Intérêt : divergence et gradient deviennent des différences centrées compactes,
// sans mode parasite en damier (analogue DF d'une paire stable pour Stokes).
class MacGrid2D {
public:
    MacGrid2D(int nx, int ny, real h);

    int  nx() const { return nx_; } //Nombre de cases horizontale
    int  ny() const { return ny_; } //Nombre de cases verticales
    real h()  const { return h_; } //pas de notre grille
    Vec2 size() const { return {nx_ * h_, ny_ * h_}; }
    Vec2 cellCenter(int i, int j) const { return {(i + real(0.5)) * h_, (j + real(0.5)) * h_}; }

    Field2D u, v, p; //champs scalaire verticale, horizontale et au centre de pression

    // Vitesse en un point quelconque : chaque composante est interpolée sur
    // sa propre grille (décalages différents).
    Vec2 velocity(Vec2 x) const { return {u.sample(x), v.sample(x)}; }

    // Remplit u et v à partir d'un champ analytique : u est évalué aux faces u,
    // v aux faces v .
    void setVelocity(const std::function<Vec2(Vec2)>& f);

    // Plus grande composante de vitesse stockée : sert au pas de temps CFL.
    real maxSpeed() const;

    // Divergence discrète dans la cellule (i, j) :
    //   (u(i+1,j) - u(i,j) + v(i,j+1) - v(i,j)) / h
    real divergence(int i, int j) const;
    real maxAbsDivergence() const;

    // Champ scalaire vierge aux centres des cellules (même disposition que p).
    Field2D makeCellField(real init = 0) const { return Field2D(nx_, ny_, h_, 0.5, 0.5, init); }

private:
    int  nx_, ny_;
    real h_;
};
