#pragma once
#include <cmath>

// Réel utilisé côté CPU : double = référence précise (cf conventions.md).
// (Le GPU passera en float ; on comparera les deux plus tard.)
using real = double;

constexpr real kPi = 3.14159265358979323846;

// Petit vecteur 2D pour le coeur numérique : on utilise pas la bibliothèque de C++ (GLM) pcq les types sont en float.
struct Vec2 {
    real x = 0, y = 0;
};

inline Vec2 operator+(Vec2 a, Vec2 b) { return {a.x + b.x, a.y + b.y}; }
inline Vec2 operator-(Vec2 a, Vec2 b) { return {a.x - b.x, a.y - b.y}; }
inline Vec2 operator-(Vec2 a)         { return {-a.x, -a.y}; }
inline Vec2 operator*(real s, Vec2 a) { return {s * a.x, s * a.y}; }
inline Vec2 operator*(Vec2 a, real s) { return {s * a.x, s * a.y}; }
inline Vec2 operator/(Vec2 a, real s) { return {a.x / s, a.y / s}; }
inline Vec2& operator+=(Vec2& a, Vec2 b) { a.x += b.x; a.y += b.y; return a; }

inline real dot(Vec2 a, Vec2 b) { return a.x * b.x + a.y * b.y; }
inline real norm(Vec2 a)        { return std::sqrt(dot(a, a)); }

// Rotation d'angle theta (sens trigonométrique).
inline Vec2 rotate(Vec2 a, real theta) {
    real c = std::cos(theta), s = std::sin(theta);
    return {c * a.x - s * a.y, s * a.x + c * a.y};
}
