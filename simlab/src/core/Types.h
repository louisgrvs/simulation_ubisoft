#pragma once
// Types de base du projet (voir CONVENTIONS.md, section 8).
#include <glm/glm.hpp>

namespace sim {

using Real = float;                 // passer à double pour les tests de convergence
using Vec2 = glm::vec<2, Real>;
using Mat2 = glm::mat<2, 2, Real>;  // rangée par colonnes : M[colonne][ligne]

inline constexpr Real kPi = Real(3.14159265358979323846);

} // namespace sim
