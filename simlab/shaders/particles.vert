#version 430 core
// Une particule = un sommet dessiné en GL_POINTS.
// Position physique (m) -> coordonnées normalisées [-1, 1] du viewport,
// lequel couvre exactement le domaine (même rectangle que le champ de fond).
layout(location = 0) in vec2  aPos;     // position physique
layout(location = 1) in float aValue;   // grandeur à colorer (norme de la vitesse)

uniform vec2  uDomain;      // taille du domaine (Lx, Ly) en mètres
uniform float uPointSize;   // diamètre en pixels

out float vValue;

void main()
{
    vValue = aValue;
    gl_Position = vec4(2.0 * aPos / uDomain - 1.0, 0.0, 1.0);   // y vers le haut, comme CONVENTIONS.md
    gl_PointSize = uPointSize;   // pris en compte grâce à glEnable(GL_PROGRAM_POINT_SIZE)
}
