#pragma once
// Affiche un Field2D comme une image colorée (jour 2, adapté au Field2D du jour 3).
// C'est le seul code, avec main.cpp, qui parle à OpenGL.
//
//   upload() : copie le champ (double) dans une texture GL_R32F (float).
//              La texture n'est réallouée que si la taille change.
//   draw()   : dessine la texture dans un rectangle de l'écran, en gardant les
//              proportions du domaine (cellules carrées), et renvoie le rectangle
//              réellement occupé (pour convertir la souris en coordonnées physiques).
#include <string>
#include <vector>
#include "grid/Field2D.h"

class FieldRenderer {
public:
    enum class Palette { Sequential, Diverging };

    // Rectangle en pixels du framebuffer, origine en BAS à gauche (convention OpenGL).
    struct Rect { int x = 0, y = 0, w = 0, h = 0; };

    // shaderDir : dossier contenant field.vert et field.frag (avec le / final).
    explicit FieldRenderer(const std::string& shaderDir);
    ~FieldRenderer();
    FieldRenderer(const FieldRenderer&) = delete;
    FieldRenderer& operator=(const FieldRenderer&) = delete;

    void setSmooth(bool smooth);       // GL_LINEAR (vrai) ou GL_NEAREST (faux)
    void upload(const Field2D& field);
    Rect draw(Rect area, float vmin, float vmax, Palette palette);

private:
    unsigned program_ = 0, vao_ = 0, texture_ = 0;
    int texW_ = 0, texH_ = 0;
    bool smooth_ = true;
    int locVmin_ = -1, locVmax_ = -1, locPalette_ = -1, locField_ = -1;
    std::vector<float> staging_;       // conversion double -> float avant l'envoi
};
