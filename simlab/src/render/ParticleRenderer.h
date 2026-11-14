#pragma once
// ============================================================================
//  Affichage des particules en points (jour 7)
// ============================================================================
//  upload() : copie positions et valeurs (double -> float) dans UN tampon de
//             sommets entrelacé [x, y, valeur] par particule. Le tampon est
//             alloué une fois (glBufferData), puis seulement mis à jour
//             (glBufferSubData) tant que le nombre de particules ne dépasse
//             pas sa capacité : cf. LearnOpenGL, « Advanced Data ».
//  draw()   : dessine en GL_POINTS dans le rectangle donné, qui doit couvrir
//             exactement le domaine (le même que celui du champ de fond).
// ============================================================================
#include <string>
#include <vector>
#include "particles/ParticleSet2D.h"
#include "render/FieldRenderer.h"

class ParticleRenderer {
public:
    using Rect = FieldRenderer::Rect;

    explicit ParticleRenderer(const std::string& shaderDir);   // dossier de particles.vert/.frag
    ~ParticleRenderer();
    ParticleRenderer(const ParticleRenderer&) = delete;
    ParticleRenderer& operator=(const ParticleRenderer&) = delete;

    // values[k] : grandeur colorée de la particule k (ex. norme de la vitesse).
    // Si values est vide, la valeur est 0 pour toutes.
    void upload(const ParticleSet2D& particles, const std::vector<float>& values);

    // Palette « inferno » sur [vmin, vmax] (uniform = false), ou couleur unie.
    void draw(Rect viewport, Vec2 domainSize, float pointSizePx, float vmin, float vmax,
              bool uniform = false, const float rgb[3] = nullptr);

    std::size_t count() const { return count_; }

private:
    unsigned program_ = 0, vao_ = 0, vbo_ = 0;
    std::size_t capacity_ = 0, count_ = 0;     // en particules
    int locDomain_ = -1, locSize_ = -1, locVmin_ = -1, locVmax_ = -1, locMode_ = -1, locColor_ = -1;
    std::vector<float> staging_;
};
