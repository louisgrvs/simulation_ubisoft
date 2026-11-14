#include "render/FieldRenderer.h"

#include <glad/glad.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

std::string readFile(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("FieldRenderer : impossible d'ouvrir " + path);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

unsigned compileShader(GLenum type, const std::string& source, const std::string& name) {
    const unsigned s = glCreateShader(type);
    const char* src = source.c_str();
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    int ok = 0;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[2048];
        glGetShaderInfoLog(s, GLsizei(sizeof log), nullptr, log);
        glDeleteShader(s);
        throw std::runtime_error("Compilation de " + name + " :\n" + log);
    }
    return s;
}

} // namespace

FieldRenderer::FieldRenderer(const std::string& shaderDir) {
    const unsigned vs = compileShader(GL_VERTEX_SHADER, readFile(shaderDir + "field.vert"), "field.vert");
    const unsigned fs = compileShader(GL_FRAGMENT_SHADER, readFile(shaderDir + "field.frag"), "field.frag");
    program_ = glCreateProgram();
    glAttachShader(program_, vs);
    glAttachShader(program_, fs);
    glLinkProgram(program_);
    glDeleteShader(vs);
    glDeleteShader(fs);
    int ok = 0;
    glGetProgramiv(program_, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[2048];
        glGetProgramInfoLog(program_, GLsizei(sizeof log), nullptr, log);
        throw std::runtime_error(std::string("Édition de liens du programme de champ :\n") + log);
    }
    locVmin_    = glGetUniformLocation(program_, "uVmin");
    locVmax_    = glGetUniformLocation(program_, "uVmax");
    locPalette_ = glGetUniformLocation(program_, "uPalette");
    locField_   = glGetUniformLocation(program_, "uField");

    // Le profil Core exige un VAO lié pour dessiner, même sans aucun attribut.
    glGenVertexArrays(1, &vao_);

    glGenTextures(1, &texture_);
    glBindTexture(GL_TEXTURE_2D, texture_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    setSmooth(true);
}

FieldRenderer::~FieldRenderer() {
    glDeleteTextures(1, &texture_);
    glDeleteVertexArrays(1, &vao_);
    glDeleteProgram(program_);
}

void FieldRenderer::setSmooth(bool smooth) {
    smooth_ = smooth;
    const GLint filter = smooth ? GL_LINEAR : GL_NEAREST;
    glBindTexture(GL_TEXTURE_2D, texture_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filter);
}

void FieldRenderer::upload(const Field2D& field) {
    const int w = field.nx(), h = field.ny();
    staging_.resize(field.data().size());
    // Même ordre mémoire des deux côtés (i rapide, j = 0 en bas) : simple copie.
    std::transform(field.data().begin(), field.data().end(), staging_.begin(),
                   [](real a) { return float(a); });

    glBindTexture(GL_TEXTURE_2D, texture_);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    if (w != texW_ || h != texH_) {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, staging_.data());
        texW_ = w;
        texH_ = h;
    } else {
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, GL_RED, GL_FLOAT, staging_.data());
    }
}

FieldRenderer::Rect FieldRenderer::draw(Rect area, float vmin, float vmax, Palette palette) {
    Rect r;
    if (texW_ == 0 || area.w <= 0 || area.h <= 0) return r;

    // Plus grand rectangle aux proportions du domaine, centré dans la zone.
    const float aspect = float(texW_) / float(texH_);
    r.w = area.w;
    r.h = int(float(area.w) / aspect);
    if (r.h > area.h) {
        r.h = area.h;
        r.w = int(float(area.h) * aspect);
    }
    r.x = area.x + (area.w - r.w) / 2;
    r.y = area.y + (area.h - r.h) / 2;

    glViewport(r.x, r.y, r.w, r.h);
    glUseProgram(program_);
    glUniform1f(locVmin_, vmin);
    glUniform1f(locVmax_, vmax);
    glUniform1i(locPalette_, palette == Palette::Sequential ? 0 : 1);
    glUniform1i(locField_, 0);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texture_);
    glBindVertexArray(vao_);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    glBindVertexArray(0);
    return r;
}
