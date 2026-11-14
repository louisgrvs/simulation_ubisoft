#include "render/ParticleRenderer.h"

#include <glad/glad.h>

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

std::string readFile(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("ParticleRenderer : impossible d'ouvrir " + path);
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

constexpr int kFloatsPerParticle = 3;   // x, y, valeur

} // namespace

ParticleRenderer::ParticleRenderer(const std::string& shaderDir) {
    const unsigned vs = compileShader(GL_VERTEX_SHADER, readFile(shaderDir + "particles.vert"), "particles.vert");
    const unsigned fs = compileShader(GL_FRAGMENT_SHADER, readFile(shaderDir + "particles.frag"), "particles.frag");
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
        throw std::runtime_error(std::string("Édition de liens du programme de particules :\n") + log);
    }
    locDomain_ = glGetUniformLocation(program_, "uDomain");
    locSize_   = glGetUniformLocation(program_, "uPointSize");
    locVmin_   = glGetUniformLocation(program_, "uVmin");
    locVmax_   = glGetUniformLocation(program_, "uVmax");
    locMode_   = glGetUniformLocation(program_, "uColorMode");
    locColor_  = glGetUniformLocation(program_, "uColor");

    // VAO : décrit une fois pour toutes la disposition du tampon entrelacé.
    glGenVertexArrays(1, &vao_);
    glGenBuffers(1, &vbo_);
    glBindVertexArray(vao_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    const GLsizei stride = kFloatsPerParticle * sizeof(float);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(2 * sizeof(float)));
    glBindVertexArray(0);
}

ParticleRenderer::~ParticleRenderer() {
    glDeleteBuffers(1, &vbo_);
    glDeleteVertexArrays(1, &vao_);
    glDeleteProgram(program_);
}

void ParticleRenderer::upload(const ParticleSet2D& p, const std::vector<float>& values) {
    count_ = p.size();
    staging_.resize(count_ * kFloatsPerParticle);
    for (std::size_t k = 0; k < count_; ++k) {
        staging_[3 * k + 0] = float(p.x[k].x);
        staging_[3 * k + 1] = float(p.x[k].y);
        staging_[3 * k + 2] = values.empty() ? 0.0f : values[k];
    }
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    if (count_ > capacity_) {
        // Réallocation (rare) : on prévoit 50 % de marge pour les fois suivantes.
        capacity_ = count_ + count_ / 2;
        glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(capacity_ * kFloatsPerParticle * sizeof(float)), nullptr,
                     GL_DYNAMIC_DRAW);
    }
    if (count_ > 0)
        glBufferSubData(GL_ARRAY_BUFFER, 0, GLsizeiptr(staging_.size() * sizeof(float)), staging_.data());
}

void ParticleRenderer::draw(Rect r, Vec2 domain, float pointSizePx, float vmin, float vmax, bool uniform,
                            const float rgb[3]) {
    if (count_ == 0 || r.w <= 0 || r.h <= 0) return;
    glViewport(r.x, r.y, r.w, r.h);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glUseProgram(program_);
    glUniform2f(locDomain_, float(domain.x), float(domain.y));
    glUniform1f(locSize_, pointSizePx);
    glUniform1f(locVmin_, vmin);
    glUniform1f(locVmax_, vmax);
    glUniform1i(locMode_, uniform ? 1 : 0);
    const float white[3] = {0.85f, 0.9f, 1.0f};
    const float* c = rgb ? rgb : white;
    glUniform3f(locColor_, c[0], c[1], c[2]);
    glBindVertexArray(vao_);
    glDrawArrays(GL_POINTS, 0, GLsizei(count_));
    glBindVertexArray(0);
}
