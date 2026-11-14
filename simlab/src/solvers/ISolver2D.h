#pragma once
// Contrat commun à tous les solveurs (Stable Fluids, FLIP, MPM).
// L'application ne connaît que cette interface pour piloter, afficher et mesurer.
//
// J5 : version alignée sur le Field2D des jours 3-4 (réels double, décalage
// (ox, oy)), qui remplace la version modèle Field2D<Real> du jour 2.
#include <string>
#include <vector>
#include "grid/Field2D.h"

// Vue en lecture seule sur un champ à afficher.
struct FieldView {
    std::string name;
    const Field2D* field = nullptr;
    bool signedField = false;   // vrai : palette divergente centrée en 0
};

// Une grandeur de validation affichée dans le panneau et exportable en CSV.
struct Diagnostic {
    std::string name;
    double value = 0.0;
};

class ISolver2D {
public:
    virtual ~ISolver2D() = default;

    virtual std::string name() const = 0;
    virtual void reset() = 0;
    virtual void step(real dt) = 0;            // avance de dt secondes
    virtual real stableDt() const = 0;         // pas conseillé (CFL cible, précision)
    virtual real time() const = 0;             // temps simulé, en secondes
    virtual std::vector<Diagnostic> diagnostics() const = 0;
    virtual std::vector<FieldView> fields() const = 0;
};
