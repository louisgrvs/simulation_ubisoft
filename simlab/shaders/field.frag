#version 430 core
// Palette appliquée à un champ scalaire stocké dans une texture GL_R32F.
in vec2 vUV;
out vec4 FragColor;

uniform sampler2D uField;
uniform float uVmin;
uniform float uVmax;
uniform int   uPalette;    // 0 : séquentielle (fumée, norme), 1 : divergente (signé)

// Approximation polynomiale de la palette « inferno » (matplotlib),
// d'après Matt Zucker (domaine public). Noir -> violet -> orange -> jaune pâle.
vec3 inferno(float t)
{
    const vec3 c0 = vec3( 0.0002189403691192265,  0.001651004631001012, -0.01948089843709184);
    const vec3 c1 = vec3( 0.1065134194856116,     0.5639564367884091,    3.932712388889277);
    const vec3 c2 = vec3( 11.60249308247187,     -3.972853965665698,   -15.9423941062914);
    const vec3 c3 = vec3(-41.70399613139459,     17.43639888205313,     44.35414519872813);
    const vec3 c4 = vec3( 77.162935699427,      -33.40235894210092,    -81.80730925738993);
    const vec3 c5 = vec3(-71.31942824499214,     32.62606426397723,     73.20951985803202);
    const vec3 c6 = vec3( 25.13112622477341,    -12.24266895238567,    -23.07032500287172);
    return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

// Palette divergente pour fond sombre : bleu (négatif) <- noir (0) -> orange (positif).
vec3 diverging(float s)
{
    const vec3 neg = vec3(0.30, 0.62, 1.00);
    const vec3 pos = vec3(1.00, 0.55, 0.20);
    float a = pow(abs(s), 0.75);       // un peu plus de contraste près de 0
    return s < 0.0 ? a * neg : a * pos;
}

void main()
{
    float v = texture(uField, vUV).r;

    // NaN ou infini : magenta, impossible à rater quand une simulation explose.
    if (isnan(v) || isinf(v)) { FragColor = vec4(1.0, 0.0, 1.0, 1.0); return; }

    float range = max(uVmax - uVmin, 1e-30);
    float t = clamp((v - uVmin) / range, 0.0, 1.0);
    vec3 color = (uPalette == 0) ? inferno(t) : diverging(2.0 * t - 1.0);
    FragColor = vec4(clamp(color, 0.0, 1.0), 1.0);
}
