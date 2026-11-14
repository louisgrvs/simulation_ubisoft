#version 430 core
// Chaque point est un carré de uPointSize pixels ; on n'en garde que le disque.
in float vValue;
out vec4 FragColor;

uniform float uVmin;
uniform float uVmax;
uniform int   uColorMode;   // 0 : palette sur vValue, 1 : couleur unie
uniform vec3  uColor;

// Même palette « inferno » que field.frag (approximation de Matt Zucker).
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

void main()
{
    // gl_PointCoord va de (0,0) à (1,1) sur le carré du point.
    vec2 d = 2.0 * gl_PointCoord - 1.0;
    float r2 = dot(d, d);
    if (r2 > 1.0) discard;                       // hors du disque

    if (isnan(vValue) || isinf(vValue)) { FragColor = vec4(1.0, 0.0, 1.0, 1.0); return; }

    vec3 color;
    if (uColorMode == 0) {
        float t = clamp((vValue - uVmin) / max(uVmax - uVmin, 1e-30), 0.0, 1.0);
        // On évite le noir pur du bas de la palette : invisible sur fond sombre.
        color = inferno(0.15 + 0.85 * t);
    } else {
        color = uColor;
    }
    color *= mix(1.0, 0.65, smoothstep(0.55, 1.0, r2));   // liseré plus sombre : on distingue les voisines
    FragColor = vec4(clamp(color, 0.0, 1.0), 1.0);
}
