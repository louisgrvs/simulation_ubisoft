#version 430 core
// Un rectangle plein écran sans aucun tampon de sommets : les 4 coins sont
// déduits de gl_VertexID (dessin en GL_TRIANGLE_STRIP de 4 sommets).
// Le rectangle couvre tout le viewport ; c'est glViewport qui le place à l'écran.
out vec2 vUV;

void main()
{
    vec2 corner = vec2(float(gl_VertexID & 1), float((gl_VertexID >> 1) & 1));  // (0,0) (1,0) (0,1) (1,1)
    vUV = corner;                      // v = 0 en bas : ligne j = 0 du champ (CONVENTIONS.md)
    gl_Position = vec4(2.0 * corner - 1.0, 0.0, 1.0);
}
