#include <glad/glad.h>   // toujours AVANT GLFW
#include <GLFW/glfw3.h>
#include <iostream>


// Fonctions définies HORS de main, je crois que c'est pour actualiser la taille d'une fenêtre
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}
// ça ça nous dis le bouton échap elle est appuyer ou non
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

//Code source pour le vertex shader (demande par le tuto), sert a creer un shader le plus simple possible
// Vertex shader : appelé une fois par sommet
const char* vertexShaderSource = R"(
#version 430 core
layout (location = 0) in vec3 aPos;
void main()
{
    gl_Position = vec4(aPos, 1.0);
}
)";

//Code source pour le fragment shader (les couleurs)
const char* fragmentShaderSource = R"(
#version 430 core
out vec4 FragColor;
void main()
{
    FragColor = vec4(0.5, 0.0, 0.0, 1.0);  // orange
}
)";

int main()
{
    // Initialiser GLFW et demander un contexte OpenGL 4.3 Core
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);


    // Creer la fenetre
    GLFWwindow* window = glfwCreateWindow(800, 600, "Simulation fluides", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Échec de la création de la fenêtre GLFW" << std::endl;
        glfwTerminate();
        return -1;
    }

    // Rendre le contexte courant
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Charger GLAD (uniquement maintenant que le contexte existe)
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Échec de l'initialisation de GLAD" << std::endl;
        return -1;
    }

    // Appels OpenGL c'est bon !
    std::cout << "OpenGL " << glGetString(GL_VERSION) << std::endl;
    glViewport(0, 0, 800, 600);

    //Et ca ca cree un shader : 
    unsigned int vertexShader;
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    //utilise le code source pour le shader
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);

    // Un ptit code pour expliquer si le shader a bien ete creer ou pas (etape precedente)
    int  success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);

    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    }

    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success) {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cout << "Erreur fragment shader :\n" << infoLog << std::endl;
    }

    //La on cree un shader program : 
    unsigned int shaderProgram;
    shaderProgram = glCreateProgram();

    //la on attache nos trucs compile a notre shader
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);

    //est-ce que le shader prgram est bien attache ? 
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::LINK_FAILED\n" << infoLog << std::endl;
    }

    glUseProgram(shaderProgram);

    //la on les delete pcq on en a plus besoin 
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    //La on dit a openGL comment il devrait interpreter les donnees
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);


    // Definition d'un triangle (c'est en 3D mais on a une projection sur l'hyperplan z=0)
    float vertices[] = {
    -0.5f, -0.5f, 0.0f,
     0.5f, -0.5f, 0.0f,
     0.0f,  0.5f, 0.0f
    };

    //Definition et initialisation d'un Vertex Buffer Object (VBO),
    // qui est en gros liste des données sur laquelle on dessine nos triangles
    unsigned int VBO;
    glGenBuffers(1, &VBO);

    //maintenant on creer notre VAO : qui dit comment lire le VBO
    unsigned int VAO;
    glGenVertexArrays(1, &VAO);

    // 1. bind Vertex Array Object
    glBindVertexArray(VAO);
    // 2. copy our vertices array in a buffer for OpenGL to use
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    // 3. then set our vertex attributes pointers
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);


    // Boucle de rendu
    while (!glfwWindowShouldClose(window))
    {
        processInput(window);

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glUseProgram(shaderProgram);
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        //someOpenGLFunctionThatDrawsOurTriangle();

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
