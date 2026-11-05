#include <glad/glad.h>   // toujours AVANT GLFW
#include <GLFW/glfw3.h>
#include <iostream>
//la bibliotheque des textures : 
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <array>
//la bibliotheque math
#define _USE_MATH_DEFINES
#include <cmath>
#include <corecrt_math_defines.h>
#include <algorithm>

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
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec2 aTexCoord;

out vec3 ourColor;
out vec2 TexCoord;

void main()
{
    gl_Position = vec4(aPos, 1.0);
    ourColor = aColor;
    TexCoord = aTexCoord;
}
)";

//Code source pour le fragment shader (les couleurs), on peut utiliser une variable de type 'in' ou 'uniform' qui sert à passer CPU -> GPU
const char* fragmentShaderSource = R"(
#version 430 core
out vec4 FragColor;

in vec3 ourColor;
in vec2 TexCoord;

uniform sampler2D ourTexture;
uniform float vMin ;
uniform float vMax ;

vec3 palette(float t)   // t dans [0, 1]
{
    vec3 bleu  = vec3(0.23, 0.30, 0.75);
    vec3 blanc = vec3(1.0);
    vec3 rouge = vec3(0.71, 0.02, 0.15);
    return t < 0.5 ? mix(bleu, blanc, 2.0 * t)
                   : mix(blanc, rouge, 2.0 * t - 1.0);
}

void main()
{
    float t = (texture(ourTexture, TexCoord).r - vMin)/ (vMax-vMin) ;
    t = clamp(t, 0.0, 1.0);
    FragColor = vec4(palette(t), 1.0); 
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
    glfwSwapInterval(1);   // synchronisation verticale : 60 images/s maximum, sans ca j'ai des gel de fou :(
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

    //la on les delete pcq on en a plus besoin 
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    // Definition d'un triangle (c'est en 3D mais on a une projection sur l'hyperplan z=0)
    // et la c'est un triangle mais avec des couleurs en plus sur chaque coin
    //float vertices[] = {
        // positions         // colors
       //  0.5f, -0.5f, 0.0f,  1.0f, 0.0f, 0.0f,   // bottom right
       // -0.5f, -0.5f, 0.0f,  0.0f, 1.0f, 0.0f,   // bottom left
       //  0.0f,  0.5f, 0.0f,  0.0f, 0.0f, 1.0f    // top 
    //};

    //Pour un rectangle : 

    float vertices[] = {
        // positions          // colors           // texture coords
         0.5f,  0.5f, 0.0f,   1.0f, 0.0f, 0.0f,   1.0f, 1.0f,   // top right
         0.5f, -0.5f, 0.0f,   0.0f, 1.0f, 0.0f,   1.0f, 0.0f,   // bottom right
        -0.5f, -0.5f, 0.0f,   0.0f, 0.0f, 1.0f,   0.0f, 0.0f,   // bottom left
        -0.5f,  0.5f, 0.0f,   1.0f, 1.0f, 0.0f,   0.0f, 1.0f    // top left 
    };

    // tableau de 64 element qui nous servira pour plus tard 
    std::array<float, 64 * 64> tab{};

    for (int i = 0; i < 64; i++) {
        for (int j = 0; j < 64; j++) {

            float x = (float)i / 64.0f;
            float y = (float)j / 64.0f;
             
            tab[j * 64+ i] = std::sin(2.0f * M_PI * x)
                * std::sin(2.0f * M_PI * y);
        }
    }

    float texCoords[] = {
    0.0f, 0.0f,  // lower-left corner  
    1.0f, 0.0f,  // lower-right corner
    0.5f, 1.0f   // top-center corner
    };

    unsigned int texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    //La on precise la bordure du rectangle et a comment remplir la texture en dehors du triangle
    //float borderColor[] = { 1.0f, 1.0f, 0.0f, 1.0f };
    //glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

    // set the texture wrapping/filtering options (on the currently bound texture object)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    //la texture qu'on choisit la : 
    stbi_set_flip_vertically_on_load(true); //pour remettre a l'endroit
    //int width, height, nrChannels;
    //unsigned char* data = stbi_load("container.jpg", &width, &height, &nrChannels, 0);
    //if (data)
    //{
        //genere la texture en utilisant l'image :
      //  glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, width, height, 0, GL_R32F, GL_UNSIGNED_BYTE, vertices);
       // glGenerateMipmap(GL_TEXTURE_2D);
       // std::cout << "Texture chargée : " << width << " x " << height
        //    << ", " << nrChannels << " canaux" << std::endl;
    //}
    //else
    //{
      //  std::cout << "Failed to load texture" << std::endl;
    //}

    //stbi_image_free(data); //enleve la memoire de l'image


    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, 64, 64, 0, GL_RED, GL_FLOAT, tab.data());

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

    unsigned int indices[] = {
    0, 1, 3,   // premier triangle : haut droit, bas droit, haut gauche
    1, 2, 3    // second triangle  : bas droit, bas gauche, haut gauche
    };

    unsigned int EBO;
    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    //La on dit a openGL comment il devrait interpreter les donnees (les pointeurs quoi)
    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);


    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

    auto mm = std::minmax_element(tab.begin(), tab.end());
    float vMin = *mm.first;    // plus petite valeur de tab
    float vMax = *mm.second;   // plus grande valeur de tab

    // Boucle de rendu
    while (!glfwWindowShouldClose(window))
    {
        // input
        processInput(window);

        // render
        // clear the colorbuffer
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // be sure to activate the shader
        glUseProgram(shaderProgram);

        // update the uniform color
        float timeValue = glfwGetTime();
        float greenValue = sin(timeValue) / 2.0f + 0.5f;

        // 1. Trouver les adresses des uniforms (une seule fois)
        int locVMin = glGetUniformLocation(shaderProgram, "vMin");
        int locVMax = glGetUniformLocation(shaderProgram, "vMax");

        // 2. Les remplir (le programme doit être actif)
        glUseProgram(shaderProgram);
        glUniform1f(locVMin, vMin);   // minimum de sin·sin
        glUniform1f(locVMax, vMax);   // maximum de sin·sin

        //attache la structure au VAO
        glBindTexture(GL_TEXTURE_2D, texture);
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);

        // swap buffers and poll IO events
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}
