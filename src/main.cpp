// GLEW
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <iomanip> 
#include <map>
#include <string>
#include <iterator>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Classes/Shader.hpp"
#include "Classes/Texture.hpp"
#include "Classes/Camera.hpp"
#include "Classes/Rainbow.hpp"

#include "MarchingCubes/MarchingCubes.hpp"

#include "constants.hpp"

void key_callback(GLFWwindow *window, int key, int scancode, int action, int mode);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void do_movement();
void setGravity();

void setBunny();
int generateAttachmentTexture(GLboolean depth, GLboolean stencil);
void createFBO(GLuint& fb, GLuint& tex);

const GLuint WIDTH = 1200, HEIGHT = 1000;

bool keys[1024]; 
bool firstMouse = true;
GLfloat lastX = WIDTH / 2.0;
GLfloat lastY = HEIGHT / 2.0;

GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;

glm::vec3 lightPos(1.2f, 1.0f, 2.0f);

std::vector<float> data;
int numVertices = 0;
int strideSize = 3;

std::vector<glm::vec3> vertices;
std::vector<glm::vec3> normals;

Camera camera(glm::vec3(0.0f));
MarchingCubes marchingCubes;

bool doForce = false;
bool useNormal = false;

GLFWwindow* window;

int main()
{
    glfwInit();

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    window = glfwCreateWindow(WIDTH, HEIGHT, "LearnOpenGL", nullptr, nullptr);    
    if (window == nullptr)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glViewport(0, 0, WIDTH, HEIGHT);
    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, mouse_callback);

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    
    glewExperimental = GL_TRUE;

    if (glewInit() != GLEW_OK)
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return -1;
    }    

    camera.Position = glm::vec3(0, 10, 10);

    // Bunny:
    setBunny();


    Shader *shader = new Shader(SHADER_PATH "lamp.vert", SHADER_PATH "lamp.frag");
    Shader *screenShader = new Shader(SHADER_PATH "screen.vert", SHADER_PATH "screen.frag");
    // set up quad
    // Vertex attributes for a quad that fills the entire screen in Normalized Device Coordinates.
    GLfloat quadVertices[] = {
        // Positions   // TexCoords
        -1.0f,  1.0f,  0.0f, 1.0f,
        -1.0f, -1.0f,  0.0f, 0.0f,
         1.0f, -1.0f,  1.0f, 0.0f,

        -1.0f,  1.0f,  0.0f, 1.0f,
         1.0f, -1.0f,  1.0f, 0.0f,
         1.0f,  1.0f,  1.0f, 1.0f
    };  

    // Setup screen VAO
    GLuint quadVAO, quadVBO;
    glGenVertexArrays(1, &quadVAO);
    glGenBuffers(1, &quadVBO);
    glBindVertexArray(quadVAO);
    glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (GLvoid*)(2 * sizeof(GLfloat)));
    glBindVertexArray(0);



    GLuint frameBuf[2];
    GLuint texBuf[2];

    glGenFramebuffers(2, frameBuf);
    glGenTextures(2, texBuf);
    createFBO(frameBuf[0], texBuf[0]);
        /* not sure what this does... */
        GLenum DrawBuffers[1] = {GL_COLOR_ATTACHMENT0};
        glDrawBuffers(1, DrawBuffers);
    createFBO(frameBuf[1], texBuf[1]);

    bool toggle = true;

    double lastTime = glfwGetTime();
    int nbFrames = 0;
    while (!glfwWindowShouldClose(window))
    {
        // Calculate deltatime of current frame
        GLfloat currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        double currentTime = glfwGetTime();
        nbFrames++;
        if ( currentTime - lastTime >= 1.0 ){ // If last prinf() was more than 1 sec ago
            // printf and reset timer
            // printf("%f ms/frame\n", 1000.0/double(nbFrames));
            nbFrames = 0;
            lastTime += 1.0;
        }
        // Check if any events have been activiated (key pressed, mouse moved etc.) and call corresponding response functions
        glfwPollEvents();
        do_movement();



        /////////////////////////////////////////////////////
        // Bind to framebuffer and draw to color texture 
        // as we normally would.
        ////////////////////////////////////////////////////
        glBindFramebuffer(GL_FRAMEBUFFER, frameBuf[toggle]);
        // Clear all attached buffers        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // We're not using stencil buffer so why bother with clearing?

        shader->use();
        glBindTexture(GL_TEXTURE_2D, texBuf[!toggle]);

        // render scene to frame buffer




            /////////////////////////////////////////////////////
            // Bind to default framebuffer again and draw the 
            // quad plane with attched screen texture.
            // //////////////////////////////////////////////////
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
            // Clear all relevant buffers
            glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set clear color to white (not really necessery actually, since we won't be able to see behind the quad anyways)
            glClear(GL_COLOR_BUFFER_BIT);
            glDisable(GL_DEPTH_TEST); // We don't care about depth information when rendering a single quad

            // Draw Screen
            screenShader->use();
            glBindVertexArray(quadVAO);
            glBindTexture(GL_TEXTURE_2D, texBuf[toggle]);   // Use the color attachment texture as the texture of the quad plane
            glDrawArrays(GL_TRIANGLES, 0, 6);
            glBindVertexArray(0);

            toggle = !toggle;
            // Swap the screen buffers
            glfwSwapBuffers(window);
    }


    glDeleteFramebuffers(1, &frameBuf[0]);
    glDeleteFramebuffers(1, &frameBuf[1]);
    return 0;
}

/*
Helper function to create the framebuffer object and associated texture to write to
*/
void createFBO(GLuint& fb, GLuint& tex) {
    //initialize FBO (global memory)
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    //set up framebuffer
    glBindFramebuffer(GL_FRAMEBUFFER, fb);
    //set up texture
    glBindTexture(GL_TEXTURE_2D, tex);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex, 0);

    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        std::cout << "Error setting up frame buffer - exiting" << std::endl;
        exit(0);
    }
}

int generateAttachmentTexture(GLboolean depth, GLboolean stencil)
{
    // What enum to use?
    GLenum attachment_type;
    if(!depth && !stencil)
        attachment_type = GL_RGB;
    else if(depth && !stencil)
        attachment_type = GL_DEPTH_COMPONENT;
    else if(!depth && stencil)
        attachment_type = GL_STENCIL_INDEX;

    //Generate texture ID and load texture data 
    GLuint textureID;
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);
    if(!depth && !stencil)
        glTexImage2D(GL_TEXTURE_2D, 0, attachment_type, WIDTH, HEIGHT, 0, attachment_type, GL_UNSIGNED_BYTE, NULL);
    else // Using both a stencil and depth test, needs special format arguments
        glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, WIDTH, HEIGHT, 0, GL_DEPTH_STENCIL, GL_UNSIGNED_INT_24_8, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    return textureID;
}


void setBunny() {
    const char *path = "../resources/bunny.obj";

    FILE * file = fopen(path, "r");
    if( file == NULL ){
        printf("Impossible to open the file !\n");
        exit(0);
    }
    while( 1 ){
        char lineHeader[128];
        // read the first word of the line
        int res = fscanf(file, "%s", lineHeader);
        if (res == EOF)
            break; // EOF = End Of File. Quit the loop.

        if ( strcmp( lineHeader, "v" ) == 0 ){
            glm::vec3 vertex;
            fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z );
            vertex *=20;
            data.push_back(vertex.x);
            data.push_back(vertex.y);
            data.push_back(vertex.z);

            glm::vec3 color(130, 170, 255);
            data.push_back(color.x);
            data.push_back(color.y);
            data.push_back(color.z);

            if (!useNormal) {
                glm::vec3 force = glm::normalize(glm::vec3(-0.466541, 2.08864, 0.117601) - vertex);
                data.push_back(force.x);
                data.push_back(force.y);
                data.push_back(force.z);
            }
            numVertices++;
        }
        if ( strcmp( lineHeader, "vn" ) == 0 ){
            glm::vec3 normal;
            fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z );

            // glm::vec3 force = glm::normalize(glm::vec3(-0.466541, 2.08864, 0.117601) - vertex);
            if (useNormal) {

                data.push_back(normal.x);
                data.push_back(normal.y);
                data.push_back(normal.z);                
            }
        }
        camera.Position = glm::vec3(0, 2, 50);
        // if ( strcmp( lineHeader, "f" ) == 0 ){
        //     glm::vec3 triangle;
        //     fscanf(file, "%f %f %f\n", &triangle.x, &triangle.y, &triangle.z );

        //     std::cout << std::endl << triangle.x << " ";
        //     std::cout << triangle.y << " ";
        //     std::cout << triangle.z << " ";

        //     glm::vec3 vert1, vert2, vert3;

        //     float offset;
        //     offset = (triangle.x -1) * 3;
        //     vert1.x = data[offset + 0];
        //     vert1.y = data[offset + 1];
        //     vert1.z = data[offset + 2];

        //     offset = (triangle.y -1) * 3;
        //     vert2.x = data[offset + 0];
        //     vert2.y = data[offset + 1];
        //     vert2.z = data[offset + 2];

        //     offset = (triangle.z -1) * 3;
        //     vert3.x = data[offset + 0];
        //     vert3.y = data[offset + 1];
        //     vert3.z = data[offset + 2];


        //     addLine(vert1, vert2);
        //     addLine(vert1, vert3);
        //     addLine(vert2, vert3);

        //     numFaces++;
        // }

    }
    std::cout << "numVertices: " << numVertices << std::endl;
}

void setGravity() {

    // for (glm::vec3 vert : vertices) {
    //     data.push_back(vert.x);
    //     data.push_back(vert.y);
    //     data.push_back(vert.z);

    //     numVertices++;

    //     // float dir = vert.x - marchingCubes.getCenter().x;
    //     glm::vec3 col = rainbow.getColor(vert.y);
    //     data.push_back(col.x);
    //     data.push_back(col.y);
    //     data.push_back(col.z);


    //     glm::vec3 force = glm::normalize(marchingCubes.getCenter()-vert);
    //     data.push_back(randFloat(-1, 1) * force.x);
    //     data.push_back(randFloat(-1, 1) * force.y);
    //     data.push_back(randFloat(-1, 1) * force.z);
    // }

    // for (int j = 0; j < data.size()/3; j++) {
    //     // printf("%f ", data[j]);
    //     // data[6*j] = 0;
    //     data[0*j+1] = -9.8;
    //     // data[6*j+2] =0;
    // }
}

void key_callback(GLFWwindow *window, int key, int scancode, int action, int mode) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }    
    if (key == GLFW_KEY_SLASH && action == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }
    if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
        // setGravity();
        doForce = !doForce;
    }
    if(action == GLFW_PRESS)
        keys[key] = true;
    else if(action == GLFW_RELEASE)
        keys[key] = false; 
}

// Moves/alters the camera positions based on user input
void do_movement() {
    // Camera controls
    if(keys[GLFW_KEY_W])
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if(keys[GLFW_KEY_S])
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if(keys[GLFW_KEY_A])
        camera.ProcessKeyboard(LEFT, deltaTime);
    if(keys[GLFW_KEY_D])
        camera.ProcessKeyboard(RIGHT, deltaTime);
    if(keys[GLFW_KEY_J])
        camera.ProcessKeyboard(SPEEDUP, deltaTime);
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos){
    if(firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }
  
    GLfloat xoffset = xpos - lastX;
    GLfloat yoffset = lastY - ypos; 
    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}  