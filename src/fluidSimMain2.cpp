// GLEW
#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <unistd.h>
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

#define TIME 2000000
#define SIZE 15
#define START 1
#define END SIZE + 1
#define at(i,j) ((i) + (SIZE + 2)*(j))
#define SWAP(x0,x) {float *tmp=x0;x0=x;x=tmp;}

#define NUMCELLS (SIZE) * (SIZE)
#define SIZEOFCELL sizeof(GLfloat) * 4 * 6

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);

void printVec(glm::vec3 toPrint, std::string vecName);
void createVertices(GLfloat *ptr);
int getIndex(int i, int j);
void setDensity(float *density);


void dens_step(int N, float * x, float * x0, float * u, float * v, float diff, float dt);
void vel_step(int N, float * u, float * v, float * u0, float * v0, float visc, float dt);

void drawDensity(float *density);
void addDensity(float *density);
void addForce(float *velocityX, float *velocityY, float forceX, float forceY);
void addForceAwayFromCenter(float *velocityX, float *velocityY);
void setCellDensity(glm::vec3 cell, float densityVal);

const GLuint WIDTH = 800, HEIGHT = 600;
bool keys[1024]; 
bool firstMouse = true;
GLfloat lastX = WIDTH / 2.0;
GLfloat lastY = HEIGHT / 2.0;

int numRects = 0;
float *density;
float *densityPrev;
float *velocityX;
float *velocityY;
float *velocityPrevX;
float *velocityPrevY;


int sizeOfData = (int) (24 * SIZE * SIZE * sizeof(GLfloat));
GLfloat *vertices = (GLfloat *) malloc(sizeOfData);
GLuint VBO, VAO;

// Shaders
const GLchar* vertexShaderSource = "#version 330 core\n"
    "layout (location = 0) in vec3 position;\n"
    "layout (location = 1) in float density;\n"

    "out float densityColor;\n"
    
    "void main()\n"
    "{\n"
        "gl_Position = vec4(position.x, position.y, position.z, 1.0);\n"
        "densityColor = density;\n"
    "}\0";
const GLchar* fragmentShaderSource = "#version 330 core\n"
    "in float densityColor;\n"
    "out vec4 color;\n"
    "void main()\n"
    "{\n"
        "vec3 density = vec3(densityColor);\n"
        "color = vec4(density, 1.0f);\n"
    "}\n\0";

// The MAIN function, from here we start the application and run the game loop
int main()
{
        density = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
    densityPrev = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

    velocityX = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
    velocityY = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

    velocityPrevX = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
    velocityPrevY = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

    glfwInit();

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "LearnOpenGL", nullptr, nullptr);    
    if (window == nullptr)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    glewExperimental = GL_TRUE;

    if (glewInit() != GLEW_OK)
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
        return -1;
    }    

    // Define the viewport dimensions
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);  
    glViewport(0, 0, width, height);
    // glEnable(GL_DEPTH_TEST);

    // Build and compile our shader program
    // Vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    // Check for compile time errors
    GLint success;
    GLchar infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    }
    // Fragment shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    // Check for compile time errors
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
    }
    // Link shaders
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    // Check for linking errors
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    }
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);


    //each vertex has 4 points (x,y,z) density
    //each cell has 6 vertexes

    //each cell needs to be size: sizeof(GLfloat) * 4 * 6
    createVertices(vertices);

 
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    // Bind the Vertex Array Object first, then bind and set vertex buffer(s) and attribute pointer(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeOfData, vertices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 4 * sizeof(GLfloat), (GLvoid*)(3* sizeof(GLfloat)));
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0); // Note that this is allowed, the call to glVertexAttribPointer registered VBO as the currently bound vertex buffer object so afterwards we can safely unbind

    glBindVertexArray(0); // Unbind VAO (it's always a good thing to unbind any buffer/array to prevent strange bugs), remember: do NOT unbind the EBO, keep it bound to this VAO


    // Uncommenting this call will result in wireframe polygons.
    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    float visc = 100.0f;
    float diff = 1.0f;
    float dt = 0.0f;

    // addDensity(density);
    // density[at(0, 0)] = 1.0f;
    bool toggle = false;
    for (int i = 1; i < SIZE-1; i++) {
        for (int j = 1; j < SIZE-1; j++) {
            if (toggle) {
                densityPrev[at(i,j)] = 1.0f;
            }
            else {
                densityPrev[at(i,j)] = 0.0f;
            }
            density[at(i,j)] = 0.0f;

            velocityPrevX[at(i,j)] = 0.0f;
            velocityPrevY[at(i,j)] = 0.0f;
            
            toggle = !toggle;
        }        
    }
    // setDensity(densityPrev);

    // Game loop
    while (!glfwWindowShouldClose(window))
    {
        // Check if any events have been activiated (key pressed, mouse moved etc.) and call corresponding response functions
        glfwPollEvents();

        // Render
        // Clear the colorbuffer
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        vel_step(SIZE, velocityX, velocityY, velocityPrevX, velocityPrevY, visc, dt);
        dens_step(SIZE, density, densityPrev, velocityX, velocityY, diff, dt);
        setDensity(density);
        dt++;

        // Draw our first triangle
        glUseProgram(shaderProgram);
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, numRects * 3 * 2);
        // glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        glBindVertexArray(0);

        // Swap the screen buffers
        glfwSwapBuffers(window);
    }
    // Properly de-allocate all resources once they've outlived their purpose
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    // Terminate GLFW, clearing any resources allocated by GLFW.
    glfwTerminate();
    return 0;
}

void addDensity(float *density) {
    // for (int i = 100; i < 250; i++) {
    //     for (int j = 100; j < 250; j++) {
    //         density[at(i, j)] = 1.0f;
    //     }   
    // }
    // for (int i = 50; i < 60; i++) {
    //     for (int j = 0; j < 30; j++) {
    //         density[at(2+i, 70 + j)] = 1.0f;
    //         density[at(3+i, 70 + j)] = 1.0f;
    //         density[at(4+i, 70 + j)] = 1.0f;
    //         density[at(3+i, 70 + j)] = 1.0f;
    //         density[at(3+i, 70 + j)] = 1.0f;        
    //     }
    // }
    density[at(1,1)] = 1.0f;
    // density[at(3,5)] = 1.0f;
    // density[at(4,5)] = 1.0f;
    // density[at(5,3)] = 1.0f;
    // density[at(5,4)] = 1.0f;
    // density[at(5,5)] = 1.0f;
}

void createVertices(GLfloat *ptr) {
    bool useDensity = 0.0f;
    float cellSize = 2.0f / SIZE;

    std::cout << cellSize << std::endl;

    int counter = 0;
    float k = 0.0f;
    float x = 0.0f;
    float y = 1.0f;
    cellSize = 2.0f / SIZE;
    for (int j = 0; j < SIZE; j++) {
        x = -1.0f;
        // std::cout << std::endl;
        for (int i = 0; i < SIZE; i++) {
            // std::cout << y << " ";
                ptr[counter++] = x;
                ptr[counter++] = y;
                ptr[counter++] = k;
                ptr[counter++] = 0.0f;

                // lowerLeftVert
                ptr[counter++] = x;
                ptr[counter++] = y - cellSize;
                ptr[counter++] = k;
                ptr[counter++] = 0.0f;

                // upperRightVert
                ptr[counter++] = x + cellSize;
                ptr[counter++] = y;
                ptr[counter++] = k;
                ptr[counter++] = 0.0f;

            //right triangle
                // lowerLeftVert
                ptr[counter++] = x;
                ptr[counter++] = y - cellSize;
                ptr[counter++] = k;
                ptr[counter++] = 0.0f;

                // upperRightVert
                ptr[counter++] = x + cellSize;
                ptr[counter++] = y;
                ptr[counter++] = k;
                ptr[counter++] = 0.0f;

                // lowerRightVert
                ptr[counter++] = x + cellSize;
                ptr[counter++] = y - cellSize;
                ptr[counter++] = k;
                ptr[counter++] = 0.0f;
                numRects++;


            x += cellSize;
        }
        y -= cellSize;
    }
}

/* returns offset into ptr of the upper left hand corner
    ex: 0,0 -> 0
    ex: 1,0 -> SIZEOFCELL;
*/
int getIndex(int i, int j) {
    return ((i) + (SIZE)*j) * 24;
}

glm::vec3 getCell(glm::vec3 mousePos) {
    glm::vec3 cell;
    float screenCellSizeX = WIDTH / SIZE;
    float screenCellSizeY = HEIGHT / SIZE;

    int xCoord = mousePos.x / screenCellSizeX;
    int yCoord = mousePos.y / screenCellSizeY;
    
    cell.x = xCoord;
    cell.y = yCoord;
    
    return cell;
}

// Is called whenever a key is pressed/released via GLFW
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
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

} 

void setCellDensity(glm::vec3 cell, float densityVal) {
    int index1 = getIndex(cell.x, cell.y) + 3;
    int index2 = index1 + 4;
    int index3 = index1 + 8;
    int index4 = index1 + 12;
    int index5 = index1 + 16;
    int index6 = index1 + 20;
    vertices[index1] = densityVal;
    vertices[index2] = densityVal;
    vertices[index3] = densityVal;
    vertices[index4] = densityVal;
    vertices[index5] = densityVal;
    vertices[index6] = densityVal;
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferSubData(GL_ARRAY_BUFFER, index1*sizeof(float), sizeof(float), &vertices[index1]);
    glBufferSubData(GL_ARRAY_BUFFER, index2*sizeof(float), sizeof(float), &vertices[index2]);
    glBufferSubData(GL_ARRAY_BUFFER, index3*sizeof(float), sizeof(float), &vertices[index3]);
    glBufferSubData(GL_ARRAY_BUFFER, index4*sizeof(float), sizeof(float), &vertices[index4]);
    glBufferSubData(GL_ARRAY_BUFFER, index5*sizeof(float), sizeof(float), &vertices[index5]);
    glBufferSubData(GL_ARRAY_BUFFER, index6*sizeof(float), sizeof(float), &vertices[index6]);

    printVec(cell, "cell");
    std::cout << "   " << densityVal << std::endl;
}

void setDensity(float *density) {
    for (int i = 1; i < SIZE-1; i++) {
        for (int j = 1; j < SIZE-1; j++) {
            setCellDensity(glm::vec3(i, j, 0.0), density[at(i,j)]);
        }        
    }
}

void addDensity(glm::vec3 cell) {
    densityPrev[at((int)cell.x, (int)cell.y)] = 1.0f;
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    if (state == GLFW_PRESS) {
        std::cout << "x: " << lastX << " y: " << lastY << std::endl;
        glm::vec3 cellClicked = getCell(glm::vec3(lastX, lastY, 0.0f));
        setCellDensity(getCell(glm::vec3(lastX, lastY, 0.0f)), 1.0f);
        addDensity(getCell(glm::vec3(lastX, lastY, 0.0f)));
        printVec(cellClicked, "cellClicked");
    }    
}

void printVec(glm::vec3 toPrint, std::string vecName) {
    std::cout << vecName << ": ";
    std::cout << toPrint.x << " ";
    std::cout << toPrint.y << " ";
    std::cout << toPrint.z << " ";
    std::cout << std::endl;
}


// Fluid solver code

void add_source(int N, float *x, float *s, float dt) {
    int size = (N + 2) * (N + 2);
    for (int i = 0; i < size; i++) {
        x[i] += dt * s[i];
    }
}

void set_bnd(int N, int b, float * x ) {
    int i;
    for ( i=1 ; i<=N ; i++ ) {
        x[at(0, i)]      = b==1 ? -x[at(1, i)] : x[at(1,i)];
        x[at(N+1, i)]    = b==1 ? -x[at(N, i)] : x[at(N,i)];
        x[at(i, 0)]      = b==2 ? -x[at(i, 1)] : x[at(i, 1)];
        x[at(i, N+1)]    = b==2 ? -x[at(i, N)] : x[at(i,N)];
    }
    x[at(0 ,0 )] = 0.5*(x[at(1,0 )]+x[at(0 ,1)]);
    x[at(0 ,N+1)] = 0.5*(x[at(1,N+1)]+x[at(0 ,N )]);
    x[at(N+1,0 )] = 0.5*(x[at(N,0 )]+x[at(N+1,1)]);
    x[at(N+1,N+1)] = 0.5*(x[at(N,N+1)]+x[at(N+1,N)]);
}

void diffuse(int N, int b, float *x, float *x0, float diff, float dt) {
    int i, j, k;
    float a = dt * diff * N * N;

    for (k = 0; k < 20; k++) {
        for(int i = 1; i <= N; i++) {
            for (int j = 0; j <= N; j++) {
                float surrDensity = x[at(i-1,j)] + x[at(i+1,j)] + x[at(i,j-1)] + x[at(i,j+1)];
                x[at(i,j)] = (x0[at(i,j)] + a * surrDensity) / (1 + 4*a);
            }
        }
        set_bnd(N,b,x);
    }
}

void advect(int N, int b, float *d, float *d0, float *u, float *v, float dt) {
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;

    dt0 = dt * N;
    for(int i = 1; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            x = i-dt0*u[at(i,j)]; y = j-dt0*v[at(i,j)];
            if (x<0.5) x=0.5; if (x>N+0.5) x=N+0.5; i0=(int)x; i1=i0+1;
            if (y<0.5) y=0.5; if (y>N+0.5) y=N+0.5; j0=(int)y; j1=j0+1; 
            s1 = x-i0; s0 = 1-s1; 
            t1 = y-j0; t0 = 1-t1;
            d[at(i,j)] = s0*(t0*d0[at(i0,j0)]+t1*d0[at(i0,j1)])+s1*(t0*d0[at(i1,j0)]+t1*d0[at(i1,j1)]);
        }
    }
    set_bnd(N, b, d);
}

void project ( int N, float * u, float * v, float * p, float * div ) {
    int i, j, k;
    float h;

    h = 1.0 / N;
    for(int i = 1; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            div[at(i,j)] = -0.5*h*(u[at(i+1,j)]-u[at(i-1,j)]+v[at(i,j+1)]-v[at(i,j-1)]);
            p[at(i,j)] = 0;
        }
    }

    set_bnd(N, 0, div);
    set_bnd(N, 0, p);

    for (k = 0; k < 20; k++) {
        for(int i = 1; i <= N; i++) {
            for (int j = 0; j <= N; j++) {
                p[at(i,j)] = (div[at(i,j)]+p[at(i-1,j)]+p[at(i+1,j)]+p[at(i,j-1)]+p[at(i,j+1)])/4;
            }
        }
        set_bnd(N, 0, p);
    }

    for(int i = 1; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
        u[at(i,j)] -= 0.5*(p[at(i+1,j)]-p[at(i-1,j)])/h;
        v[at(i,j)] -= 0.5*(p[at(i,j+1)]-p[at(i,j-1)])/h;
    }
}
    set_bnd(N, 1, u); 
    set_bnd(N, 2, v);
}

void dens_step(int N, float * x, float * x0, float * u, float * v, float diff, float dt) {
    add_source(N, x, x0, dt);
    SWAP(x0, x); 
    diffuse(N, 0, x, x0, diff, dt);
    SWAP(x0, x);
    advect(N, 0, x, x0, u, v, dt);
}

void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt ) {
    add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
    SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt );
    SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt );
    project ( N, u, v, u0, v0 );
    SWAP ( u0, u ); SWAP ( v0, v );
    advect ( N, 1, u, u0, u0, v0, dt ); advect ( N, 2, v, v0, u0, v0, dt ); project ( N, u, v, u0, v0 );
}