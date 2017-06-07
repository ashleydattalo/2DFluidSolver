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

#include "FluidSim/Image.h"
#include "FluidSim/types.h"

#include "constants.hpp"

#define TIME 2000000
#define SIZE 200
#define at(i,j) ((i) + (SIZE + 2)*(j))
// #define SWAP(x0,x) {float *tmp=x0;x0=x;x=tmp;}
#define START 1
#define END SIZE + 1

// #define IX(i,j) ((i)+(N+2)*(j))
// #define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
// #define END_FOR }}
#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) {\
for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

void dens_step(int N, float * x, float * x0, float * u, float * v, float diff, float dt);
void vel_step(int N, float * u, float * v, float * u0, float * v0, float visc, float dt);

void printArr(float *arr, std::string arrName);
void drawDensity(float *density);
void addDensity(float *density);
void addForce(float *velocityX, float *velocityY, float forceX, float forceY);
void addForceAwayFromCenter(float *velocityX, float *velocityY);
Image img(SIZE + 2, SIZE + 2);

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);

glm::vec3 getCell(glm::vec3 mousePos);
void printVec(glm::vec3 toPrint, std::string vecName);
void createGrid();
void addCellDensity(glm::vec3 cellClicked);
void addCellVelocity(glm::vec3 cellClicked, glm::vec2 offset);

void getFromUI(float *densityPrev, float *velocityPrevX, float *velocityPrevY);
void setDensity(float *density);

const GLuint WIDTH = 800, HEIGHT = 600;
bool keys[1024]; 
bool firstMouse = true;
GLfloat lastX = WIDTH / 2.0;
GLfloat lastY = HEIGHT / 2.0;
GLfloat currX, currY;
bool go;

GLuint VBO, VAO;
std::vector<glm::vec3> vertices;

struct Triangle {
    glm::vec3 p1;
    glm::vec3 p2;
    glm::vec3 p3;
};

struct Rect {
    int id;
    Triangle leftTri;
    Triangle rightTri;
    bool hasVelocity;
    glm::vec2 velocity;
};

Rect data[(SIZE + 2) * (SIZE + 2)];
int numRects = 0;

bool step = true;
bool pauseSim = false;

int main()
{
    //allocated data
    float *density = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
    float *velocityX = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
    float *velocityY = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

    float *densityPrev = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
    float *velocityPrevX = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
    float *velocityPrevY = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

    memset(densityPrev, 0, sizeof(float) * (SIZE + 2) * (SIZE + 2));
    memset(velocityPrevX, 0, sizeof(float) * (SIZE + 2) * (SIZE + 2));
    memset(velocityPrevY, 0, sizeof(float) * (SIZE + 2) * (SIZE + 2));

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

    Shader *shader = new Shader(SHADER_PATH "waterSim.vert", SHADER_PATH "waterSim.frag");
    
    createGrid();

    // Define the viewport dimensions
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);  
    glViewport(0, 0, width, height);

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), &vertices[0], GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    glBindVertexArray(0);


    // Uncommenting this call will result in wireframe polygons.
    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    float visc = 1.0f;
    float diff = 1.0f;
    float dt = 3.0f;
        
    

    // addForceAwayFromCenter(velocityPrevX, velocityPrevY);
    // void addForce(float *velocityX, float *velocityY, float forceX, float forceY);
    addDensity(densityPrev);

    while (!glfwWindowShouldClose(window))
    {
        // Check if any events have been activiated (key pressed, mouse moved etc.) and call corresponding response functions
        glfwPollEvents();

        // Render
        // Clear the colorbuffer
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // addForce(velocityPrevX, velocityPrevY, 2, 2);

        if (!pauseSim) {
            getFromUI(densityPrev, velocityPrevX, velocityPrevY);
            if (step) {
                vel_step(SIZE, velocityX, velocityY, velocityPrevX, velocityPrevY, visc, dt);
                dens_step(SIZE, density, densityPrev, velocityX, velocityY, diff, dt);            
            }
            setDensity(density);
        }


        // memcpy(density, densityPrev, sizeof(float) * (SIZE + 2) * (SIZE + 2));
        // memcpy(velocityX, velocityPrevX, sizeof(float) * (SIZE + 2) * (SIZE + 2));
        // memcpy(velocityY, velocityPrevY, sizeof(float) * (SIZE + 2) * (SIZE + 2));

        // Draw our first triangle
        shader->use();
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, numRects * 6);
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

void createGrid() {
    float cellSize = 2.0f / SIZE;
    float x = 0.0;
    float y = 1.0;
    for (int i = 0; i <= SIZE + 1; i++) {
        x = -1.0f;
        for (int j = 0; j <= SIZE + 1; j++) {
            data[at(i, j)].id = numRects++;
            data[at(i, j)].hasVelocity = false;
            data[at(i, j)].velocity = glm::vec2(0.0f);
            data[at(i, j)].leftTri.p1 = glm::vec3(x, y, 0.0f);
            data[at(i, j)].leftTri.p2 = glm::vec3(x, y-cellSize, 0.0f);
            data[at(i, j)].leftTri.p3 = glm::vec3(x+cellSize, y, 0.0f);

            data[at(i, j)].rightTri.p1 = glm::vec3(x, y-cellSize, 0.0f);
            data[at(i, j)].rightTri.p2 = glm::vec3(x+cellSize, y, 0.0f);
            data[at(i, j)].rightTri.p3 = glm::vec3(x+cellSize, y-cellSize, 0.0f);

            vertices.push_back(glm::vec3(x, y, 1.0f));
            vertices.push_back(glm::vec3(x, y-cellSize, 1.0f));
            vertices.push_back(glm::vec3(x+cellSize, y, 1.0f));

            vertices.push_back(glm::vec3(x, y-cellSize, 0.0f));
            vertices.push_back(glm::vec3(x+cellSize, y, 0.0f));
            vertices.push_back(glm::vec3(x+cellSize, y-cellSize, 0.0f));

            x += cellSize;
        }
        y -= cellSize;
    }
    std::cout << "\n" << numRects << std::endl << std::endl;
}

void setDensity(float *density) {
    int i, j;
    int N = SIZE;
    FOR_EACH_CELL
        float densityVal = density[at(i,j)];
        vertices[at(i,j)*6].z = densityVal;
        vertices[at(i,j)*6+1].z = densityVal;
        vertices[at(i,j)*6+2].z = densityVal;
        vertices[at(i,j)*6+3].z = densityVal;
        vertices[at(i,j)*6+4].z = densityVal;
        vertices[at(i,j)*6+5].z = densityVal;

        // data[at(i, j)].leftTri.p1.z = 0.0f;
        // data[at(i, j)].leftTri.p2.z = 0.0f;
        // data[at(i, j)].leftTri.p3.z = 0.0f;
        // data[at(i, j)].rightTri.p1.z = 0.0f;
        // data[at(i, j)].rightTri.p2.z = 0.0f;
        // data[at(i, j)].rightTri.p3.z = 0.0f;
    END_FOR
    
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), &vertices[0], GL_STATIC_DRAW);
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

void printVec(glm::vec3 toPrint, std::string vecName) {
    std::cout << vecName << ": ";
    std::cout << toPrint.x << " ";
    std::cout << toPrint.y << " ";
    std::cout << toPrint.z << " ";
    std::cout << std::endl;
}

void addCellDensity(glm::vec3 cellClicked) {
    int index = at(cellClicked.x, cellClicked.y);
    data[index].leftTri.p1.z = 1.0f;
    data[index].leftTri.p2.z = 1.0f;
    data[index].leftTri.p3.z = 1.0f;
    data[index].rightTri.p1.z = 1.0f;
    data[index].rightTri.p2.z = 1.0f;
    data[index].rightTri.p3.z = 1.0f;
}

void addCellVelocity(glm::vec3 cellClicked, glm::vec2 offset) {
    int index = at(cellClicked.x, cellClicked.y);
    for(int i = 0 ; i < SIZE; i++) {
        for(int j = 0 ; j < SIZE; j++) {
            // index = at(i,j);
        }
    }
    // data[index].hasVelocity = true;
    // data[index].velocity = offset;  

    int x = cellClicked.x; 
    int y = cellClicked.y; 
    if (x - 5 > 1 && x + 5 < WIDTH -1) {
        if (y - 5 > 1 && y + 5 < HEIGHT -1) {
            for (int i = x-5; i < x+5; i++) {
                for (int j = y-5; j < y+5; j++) {
                    index = at(i,j);
                    data[index].hasVelocity = true;
                    data[index].velocity = offset;  
                }
            }
        }
    }
}

void addDensity(float *density) {
    for (int i = (SIZE/4); i < SIZE/2 + (SIZE/4); i++) {
        for (int j = (SIZE/4); j < SIZE/2 + (SIZE/4); j++) {
            density[at(i, j)] = 1.0f;
        }   
    }
    // for (int i = 0; i < 40; i++) {
    //     for (int j = 0; j < 40; j++) {
    //         density[at(i, j)] = 1.0f;
    //     }   
    // }
}

void getFromUI(float *densityPrev, float *velocityPrevX, float *velocityPrevY) {

    // std::cout << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;

    for (int i = 1; i <= SIZE; i++) {
        for (int j = 1; j <= SIZE; j++) {
            float right = (data[at(i,j)].rightTri.p1.z + 
                data[at(i,j)].rightTri.p2.z + data[at(i,j)].rightTri.p3.z)/3;
            float left = (data[at(i,j)].leftTri.p1.z + 
                data[at(i,j)].rightTri.p2.z + data[at(i,j)].rightTri.p3.z)/3;
            float densityVal = (left + right) / 2;

            densityPrev[at(i,j)] += densityVal;
            // std::cout << densityVal << " ";
        }        
        // std::cout << std::endl;
    }

    for (int i = 1; i <= SIZE; i++) {
        for (int j = 1; j <= SIZE; j++) {
            if (data[at(i,j)].hasVelocity) {
                velocityPrevX[at(i,j)] = data[at(i,j)].velocity.x;
                velocityPrevY[at(i,j)] = data[at(i,j)].velocity.y;
            }
        }        
    }
}

void resetDensity() {
    for (int i = 1; i <= SIZE; i++) {
        for (int j = 1; j <= SIZE; j++) {
            data[at(i,j)].leftTri.p1.z = 0.0f;
            data[at(i,j)].leftTri.p2.z = 0.0f;
            data[at(i,j)].leftTri.p3.z = 0.0f;
            data[at(i,j)].rightTri.p1.z = 0.0f;
            data[at(i,j)].rightTri.p2.z = 0.0f;
            data[at(i,j)].rightTri.p3.z = 0.0f;
        }        
    }
}




// MOUSE
// Is called whenever a key is pressed/released via GLFW
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
    if (key == GLFW_KEY_S && action == GLFW_PRESS) {
        step = !step;
    }
    if (key == GLFW_KEY_P && action == GLFW_PRESS) {
        pauseSim = !pauseSim;
    }
    if (key == GLFW_KEY_R && action == GLFW_PRESS) {
        resetDensity();
    }
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos){
    if(firstMouse)
    {
        currX = xpos;
        currY = ypos;
        firstMouse = false;
    }
    // else {
    //     lastX = currX;
    //     lastY = currY;
    //     currX = xpos;
    //     currY = xpos;
    // }
    else {
        if (xpos <= WIDTH && xpos >=0 && ypos <= HEIGHT && ypos >=0) {
            GLfloat xoffset = xpos - lastX;
            GLfloat yoffset = lastY - ypos; 
            lastX = xpos;
            lastY = ypos;

            glm::vec2 offset = glm::normalize(glm::vec2(xoffset, -yoffset));
            glm::vec3 cellClicked = getCell(glm::vec3(lastX, lastY, 0.0f));
            if (go) {
                addCellDensity(cellClicked);
                printVec(cellClicked, "cellClicked");            
            }
            else {            
                addCellVelocity(cellClicked, offset);
                printVec(glm::vec3(offset, 0.0f), "offset");            
            }
        }
    }
} 

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    if (state == GLFW_PRESS) {
        go = true;

    }
    else {
        go = false;
    }   
}









void add_source(int N, float *x, float *s, float dt)
{
int i, size=(N+2)*(N+2);
for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i];
}

void set_bnd(int N, int b, float *x)
{
int i;
for ( i=1 ; i<=N ; i++ ) {
x[IX(0 ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
x[IX(i,0 )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
}
x[IX(0 ,0 )] = 0.5f*(x[IX(1,0 )]+x[IX(0 ,1)]);
x[IX(0 ,N+1)] = 0.5f*(x[IX(1,N+1)]+x[IX(0 ,N)]);
x[IX(N+1,0 )] = 0.5f*(x[IX(N,0 )]+x[IX(N+1,1)]);
x[IX(N+1,N+1)] = 0.5f*(x[IX(N,N+1)]+x[IX(N+1,N)]);
}
void lin_solve(int N, int b, float *x, float *x0,
float a, float c)
{
int i, j, n;
for ( n=0 ; n<20 ; n++ ) {
FOR_EACH_CELL
x[IX(i,j)] = (x0[IX(i,j)]+a*(x[IX(i-1,j)]+
x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;
END_FOR
set_bnd ( N, b, x );
}
}
void diffuse(int N, int b, float *x, float *x0,
float diff, float dt)
{
float a=dt*diff*N*N;
lin_solve ( N, b, x, x0, a, 1+4*a );
}

void advect(int N, int b, float *d, float *d0, float *u, float *v, float dt)
{
int i, j, i0, j0, i1, j1;
float x, y, s0, t0, s1, t1, dt0;
dt0 = dt*N;
FOR_EACH_CELL
x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];
if (x<0.5f) x=0.5f; if (x>N+0.5f) x=N+0.5f; i0=(int)x; i1=i0+1;
if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
END_FOR
set_bnd ( N, b, d );
}
void project(int N, float * u, float * v, float * p, float * div)
{
int i, j;
FOR_EACH_CELL
div[IX(i,j)] = -0.5f*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
p[IX(i,j)] = 0;
END_FOR
set_bnd ( N, 0, div ); set_bnd ( N, 0, p );
lin_solve ( N, 0, p, div, 1, 4 );
FOR_EACH_CELL
u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
END_FOR
set_bnd ( N, 1, u ); set_bnd ( N, 2, v );
}
void dens_step(int N, float *x, float *x0, float *u, float *v, float diff, float dt)
{
add_source ( N, x, x0, dt );
SWAP ( x0, x ); diffuse ( N, 0, x, x0, diff, dt );
SWAP ( x0, x ); advect ( N, 0, x, x0, u, v, dt );
}
void vel_step(int N, float *u, float *v, float *u0, float *v0, float visc, float dt)
{
add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt );
SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt );
project ( N, u, v, u0, v0 );
SWAP ( u0, u ); SWAP ( v0, v );
advect ( N, 1, u, u0, u0, v0, dt ); advect ( N, 2, v, v0, u0, v0, dt );
project ( N, u, v, u0, v0 );
}

// void add_source ( int N, float * x, float * s, float dt )
// {
//     int i, size=(N+2)*(N+2);
//     for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i];
// }

// void set_bnd ( int N, int b, float * x )
// {
//     int i;

//     for ( i=1 ; i<=N ; i++ ) {
//         x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
//         x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
//         x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
//         x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
//     }
//     x[IX(0  ,0  )] = 0.5f*(x[IX(1,0  )]+x[IX(0  ,1)]);
//     x[IX(0  ,N+1)] = 0.5f*(x[IX(1,N+1)]+x[IX(0  ,N)]);
//     x[IX(N+1,0  )] = 0.5f*(x[IX(N,0  )]+x[IX(N+1,1)]);
//     x[IX(N+1,N+1)] = 0.5f*(x[IX(N,N+1)]+x[IX(N+1,N)]);
// }

// void lin_solve ( int N, int b, float * x, float * x0, float a, float c )
// {
//     int i, j, k;

//     for ( k=0 ; k<20 ; k++ ) {
//         FOR_EACH_CELL
//             x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;
//         END_FOR
//         set_bnd ( N, b, x );
//     }
// }

// void diffuse ( int N, int b, float * x, float * x0, float diff, float dt )
// {
//     float a=dt*diff*N*N;
//     lin_solve ( N, b, x, x0, a, 1+4*a );
// }

// void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt )
// {
//     int i, j, i0, j0, i1, j1;
//     float x, y, s0, t0, s1, t1, dt0;

//     dt0 = dt*N;
//     FOR_EACH_CELL
//         x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];
//         if (x<0.5f) x=0.5f; if (x>N+0.5f) x=N+0.5f; i0=(int)x; i1=i0+1;
//         if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
//         s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
//         d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
//                      s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
//     END_FOR
//     set_bnd ( N, b, d );
// }

// void project ( int N, float * u, float * v, float * p, float * div )
// {
//     int i, j;

//     FOR_EACH_CELL
//         div[IX(i,j)] = -0.5f*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
//         p[IX(i,j)] = 0;
//     END_FOR 
//     set_bnd ( N, 0, div ); set_bnd ( N, 0, p );

//     lin_solve ( N, 0, p, div, 1, 4 );

//     FOR_EACH_CELL
//         u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
//         v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
//     END_FOR
//     set_bnd ( N, 1, u ); set_bnd ( N, 2, v );
// }

// void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt )
// {
//     add_source ( N, x, x0, dt );
//     SWAP ( x0, x ); diffuse ( N, 0, x, x0, diff, dt );
//     SWAP ( x0, x ); advect ( N, 0, x, x0, u, v, dt );
// }

// void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt )
// {
//     add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
//     SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt );
//     SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt );
//     project ( N, u, v, u0, v0 );
//     SWAP ( u0, u ); SWAP ( v0, v );
//     advect ( N, 1, u, u0, u0, v0, dt ); advect ( N, 2, v, v0, u0, v0, dt );
//     project ( N, u, v, u0, v0 );
// }
// // FLUID SOLVER
// void add_source(int N, float *x, float *s, float dt) {
//     int size = (N + 2) * (N + 2);
//     for (int i = 0; i < size; i++) {
//         x[i] += dt * s[i];
//     }
// }

// void set_bnd(int N, int b, float * x ) {
//     int i;
//     for ( i=1 ; i<=N ; i++ ) {
//         x[at(0, i)]      = b==1 ? -x[at(1, i)] : x[at(1,i)];
//         x[at(N+1, i)]    = b==1 ? -x[at(N, i)] : x[at(N,i)];
//         x[at(i, 0)]      = b==2 ? -x[at(i, 1)] : x[at(i, 1)];
//         x[at(i, N+1)]    = b==2 ? -x[at(i, N)] : x[at(i,N)];
//     }
//     x[at(0 ,0 )] = 0.5*(x[at(1,0 )]+x[at(0 ,1)]);
//     x[at(0 ,N+1)] = 0.5*(x[at(1,N+1)]+x[at(0 ,N )]);
//     x[at(N+1,0 )] = 0.5*(x[at(N,0 )]+x[at(N+1,1)]);
//     x[at(N+1,N+1)] = 0.5*(x[at(N,N+1)]+x[at(N+1,N)]);
// }

// void diffuse(int N, int b, float *x, float *x0, float diff, float dt) {
//     int i, j, k;
//     float a = dt * diff * N * N;

//     for (k = 0; k < 20; k++) {
//         for(int i = 1; i <= N; i++) {
//             for (int j = 0; j <= N; j++) {
//                 float surrDensity = x[at(i-1,j)] + x[at(i+1,j)] + x[at(i,j-1)] + x[at(i,j+1)];
//                 x[at(i,j)] = (x0[at(i,j)] + a * surrDensity) / (1 + 4*a);
//             }
//         }
//         set_bnd(N,b,x);
//     }
// }

// void advect(int N, int b, float *d, float *d0, float *u, float *v, float dt) {
//     int i, j, i0, j0, i1, j1;
//     float x, y, s0, t0, s1, t1, dt0;

//     dt0 = dt * N;
//     for(int i = 1; i <= N; i++) {
//         for (int j = 0; j <= N; j++) {
//             x = i-dt0*u[at(i,j)]; y = j-dt0*v[at(i,j)];
//             if (x<0.5) x=0.5; if (x>N+0.5) x=N+0.5; i0=(int)x; i1=i0+1;
//             if (y<0.5) y=0.5; if (y>N+0.5) y=N+0.5; j0=(int)y; j1=j0+1; 
//             s1 = x-i0; s0 = 1-s1; 
//             t1 = y-j0; t0 = 1-t1;
//             d[at(i,j)] = s0*(t0*d0[at(i0,j0)]+t1*d0[at(i0,j1)])+s1*(t0*d0[at(i1,j0)]+t1*d0[at(i1,j1)]);
//         }
//     }
//     set_bnd(N, b, d);
// }

// void project ( int N, float * u, float * v, float * p, float * div ) {
//     int i, j, k;
//     float h;

//     h = 1.0 / N;
//     for(int i = 1; i <= N; i++) {
//         for (int j = 0; j <= N; j++) {
//             div[at(i,j)] = -0.5*h*(u[at(i+1,j)]-u[at(i-1,j)]+v[at(i,j+1)]-v[at(i,j-1)]);
//             p[at(i,j)] = 0;
//         }
//     }

//     set_bnd(N, 0, div);
//     set_bnd(N, 0, p);

//     for (k = 0; k < 20; k++) {
//         for(int i = 1; i <= N; i++) {
//             for (int j = 0; j <= N; j++) {
//                 p[at(i,j)] = (div[at(i,j)]+p[at(i-1,j)]+p[at(i+1,j)]+p[at(i,j-1)]+p[at(i,j+1)])/4;
//             }
//         }
//         set_bnd(N, 0, p);
//     }

//     for(int i = 1; i <= N; i++) {
//         for (int j = 0; j <= N; j++) {
//         u[at(i,j)] -= 0.5*(p[at(i+1,j)]-p[at(i-1,j)])/h;
//         v[at(i,j)] -= 0.5*(p[at(i,j+1)]-p[at(i,j-1)])/h;
//     }
// }
//     set_bnd(N, 1, u); 
//     set_bnd(N, 2, v);
// }

// void dens_step(int N, float * x, float * x0, float * u, float * v, float diff, float dt) {
//     add_source(N, x, x0, dt);
//     SWAP(x0, x); 
//     diffuse(N, 0, x, x0, diff, dt);
//     SWAP(x0, x);
//     advect(N, 0, x, x0, u, v, dt);
// }

// void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt ) {
//     add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
//     SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt );
//     SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt );
//     project ( N, u, v, u0, v0 );
//     SWAP ( u0, u ); SWAP ( v0, v );
//     advect ( N, 1, u, u0, u0, v0, dt ); advect ( N, 2, v, v0, u0, v0, dt ); project ( N, u, v, u0, v0 );
// }


// //adding sources to scene
void addForce(float *velocityX, float *velocityY, float forceX, float forceY) {
    // for (int i = 50; i < 100; i++) {
    //     for (int j = 50; j < 100; j++) {
    //         velocityX[at(i,j)] = 0.2f;
    //         velocityY[at(i,j)] = 0.2f;
    //     }
    // }
    int i, j;
    int N = SIZE;
    FOR_EACH_CELL

        velocityX[at(i,j)] = forceX;
        velocityY[at(i,j)] = forceY;
    END_FOR
}

void addForceAwayFromCenter(float *velocityX, float *velocityY) {
    glm::vec3 center = glm::vec3(SIZE / 2);

    for (int i = START; i < END; i++) {
        for (int j = START; j < END; j++) {
            glm::vec3 pt = glm::vec3(i, j, 0.0f);
            glm::vec3 force = glm::normalize(pt - center);
            velocityX[at(i,j)] = force.x;
            velocityY[at(i,j)] = force.y;
        }
    }
}

// // Drawing the density
// void drawDensity(float *density) {
//     color_t color;
//     color.r = 0.0f;
//     color.g = 0.0f;
//     color.b = 0.0f;
//     color.f = 1.0f;
//     std::cout << std::endl;
//     for (int j = START; j < END; j++) {
//         for (int i = START; i < END; i++) {
//             float val = density[at(i,j)];
//             color.r = val;
//             color.g = val;
//             color.b = val;
//             img.pixel(i, j, color);
//         }
//     }
//     img.WriteTga((char *) "density.tga", true);
//     system("open density.tga");
//     usleep(TIME);
//     system("rm density.tga");
// }

// //debugging
// void printArr(float *arr, std::string arrName) {
//     std::cout << arrName << ": " << std::endl;
//     for (int j = START; j < END; j++) {
//         std::cout << "j = "<< j << ": ";
//         for (int i = START; i < END; i++) {
//             std::cout << arr[at(i, j)] << " ";
//         }
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;
// }