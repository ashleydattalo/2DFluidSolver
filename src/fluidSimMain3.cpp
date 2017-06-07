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
GLuint posBufID, densBufID;
std::vector<glm::vec3> vertices;
std::vector<GLfloat> densities;

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
    bool hasDensity;
    float density;
};

Rect data[(SIZE + 2) * (SIZE + 2)];
int numRects = 0;

bool step = true;
bool pauseSim = false;
float densityColor = 1.0f;

float visc = 1.0f;
float diff = 1.0f;
float dt = 8.0f;

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

    glDisable(GL_DEPTH_TEST);

    // Define the viewport dimensions
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);  
    glViewport(0, 0, width, height);
    
    // GLuint bufs[2];
    // glGenBuffers(2, bufs);
    // posBufID = bufs[0];
    // densBufID = bufs[1];

    glGenBuffers(1, &VBO);
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), &vertices[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (GLvoid*)0);

    // glEnableVertexAttribArray(0);
    // glBindBuffer(GL_ARRAY_BUFFER, posBufID);
    // glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec2), &vertices[0], GL_STATIC_DRAW);
    // glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0);

    // glEnableVertexAttribArray(1);
    // glBindBuffer(GL_ARRAY_BUFFER, densBufID);
    // glBufferData(GL_ARRAY_BUFFER, densities.size() * sizeof(GLfloat), &densities[0], GL_DYNAMIC_DRAW);
    // glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid*)(0));

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // Uncommenting this call will result in wireframe polygons.
    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        
    addDensity(densityPrev);

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // addForce(velocityPrevX, velocityPrevY, 2, 2);
        // addForceAwayFromCenter(velocityPrevX, velocityPrevY);

        if (!pauseSim) {
            getFromUI(densityPrev, velocityPrevX, velocityPrevY);
            if (step) {
                vel_step(SIZE, velocityX, velocityY, velocityPrevX, velocityPrevY, visc, dt);
                dens_step(SIZE, density, densityPrev, velocityX, velocityY, diff, dt);            
            }
            setDensity(density);
        }

        // memcpy(densityPrev, density, sizeof(float) * (SIZE + 2) * (SIZE + 2));
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
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
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

            data[at(i, j)].hasDensity = false;
            data[at(i, j)].density = 0.0f;

            vertices.push_back(glm::vec3(x, y, 0.0f));
            vertices.push_back(glm::vec3(x, y-cellSize, 0.0f));
            vertices.push_back(glm::vec3(x+cellSize, y, 0.0f));

            vertices.push_back(glm::vec3(x, y-cellSize, 0.0f));
            vertices.push_back(glm::vec3(x+cellSize, y, 0.0f));
            vertices.push_back(glm::vec3(x+cellSize, y-cellSize, 0.0f));

            // vertices.push_back(glm::vec2(x, y));
            // vertices.push_back(glm::vec2(x, y-cellSize));
            // vertices.push_back(glm::vec2(x+cellSize, y));

            // vertices.push_back(glm::vec2(x, y-cellSize));
            // vertices.push_back(glm::vec2(x+cellSize, y));
            // vertices.push_back(glm::vec2(x+cellSize, y-cellSize));

            densities.push_back(0.0f);
            densities.push_back(0.0f);
            densities.push_back(0.0f);
            densities.push_back(0.0f);
            densities.push_back(0.0f);
            densities.push_back(0.0f);

            x += cellSize;
        }
        y -= cellSize;
    }
    std::cout << "\n" << numRects << std::endl << std::endl;
}

void getFromUI(float *densityPrev, float *velocityPrevX, float *velocityPrevY) {
    for (int i = 1; i <= SIZE; i++) {
        for (int j = 1; j <= SIZE; j++) {
            if (data[at(i,j)].hasDensity) {
                densityPrev[at(i,j)] = data[at(i,j)].density;
                // data[at(i,j)].hasDensity = false;
                // data[at(i,j)].density = 0.0f;
            }
        }        
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

void setDensity(float *density) {
    int i, j;
    int N = SIZE;
    FOR_EACH_CELL
        int index = at(i,j);
        float densityVal = density[index];
        // densities[index] = densityVal;
        // densities[index+1] = densityVal;
        // densities[index+2] = densityVal;
        // densities[index+3] = densityVal;
        // densities[index+4] = densityVal;
        // densities[index+5] = densityVal;
        vertices[at(i,j)*6].z = densityVal;
        vertices[at(i,j)*6+1].z = densityVal;
        vertices[at(i,j)*6+2].z = densityVal;
        vertices[at(i,j)*6+3].z = densityVal;
        vertices[at(i,j)*6+4].z = densityVal;
        vertices[at(i,j)*6+5].z = densityVal;

        // data[at(i,j)].hasDensity = true;
        // data[at(i,j)].density += densityVal;
        // data[at(i,j)].hasDensity = false;
        // data[at(i,j)].density = 0.0f;

        // data[at(i,j)].hasVelocity = false;
        // data[at(i,j)].velocity = glm::vec2(0.0f);
    END_FOR

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), &vertices[0], GL_DYNAMIC_DRAW);
}

void addCellDensity(glm::vec3 cellClicked) {
    int index = at(cellClicked.x, cellClicked.y);

    int x = cellClicked.x;
    int y = cellClicked.y;
    data[index].hasDensity = true;
    data[index].density += densityColor;

    if (x - 2 > 0 && x + 2 < WIDTH) {
        if (y - 2 > 0 && y + 2 < HEIGHT) {
            for (int i = -2; i < 2; i++) {
                for (int j = -2; j < 2; j++) {
                    index = at(x+i,j+y);
                    data[index].hasDensity = true;
                    data[index].density += densityColor;
                }    
            }
        }
    }
}

void addCellVelocity(glm::vec3 cellClicked, glm::vec2 offset) {
    int index = at(cellClicked.x, cellClicked.y);
    data[index].hasVelocity = true;
    data[index].velocity = offset; 

    // int x = cellClicked.x;
    // int y = cellClicked.y;

    // if (x - 5 > 0 && x + 5 < WIDTH) {
    //     if (y - 5 > 0 && y + 5 < HEIGHT) {
    //         for (int i = -5; i < 5; i++) {
    //             for (int j = -5; j < 5; j++) {
    //                 index = at(x+i,j+y);
                    
    //                 float dist = glm::sqrt((i)*(i) + (j)*(j));
    //                 glm::vec2 force = dist * offset;

    //                 data[index].hasVelocity = true;
    //                 data[index].velocity = force;
    //             }    
    //         }
    //     }
    // } 
}

void addDensity(float *density) {
    for (int i = (SIZE/4); i < SIZE/2 + (SIZE/4); i++) {
        for (int j = (SIZE/4); j < SIZE/2 + (SIZE/4); j++) {
            density[at(i, j)] = 1.0f;
        }   
    }
}

void resetDensity() {
    for (int i = 1; i <= SIZE; i++) {
        for (int j = 1; j <= SIZE; j++) {
            data[at(i,j)].hasDensity = false;
            data[at(i,j)].density = 0.0f;
        }        
    }
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
    if (key == GLFW_KEY_K && action == GLFW_PRESS) {
        densityColor += 0.1f;
    }
    if (key == GLFW_KEY_M && action == GLFW_PRESS) {
        densityColor -= 0.1f;
    }
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos){
    if(firstMouse)
    {
        currX = xpos;
        currY = ypos;
        firstMouse = false;
    }
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







// SOLVER STUFF

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

// //adding sources to scene
void addForce(float *velocityX, float *velocityY, float forceX, float forceY) {
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