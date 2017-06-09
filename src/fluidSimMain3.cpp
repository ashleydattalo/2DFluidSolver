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
#define START 1
#define END SIZE + 1
#define at(i,j) ((i) + (SIZE + 2)*(j))
#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) {\
for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}

void dens_step(int N, float * x, float * x0, float * u, float * v, float diff, float dt);
void vel_step(int N, float * u, float * v, float * u0, float * v0, float visc, float dt);

void addDensity(float *density);
void addForce(float *velocityX, float *velocityY, float forceX, float forceY);
void addForceAwayFromCenter(float *velocityX, float *velocityY);
void printArr(float *arr, std::string arrName);
void drawDensity(float *density);
Image img(SIZE + 2, SIZE + 2);

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
void processInput();
bool inBounds(glm::vec2 mousePos);

glm::vec2 getCell(glm::vec3 mousePos);
void addDensityToCell(glm::vec2 mousePos);
void addVelocityToCell(glm::vec2 newMousePos, glm::vec2 lastMousePos);

void createGrid();
void getFromUI(float *densityPrev, float *velocityPrevX, float *velocityPrevY);
void setDensity(float *density);

void initWindow();
void bindBuffers();
void drawBunny(glm::vec2 offset, float bunnySize);
void printVec(glm::vec3 toPrint, std::string vecName);
void printVec(glm::vec2 toPrint, std::string vecName);

GLFWwindow* window;
Shader *shader;

const GLuint WIDTH = 800, HEIGHT = 600;

GLuint VBO, VAO;
GLuint posBufID, densBufID;
GLuint elementbuffer;
std::vector<glm::vec3> vertices;
std::vector<unsigned int> indices;

struct Rect {
    float density;
    glm::vec2 velocity;
};

struct IndicesStruct {
    int idx1;
    int idx2;
    int idx3;
    int idx4;
};

Rect data[(SIZE + 2) * (SIZE + 2)];
IndicesStruct indicesArr[SIZE][SIZE];

bool step = true;
bool pauseSim = false;
float densityStrength = 100.0f;
float force = 10.0f;

float visc = 0.0f;
float diff = 0.0f;
float dt = 0.1f;

float *density;
float *velocityX;
float *velocityY;
float *densityPrev;
float *velocityPrevX;
float *velocityPrevY;

glm::vec3 colorScale(0.1f, 0.08, 0.42f);
bool pauseColors = false;

int main()
{
    //allocate data
    density = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
    velocityX = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
    velocityY = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

    densityPrev = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
    velocityPrevX = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
    velocityPrevY = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

    memset(densityPrev, 0, sizeof(float) * (SIZE + 2) * (SIZE + 2));
    memset(velocityPrevX, 0, sizeof(float) * (SIZE + 2) * (SIZE + 2));
    memset(velocityPrevY, 0, sizeof(float) * (SIZE + 2) * (SIZE + 2));

    initWindow();
    createGrid();
    bindBuffers();

    shader = new Shader(SHADER_PATH "waterSim.vert", SHADER_PATH "waterSim.frag");

    // Uncommenting this call will result in wireframe polygons.
    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();
        processInput();

        GLint scalePos = glGetUniformLocation(shader->getProg(), "scale");
        glUniform3f(scalePos, colorScale.x, colorScale.y, colorScale.z);  

        if (!pauseSim) {
            getFromUI(densityPrev, velocityPrevX, velocityPrevY);
            // addForceAwayFromCenter(velocityPrevX, velocityPrevY);
            if (step) {
                vel_step(SIZE, velocityX, velocityY, velocityPrevX, velocityPrevY, visc, dt);
                dens_step(SIZE, density, densityPrev, velocityX, velocityY, diff, dt);            
            }
            setDensity(density);
        }

        shader->use();
        glBindVertexArray(VAO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
        glDrawElements(GL_TRIANGLE_STRIP, indices.size(), GL_UNSIGNED_INT, 0);
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
    float left = -1.0f, top = 1.0f;
    float x = left, y = top;
    
    // create vertices
    for (int j = 0; j <= SIZE; j++) {
        x = left;
        for (int i = 0; i <= SIZE; i++) {
            vertices.push_back(glm::vec3(x, y, 1.0f));
            x += cellSize;
        }
        y -= cellSize;
    }

    // create indices
    int sizePlus1 = SIZE + 1;
    for (int j = 0; j < SIZE; j++) {
        for (int i = 0; i < SIZE; i++) {
            indices.push_back(j * sizePlus1 + i);
            indices.push_back(j * sizePlus1 + i + 1);
            indices.push_back( (j + 1) * sizePlus1 + i);

            indices.push_back(j * sizePlus1 + i + 1);
            indices.push_back((j + 1) * sizePlus1 + i);
            indices.push_back((j + 1) * sizePlus1 + i + 1);

            indicesArr[i][j].idx1 = j * sizePlus1 + i;
            indicesArr[i][j].idx2 = j * sizePlus1 + i + 1;
            indicesArr[i][j].idx3 = (j + 1) * sizePlus1 + i;
            indicesArr[i][j].idx4 = (j + 1) * sizePlus1 + i + 1;

        }
    }
}

void getFromUI(float *densityPrev, float *velocityPrevX, float *velocityPrevY) {
    for (int i = 1; i <= SIZE; i++) {
        for (int j = 1; j <= SIZE; j++) {
            densityPrev[at(i,j)] = data[at(i,j)].density;
        }        
    }

    for (int i = 1; i <= SIZE; i++) {
        for (int j = 1; j <= SIZE; j++) {
            velocityPrevX[at(i,j)] = data[at(i,j)].velocity.x;
            velocityPrevY[at(i,j)] = data[at(i,j)].velocity.y;
        }        
    }
}

void setDensity(float *density) {
    int i, j;
    int N = SIZE;
    FOR_EACH_CELL
        int index = at(i,j);
        float densityVal = density[index];

        int idx1 = indicesArr[i-1][j-1].idx1;
        int idx2 = indicesArr[i-1][j-1].idx2;
        int idx3 = indicesArr[i-1][j-1].idx3;
        int idx4 = indicesArr[i-1][j-1].idx4;

        vertices[idx1].z = densityVal;
        vertices[idx2].z = densityVal;
        vertices[idx3].z = densityVal;

        vertices[idx2].z = densityVal;
        vertices[idx3].z = densityVal;
        vertices[idx4].z = densityVal;

        data[at(i,j)].density = 0.0f;
        data[at(i,j)].velocity = glm::vec2(0.0f);
    END_FOR

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), &vertices[0], GL_DYNAMIC_DRAW);
}

glm::vec2 getCell(glm::vec2 mousePos) {
    glm::vec2 cell;
    float screenCellSizeX = WIDTH / SIZE;
    float screenCellSizeY = HEIGHT / SIZE;

    int xCoord = mousePos.x / screenCellSizeX;
    int yCoord = mousePos.y / screenCellSizeY;
    
    cell.x = xCoord + 1;
    cell.y = yCoord + 1;
    
    return cell;
}

void addDensityToCell(glm::vec2 mousePos) {
    glm::vec2 cell = getCell(mousePos);
    int index = at(cell.x, cell.y);

    int x = cell.x;
    int y = cell.y;

    if (x - 2 > 0 && x + 2 < WIDTH) {
        if (y - 2 > 0 && y + 2 < HEIGHT) {
            for (int i = -2; i < 2; i++) {
                for (int j = -2; j < 2; j++) {
                    index = at(x+i,j+y);
                    data[index].density = densityStrength;
                }    
            }
        }
    }
}

void deleteDensityFromCell(glm::vec2 mousePos) {
    glm::vec2 cell = getCell(mousePos);
    int index = at(cell.x, cell.y);

    int x = cell.x;
    int y = cell.y;

    int offset = 3;
    if (inBounds(mousePos + glm::vec2(offset)) && inBounds(mousePos - glm::vec2(offset))) {   
        for (int i = -offset; i < offset; i++) {
            for (int j = -offset; j < offset; j++) {
                if (x+i <= SIZE && j+y <= SIZE) {
                    index = at(x+i,j+y);
                    data[index].density -= densityStrength;
                }
            }    
        }
    }
}

void addVelocityToCell(glm::vec2 newMousePos, glm::vec2 lastMousePos) {
    glm::vec2 cell = getCell(newMousePos);
    int index = at(cell.x, cell.y);
    
    glm::vec2 offset = newMousePos - lastMousePos;
    data[index].velocity = force * offset;  
}

void addPointForce(glm::vec2 mousePos) {
    glm::vec2 cell = getCell(mousePos);
    int x = cell.x;
    int y = cell.y;
    int offset = 10;
    float rad = 5;
    if (inBounds(mousePos + glm::vec2(offset)) && inBounds(mousePos - glm::vec2(offset))) {   
        for (int i = -offset; i <= offset; i++) {
            for (int j = -offset; j <= offset; j++) {
                if (x+i <= SIZE && j+y <= SIZE) {
                    float pointForce = glm::distance(cell, glm::vec2(x+i,j+y));
                    glm::vec2 vel = pointForce * (glm::vec2(x+i,j+y) - cell);
                    int index = at(x+i,j+y);
                    if (pointForce < rad)
                        data[index].velocity = vel;
                }
            }    
        }
    }
}


bool keys[1024]; 
bool firstMouse = true;
bool mouseClicked = false;
bool mouseDown = false;
glm::vec2 currMousePos, oldMousePos;

// Is called whenever a key is pressed/released via GLFW
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);

    if(action == GLFW_PRESS)
        keys[key] = true;
    else if(action == GLFW_RELEASE)
        keys[key] = false; 
}

// called when ever the mouse moves
void mouse_callback(GLFWwindow* window, double xpos, double ypos){
    if(firstMouse)
    {
        currMousePos.x = xpos;
        currMousePos.y = ypos;
        firstMouse = false;
    }
    oldMousePos = currMousePos;
    currMousePos.x = xpos;
    currMousePos.y = ypos;
} 

// called when ever the mouse is clicked
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    if (state == GLFW_PRESS) {
        mouseClicked = true;
        mouseDown = true;
    }
    else if (state == GLFW_RELEASE) {
        mouseDown = false;
    }
}

void processInput() {
    if (inBounds(currMousePos) && inBounds(oldMousePos)) {
        // adds density
        if (keys[GLFW_KEY_F]) {
            addDensityToCell(currMousePos);
        }
        // adds velocity
        if (keys[GLFW_KEY_V]) {
            addVelocityToCell(currMousePos, oldMousePos);
        }
        // eraser
        if (keys[GLFW_KEY_D]) {
            deleteDensityFromCell(currMousePos);
        }
        // draws bunny
        if (keys[GLFW_KEY_B]) {
            drawBunny(getCell(currMousePos + 50.0f), 2);
        }
        // pause
        if (keys[GLFW_KEY_P]) {
            pauseSim = !pauseSim;
        }

        // changes force
        if (keys[GLFW_KEY_A]) {
            if(keys[GLFW_KEY_LEFT_SHIFT]) {
                force--;
            }
            else {
                force++;
            }
            std::cout << "force: " << force << std::endl;
        }

        // changes r
        if (keys[GLFW_KEY_Q]) {
            if(keys[GLFW_KEY_LEFT_SHIFT]) {
                colorScale.x -= 0.001f;
            }
            else {
                colorScale.x += 0.001f;
            }
            printVec(colorScale, "colorScale");
        }

        // changes g
        if (keys[GLFW_KEY_W]) {
            if(keys[GLFW_KEY_LEFT_SHIFT]) {
                colorScale.y -= 0.001f;
            }
            else {
                colorScale.y += 0.001f;
            }
            printVec(colorScale, "colorScale");
        }

        // changes b
        if (keys[GLFW_KEY_E]) {
            if(keys[GLFW_KEY_LEFT_SHIFT]) {
                colorScale.z -= 0.001f;
            }
            else {
                colorScale.z += 0.001f;
            }
            printVec(colorScale, "colorScale");
        }

        // stops changing the colors
        if (keys[GLFW_KEY_S]) {
            if(keys[GLFW_KEY_LEFT_SHIFT]) {
                pauseColors = false;
            }
            else {
                pauseColors = true;
            }
        }

        if (keys[GLFW_KEY_P]) {
            addPointForce(currMousePos);
        }
    }
}

bool inBounds(glm::vec2 mousePos) {
    int x = mousePos.x;
    int y = mousePos.y;
    return x <= WIDTH && x >=0 && y <= HEIGHT && y>=0;
}

bool inArrBounds(glm::vec2 index) {
    int x = index.x;
    int y = index.y;
    return  x <= SIZE && x >=1 && y <= SIZE && y>=1; 
}

// OBJ
void drawBunny(glm::vec2 offset, float bunnySize) {
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
            glm::vec2 vertex;
            float random;
            fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &random);
            vertex.y *= -1;
            vertex *= 100 * bunnySize;
            vertex += offset;

            if (inArrBounds(vertex)) {
                for (int i = -2; i < 2; i++) {
                    for (int j = -2; j < 2; j++) {
                        glm::vec2 newIndex = glm::vec2(vertex.x+i, vertex.y+j);
                        if (inArrBounds(newIndex)) {
                            data[at((int)newIndex.x, (int)newIndex.y)].density = 10.0f;
                        }
                    }    
                }
            }
        }
        if ( strcmp( lineHeader, "vn" ) == 0 ){
            glm::vec3 normal;
            fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z );
        }
    }
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
            data[at(i,j)].density = 0.0f;
        }        
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
void drawDensity(float *density) {
    color_t color;
    color.r = 0.0f;
    color.g = 0.0f;
    color.b = 0.0f;
    color.f = 1.0f;
    std::cout << std::endl;
    for (int j = START; j < END; j++) {
        for (int i = START; i < END; i++) {
            float val = density[at(i,j)];
            color.r = val;
            color.g = val;
            color.b = val;
            img.pixel(i, j, color);
        }
    }
    img.WriteTga((char *) "density.tga", true);
    system("open density.tga");
    usleep(TIME);
    system("rm density.tga");
}

// //debugging
void printArr(float *arr, std::string arrName) {
    std::cout << arrName << ": " << std::endl;
    for (int j = START; j < END; j++) {
        std::cout << "j = "<< j << ": ";
        for (int i = START; i < END; i++) {
            std::cout << arr[at(i, j)] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void bindBuffers() {
    glDisable(GL_DEPTH_TEST);
    glGenBuffers(1, &VBO);
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), &vertices[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (GLvoid*)0);

    // Generate a buffer for the indices
    glGenBuffers(1, &elementbuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, elementbuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void initWindow() {
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
    }
    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    glewExperimental = GL_TRUE;

    if (glewInit() != GLEW_OK)
    {
        std::cout << "Failed to initialize GLEW" << std::endl;
    }    

    // Define the viewport dimensions
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);  
    glViewport(0, 0, width, height);
}

void printVec(glm::vec3 toPrint, std::string vecName) {
    std::cout << vecName << ": ";
    std::cout << toPrint.x << " ";
    std::cout << toPrint.y << " ";
    std::cout << toPrint.z << " ";
    std::cout << std::endl;
}

void printVec(glm::vec2 toPrint, std::string vecName) {
    std::cout << vecName << ": ";
    std::cout << toPrint.x << " ";
    std::cout << toPrint.y << " ";
    std::cout << std::endl;
}


