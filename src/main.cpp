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
bool inScreenBounds(glm::vec2 mousePos);

glm::vec2 getCell(glm::vec3 mousePos);
void addDensityToCell(glm::vec2 cell, float densityAmount);
void addVelocityToCell(glm::vec2 newMousePos, glm::vec2 lastMousePos);
void addSquareDensity(glm::vec2 mousePos);
bool inScreenBounds(glm::vec2 mousePos);
bool inArrBounds(glm::vec2 index);

void createGrid();
void getFromUI(float *densityPrev, float *velocityPrevX, float *velocityPrevY);
void setDensity(float *density);

void initWindow();
void bindBuffers();

void createBunnyVertices();
void drawBunny(glm::vec2 offset, float bunnySize, float densStrength);

void printVec(glm::vec3 toPrint, std::string vecName);
void printVec(glm::vec2 toPrint, std::string vecName);
void debugMiniGrid(glm::vec2 cell, int offset, std::string name);

void processColors();

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
float densityStrength = 10.0f;
float force = 1.0f;

float visc = 0.0f;
float diff = 0.0f;
float dt = 0.1f;

float *density;
float *velocityX;
float *velocityY;
float *densityPrev;
float *velocityPrevX;
float *velocityPrevY;

std::vector<glm::vec2> bunnyVerts;

glm::vec3 densityColor(0.1f, 0.08, 0.42f);
// glm::vec3 origColorScale(0.0008f, 0.001f, 0.002f);
glm::vec3 origColorScale(0.0003f, 0.001f, 0.0008f);
glm::vec3 colorScale = glm::vec3(origColorScale.x, origColorScale.y, -origColorScale.z);
bool pauseColors = false;

bool changeRed = true;
bool changeGreen = false;
bool changeBlue = true;

float bunnyStrength = 5.0f;

// int loopCount = 0.0;
GLfloat deltaTime = 0.0f;
GLfloat lastFrame = 0.0f;

int main()
{
    createBunnyVertices();
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

    float rCol = 0.0f;
    float gCol = 0.0f;
    float bCol = 0.0f;

    // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

    
    // vel_step(SIZE, velocityX, velocityY, velocityPrevX, velocityPrevY, visc, dt);
    // dens_step(SIZE, density, densityPrev, velocityX, velocityY, diff, dt);            
    // std::cout << loopCount << std::endl;

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
            // printf("\n");
            // printf("frames per sec: %f\n", double(nbFrames));
            glfwSetWindowTitle(window,
                (std::string("Grid Size: ") 
                    + std::string(std::to_string(SIZE))
                    + std::string("x")
                    + std::string(std::to_string(SIZE))
                    + std::string(" , Frame Rate: ")
                    + std::to_string(1.0f / deltaTime).substr(0, 4) + std::string(" fps)")).c_str());
            nbFrames = 0;
            lastTime += 1.0;
            // printf("sec to render frame: %f\n", double(deltaTime));
        }

        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glfwPollEvents();
        processInput();

        processColors();
        GLint scalePos = glGetUniformLocation(shader->getProg(), "scale");
        glUniform3f(scalePos, densityColor.x, densityColor.y, densityColor.z);  

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

void processColors() {
        if (!pauseColors) {
            if (changeRed) {
                densityColor.x += colorScale.x;
            }
            if (changeGreen) {
                densityColor.y += colorScale.y;
            }
            if (changeBlue) {
                densityColor.z += colorScale.z;
            }
        }

        // switches color dir if they hit 1 or 0
        if (densityColor.x >= 0.2f) {
            colorScale.x = -origColorScale.x;
        }
        if (densityColor.x <= 0.01f) {
            colorScale.x = origColorScale.x;
        }

        if (densityColor.y >= 1.0f) {
            colorScale.y = -origColorScale.y;
        }
        if (densityColor.y <= 0.01f) {
            colorScale.y = origColorScale.y;
        }

        if (densityColor.z >= 1.0f) {
            colorScale.z = -origColorScale.z;
        }
        if (densityColor.z <= 0.0f) {
            colorScale.z = origColorScale.z;
        }

        // printVec(densityColor, "densityColor");

}

void createGrid() {
    float cellSize = 2.0f / SIZE;
    float left = -1.0f, top = 1.0f;
    float x = left, y = top;
    
    Rainbow rainbow(0, SIZE);

    // create vertices
    for (int j = 0; j <= SIZE; j++) {
        x = left;
        for (int i = 0; i <= SIZE; i++) {
            vertices.push_back(glm::vec4(x, y, 1.0f, j*1.0f/SIZE));
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

void addDensityToCell(glm::vec2 cell, float densityAmount) {
    int index = at(cell.x, cell.y);
    if (inArrBounds(cell)) {
        data[index].density = densityAmount;
    }
}

// void addSquareDensity(glm::vec2 mousePos) {
//     printVec(mousePos, "mousePos");
//     glm::vec2 cell = getCell(mousePos);
//     int x = cell.x;
//     int y = cell.y;

//     int offset = 3;
//     for (int i = -offset; i < offset; i++) {
//         for (int j = -offset; j < offset; j++) {
//             glm::vec2 newCell(x+i,j+y);
//             addDensityToCell(newCell, densityStrength);
//         }    
//     }
//     debugMiniGrid(cell, offset, "addSquareDensity");
// }

void addSquareDensity(glm::vec2 mousePos) {
    glm::vec2 cell = getCell(mousePos);

    int x = mousePos.x;
    int y = mousePos.y;

    int offset = 10;
    bool foundAmt = false;
    float amt = 100;
    for (int i = -offset; i < offset; i++) {
        for (int j = -offset; j < offset; j++) {
            glm::vec2 newMousePos(x+i,j+y);
            glm::vec2 newCell = getCell(newMousePos);

            float dist = glm::sqrt((newMousePos.x-x)*(newMousePos.x-x) + (newMousePos.y-y)*(newMousePos.y-y));

            if (inArrBounds(newMousePos) && !foundAmt) {
                amt = dist;
                foundAmt = true;
            }
            // float dens = glm::sin(dist) * 100;
            // if (dist < 200 && dens > -100) {
            if (dist < 8) {
                addDensityToCell(newCell, amt-dist);
            }
        }    
    }
    // debugMiniGrid(cell, offset, "addDensity");
}

void addSineDensity(glm::vec2 mousePos) {
    glm::vec2 cell = getCell(mousePos);

    int x = mousePos.x;
    int y = mousePos.y;

    int offset = 300;
    for (int i = -offset; i < offset; i++) {
        for (int j = -offset; j < offset; j++) {
            glm::vec2 newMousePos(x+i,j+y);
            glm::vec2 newCell = getCell(newMousePos);

            float dist = glm::sqrt((newMousePos.x-x)*(newMousePos.x-x) + (newMousePos.y-y)*(newMousePos.y-y));

            float dens = glm::sin(dist) * 100;
            if (dist < 200) {
                addDensityToCell(newCell, dist);
            }
        }    
    }
    debugMiniGrid(cell, offset, "addDensity");
}

void deleteDensityFromCell(glm::vec2 mousePos) {
    glm::vec2 cell = getCell(mousePos);

    int x = cell.x;
    int y = cell.y;

    int offset = 4;
    for (int i = -offset; i < offset; i++) {
        for (int j = -offset; j < offset; j++) {
            glm::vec2 newCell(x+i,j+y);
            addDensityToCell(cell, -densityStrength * 10);
        }    
    }
    // debugMiniGrid(cell, offset, "deleteDensityFromCell");
}

void addVelocityToCell(glm::vec2 newMousePos, glm::vec2 lastMousePos) {
    glm::vec2 cell = getCell(newMousePos);
    int index = at(cell.x, cell.y);
    
    glm::vec2 offset = newMousePos - lastMousePos;
    if (inArrBounds(cell) && inScreenBounds(newMousePos)) {
        data[index].velocity = force * offset;  
    }
}

void addPointForce(glm::vec2 mousePos) {
    glm::vec2 cell = getCell(mousePos);
    int x = cell.x;
    int y = cell.y;
    int offset = 10;
    float rad = 5;
    if (inScreenBounds(mousePos + glm::vec2(offset)) && inScreenBounds(mousePos - glm::vec2(offset))) {   
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

void addLineDensity(glm::vec2 firstClick, glm::vec2 currMousePos) {
    int minX = glm::min(firstClick.x, currMousePos.x);
    int maxX = glm::max(firstClick.x, currMousePos.x);

    int minY = glm::min(firstClick.y, currMousePos.y);
    int maxY = glm::max(firstClick.y, currMousePos.y);
    for (int i = minX; i < maxX; i++) {
        for (int j = minY; j < maxY; j++) {
            glm::vec2 cell = getCell(glm::vec2(i, j));
            addDensityToCell(cell, 50.0f);
        }
    }
}

void spawnBunnies() {
    for (int i = 0; i < 500; i ++) {
        drawBunny(glm::vec2(randFloat(0, 700), randFloat(0, 700)), randFloat(2, 5), bunnyStrength);
    }
}

void spawnUglyBunnies() {
    for (int i = 0; i < 200; i ++) {
        drawBunny(glm::vec2(randFloat(0, 700), randFloat(0, 700)), randFloat(2, 5), -bunnyStrength);
    }
}

void spawnGiantBunny(glm::vec2 mousePos) {
    drawBunny(mousePos, 30, bunnyStrength);
}


bool keys[1024]; 
bool firstMouse = true;
bool mouseClicked = false;
bool mouseReleased = false;
bool mouseDown = false;
glm::vec2 currMousePos, oldMousePos;
glm::vec2 clickedPos;

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
    // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
} 

// called when ever the mouse is clicked
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    if (state == GLFW_PRESS) {
        mouseClicked = true;
        mouseDown = true;
        clickedPos = currMousePos;
        // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }
    else if (state == GLFW_RELEASE) {
        mouseDown = false;
        mouseReleased = true;
        addLineDensity(clickedPos, currMousePos);
    }
    // if (mouseDown) {
    //     glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    // }
}

void scaleVal(int key, float *ptr, float scale, std::string varName) {
    if (keys[key]) {
        if(keys[GLFW_KEY_LEFT_SHIFT]) {
            *ptr -= scale;
        }
        else {
            *ptr += scale;
        }
        std::cout << varName << ": " << *ptr << std::endl;
    }
}
void toggleVal(int key, bool *ptr, std::string varName) {
    if (keys[key]) {
        if(keys[GLFW_KEY_LEFT_SHIFT]) {
            *ptr = true;
        }
        else {
            *ptr = false;
        }
        std::cout << varName << ": " << *ptr << std::endl;
    }
}
void toggleSim(int key, bool *ptr, std::string varName) {
    if (keys[key]) {
        if(keys[GLFW_KEY_LEFT_SHIFT]) {
            *ptr = false;
        }
        else {
            *ptr = true;
        }
        std::cout << varName << ": " << *ptr << std::endl;
    }
}

void processInput() {
    if (inScreenBounds(currMousePos) && inScreenBounds(oldMousePos)) {
        // adds density
        if (keys[GLFW_KEY_F]) {
            addSquareDensity(currMousePos);
        }
        // adds velocity
        if (keys[GLFW_KEY_V]) {
            addVelocityToCell(currMousePos, oldMousePos);
        }
        // eraser
        if (keys[GLFW_KEY_D]) {
            deleteDensityFromCell(currMousePos);
        }

        // draws pos bunny
        if (keys[GLFW_KEY_C]) {
            drawBunny(getCell(currMousePos+40.0f), 5, bunnyStrength);
        }
        // draws neg bunny
        if (keys[GLFW_KEY_X]) {
            drawBunny(getCell(currMousePos+40.0f), 5, -bunnyStrength);
        }

        // spawns multiple pos bunny
        if (keys[GLFW_KEY_L]) {
            spawnBunnies();
        }
        // spawns multiple neg bunny
        if (keys[GLFW_KEY_K]) {
            spawnUglyBunnies();
        }
        // draws large pos bunny
        if (keys[GLFW_KEY_H]) {
            spawnGiantBunny(getCell(currMousePos + 70.0f));
        }
        // draws large neg bunny
        if (keys[GLFW_KEY_G]) {
            drawBunny(getCell(currMousePos + 70.0f), 14, -bunnyStrength);
        }

        // pause sim
        toggleSim(GLFW_KEY_O, &pauseSim, "pauseSim");
        // pause step
        toggleVal(GLFW_KEY_I, &step, "step");

        // add's point force
        if (keys[GLFW_KEY_P]) {
            addPointForce(currMousePos);
        }

        // scales force
        scaleVal(GLFW_KEY_A, &force, 1.0f, "force");
        scaleVal(GLFW_KEY_Z, &dt, 0.01f, "dt");

        //changes visc
        if (keys[GLFW_KEY_B]) {
            if (keys[GLFW_KEY_LEFT_SHIFT]) {
                visc = 0.0f;
            }
            else {
                visc += 0.0001f;
            }
            std::cout << "visc: " << visc << std::endl;
        }
        
        // changes rbg
        scaleVal(GLFW_KEY_R, &densityColor.x, 0.001f, "densityColor.x");
        scaleVal(GLFW_KEY_T, &densityColor.y, 0.001f, "densityColor.y");
        scaleVal(GLFW_KEY_Y, &densityColor.z, 0.001f, "densityColor.z");
      
        // pauses colors
        toggleVal(GLFW_KEY_S, &pauseColors, "pauseColors");


        if (keys[GLFW_KEY_1]) {
            changeRed = !changeRed;
        }
        if (keys[GLFW_KEY_2]) {
            changeGreen = !changeGreen;
        }
        if (keys[GLFW_KEY_3]) {
            changeBlue = !changeBlue;
        }

        if (keys[GLFW_KEY_Q]) {
            if (keys[GLFW_KEY_LEFT_SHIFT]) {
                colorScale.x = -origColorScale.x;
            }
            else {
                colorScale.x = origColorScale.x;   
            }
        }
        if (keys[GLFW_KEY_W]) {
            if (keys[GLFW_KEY_LEFT_SHIFT]) {
                colorScale.y = -origColorScale.y;
            }
            else {
                colorScale.y = origColorScale.y;   
            }
        }
        if (keys[GLFW_KEY_E]) {
            if (keys[GLFW_KEY_LEFT_SHIFT]) {
                colorScale.z = -origColorScale.z;
            }
            else {
                colorScale.z = origColorScale.z;   
            }
        }
    }
}

bool inScreenBounds(glm::vec2 mousePos) {
    int x = mousePos.x;
    int y = mousePos.y;
    return x <= WIDTH && x >=0 && y <= HEIGHT && y>=0;
}

bool inArrBounds(glm::vec2 index) {
    int x = index.x;
    int y = index.y;
    return  x <= SIZE && x >=1 && y <= SIZE && y>=1; 
}

void debugMiniGrid(glm::vec2 cell, int offset, std::string name) {
    std::cout << name << ":" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    int x = cell.x;
    int y = cell.y;
    for (int i = -offset; i < offset; i++) {
        for (int j = -offset; j < offset; j++) {
            if (inArrBounds(cell)) {
                int index = at(cell.x, cell.y);
                if (x+i == x && y+j == y) {
                    std::cout << "~" << data[index].density << "~ ";
                }
                else {
                    std::cout << data[index].density << " ";
                }
            }
        }
        std::cout << std::endl;
    }
}

void createBunnyVertices() {
    const char *path = "../resources/bunny.obj";
    // const char *path = "/Users/ashleydattalo/graphics/471/THE_FINAL_AUG_11_2016/finalPro/finalProject/resources/voyager.obj";
    // const char *path = "/Users/ashleydattalo/graphics/471/THE_FINAL_AUG_11_2016/finalPro/finalProject/resources/superman.obj";
    // const char *path = "/Users/ashleydattalo/graphics/471/THE_FINAL_AUG_11_2016/finalPro/finalProject/resources/NewTieFighter.obj";
    // const char *path = "/Users/ashleydattalo/graphics/471/THE_FINAL_AUG_11_2016/finalPro/finalProject/resources/UFO.obj";
    // const char *path = "/Users/ashleydattalo/graphics/471/THE_FINAL_AUG_11_2016/finalPro/finalProject/resources/newE-wing.obj";
    // const char *path = "/Users/ashleydattalo/graphics/471/THE_FINAL_AUG_11_2016/finalPro/finalProject/resources/spaceShip.obj";
    // const char *path = "/Users/ashleydattalo/graphics/471/THE_FINAL_AUG_11_2016/finalPro/finalProject/resources/sphere.obj";
    // const char *path = "/Users/ashleydattalo/graphics/471/THE_FINAL_AUG_11_2016/finalPro/finalProject/resources/square.obj";

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
            bunnyVerts.push_back(vertex);
        }
        if ( strcmp( lineHeader, "vn" ) == 0 ){
            glm::vec3 normal;
            fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z );
        }
    }
}

void drawBunny(glm::vec2 offset, float bunnySize, float densStrength) {
    for (glm::vec2 vertex : bunnyVerts) {
        vertex.y *= -1;
        vertex *= 100 *bunnySize;
        glm::vec2 idx = getCell(glm::vec2(vertex.x, vertex.y));
        idx += offset;
        addDensityToCell(idx, densStrength);

        // for (int i = -2; i < 2; i++) {
            // for (int j = -2; j < 2; j++) {
                // glm::vec2 newIndex = glm::vec2(vertex.x+i, vertex.y+j);
                glm::vec2 newIndex = glm::vec2(vertex.x, vertex.y);
            // }    
        // }
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

    window = glfwCreateWindow(WIDTH, HEIGHT, "Final Project", nullptr, nullptr);    
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

// Solver
void set_bnd ( int N, int b, float * x ) {
    int i;
    for ( i=1 ; i<=N ; i++ ) {
        x[IX(0 ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)]; 
        x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
        x[IX(i,0 )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
        x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)]; 
    }
    x[IX(0 ,0 )] = 0.5*(x[IX(1,0 )]+x[IX(0 ,1)]); 
    x[IX(0 ,N+1)] = 0.5*(x[IX(1,N+1)]+x[IX(0 ,N )]); 
    x[IX(N+1,0 )] = 0.5*(x[IX(N,0 )]+x[IX(N+1,1)]); 
    x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)]+x[IX(N+1,N)]);
}

void add_source ( int N, float * x, float * s, float dt ) {
    int i, size = (N+2)*(N+2);
    for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i]; 
}

void diffuse ( int N, int b, float * x, float * x0, float diff, float dt ) {
    int i, j, k;
    float a=dt*diff*N*N;
    for ( k=0 ; k<80 ; k++ ) {
        for ( i=1 ; i<=N ; i++ ) {
            for ( j=1 ; j<=N ; j++ ) {
                x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/(1+50*a);
            }
        }
        set_bnd ( N, b, x );
    }
}

void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt ) {
    int i, j, i0, j0, i1, j1;
    float x, y, s0, t0, s1, t1, dt0;
    dt0 = dt*N;
    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=N ; j++ ) {
            x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];
            if (x<0.5) x=0.5; 
            if (x>N+0.5) x=N+0.5; 
            i0=(int)x;
            i1=i0+1; 
            if (y<0.5) y=0.5; 
            if (y>N+0.5) y=N+0.5; 
            j0=(int)y; j1=j0+1; 
            s1 = x-i0; s0 = 1-s1; 
            t1 = y-j0; t0 = 1-t1;
            d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
        } 
    }
    set_bnd ( N, b, d ); 
}

void project ( int N, float * u, float * v, float * p, float * div ) {
    int i, j, k;
    float h;
    h = 1.0/N;
    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=N ; j++ ) {
            div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)]); 
            p[IX(i,j)] = 0;
        }
    }
    set_bnd ( N, 0, div ); 
    set_bnd ( N, 0, p );

    for ( k=0 ; k<80 ; k++ ) {
        for ( i=1 ; i<=N ; i++ ) {
            for ( j=1 ; j<=N ; j++ ) {
                p[IX(i,j)] = (div[IX(i,j)]+p[IX(i-1,j)]+p[IX(i+1,j)]+p[IX(i,j-1)]+p[IX(i,j+1)])/4 ;
            } 
        }
        set_bnd ( N, 0, p ); 
    }

    for ( i=1 ; i<=N ; i++ ) {
        for ( j=1 ; j<=N ; j++ ) {
            u[IX(i,j)] -= 0.5*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
            v[IX(i,j)] -= 0.5*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
        }
    }
    set_bnd ( N, 1, u ); set_bnd ( N, 2, v ); 
}

void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt ) {
    add_source ( N, x, x0, dt );
    SWAP ( x0, x ); diffuse ( N, 0, x, x0, diff, dt ); 
    SWAP ( x0, x ); advect ( N, 0, x, x0, u, v, dt );
}


void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt ) {
    add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
    SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt );
    SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt );
    project ( N, u, v, u0, v0 );
    SWAP ( u0, u ); SWAP ( v0, v );
    advect ( N, 1, u, u0, u0, v0, dt );
    advect ( N, 2, v, v0, u0, v0, dt );
    project ( N, u, v, u0, v0 );
}