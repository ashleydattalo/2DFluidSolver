#include "ParticleSystem.hpp"
#include "Particle.hpp"
#include "../Classes/Shader.hpp"

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/string_cast.hpp>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <vector>
#include <tgmath.h>

#define SHADER_PATH "/Users/ashleydattalo/graphics/LearnOpenGL/src/Shaders/"

#define HEIGHT 100

static int numPSys = 0;
float currTime = 0;
ParticleSystem::ParticleSystem()
{
    origin = glm::vec3(0.0f, 0.0f, 0.0f);
    particleType = "general";
    nParticles = 15000;
    id = numPSys++;
    init();
}
    
ParticleSystem::ParticleSystem(glm::vec3 origin, std::string type, int numParticles) :
    origin(origin), particleType(type), nParticles(numParticles) 
    {
        id = numPSys++;
        init();
    }

ParticleSystem::~ParticleSystem() {}

void ParticleSystem::init() {
    posBuf.resize(3*nParticles);
    colBuf.resize(3*nParticles);
    alpBuf.resize(nParticles);
    scaBuf.resize(nParticles);

    radBuf.resize(nParticles);
    initTimeBuf.resize(nParticles);

    createParticles();
    createGLBuffers();
    createShader();
}
float ParticleSystem::getR(float percent) {
    if (percent <= .5) {
        return 1 - (2 * percent);
    }
    else {
        return (percent-.5) * 2;
    }
}
float ParticleSystem::getG(float percent) {
    if (percent <= .25) {
        return percent * 4;
    }
    else if (percent <= .5) {
        return 1;
    }
    else if (percent <= .75) {
        return 1 - ((percent-.5) * 4);
    } 
    else {
        return 0;
    }
}
float ParticleSystem::getB(float percent) {
    if (percent <= .25) {
        return 0;
    }
    else if (percent <= .5) {
        return (percent-.25) * 4;
    }
    else if (percent <= .75) {
        return 1;
    } 
    else {
        return 1-((percent-.75) * 4);
    }
}
glm::vec3 ParticleSystem::rgbVal(float percent) {
    return glm::vec3(getR(percent), getG(percent), getB(percent));
}

float ParticleSystem::getPercentage(float y) {
    float min = 0;
    float max = HEIGHT;

    return y / (max - min);
}

void ParticleSystem::setRainbow() {
    // float min = 0;
    // float max = HEIGHT;

    // float step = 1 / (max - min);
    // //ex: height = 4
    // //    .25
    // //    .5
    // //    .75
    // //    1

    // static float percentage = 0;

    // for (int i = 1; i <= max; i++) {
    //     //step = step + percentage;
    // }
    // std::cout << std::endl;
    for (Particle *p : particles) {
        glm::vec3 pos = p->getPosition();
        float percentage = getPercentage(pos.y);
        // std::cout << currTime << " ";
        // std::cout << percentage << " ";
        // std::cout << std::endl;
        if (percentage > 1 || percentage < 0) {
            // std::cout << "percent is outside range [0,1]" << std::endl;
            exit(0);
        }
        // std::cout << percentage<<": "<< glm::to_string(rgbVal(percentage)) << std::endl;
        p->setColor(rgbVal(percentage));
    }
}

void ParticleSystem::update(float dt) {
    // setRainbow();
    for (Particle *p : particles) {
        p->update(dt);
        // updateGLData(p);
    }
}

void ParticleSystem::render(const glm::mat4 &projection, const glm::mat4 &view) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // glDisable(GL_DEPTH_TEST);

    shader->use();

    glBindVertexArray(VAO);

    currTime = fmod(glfwGetTime(), HEIGHT);

    std::cout << glfwGetTime()/10 << std::endl;
    GLint timeLoc = glGetUniformLocation(shader->getProg(), "time");
    glUniform1f(timeLoc, glfwGetTime()/10);

    GLint projLoc = glGetUniformLocation(shader->getProg(), "P");
    glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(projection));

    GLint viewLoc = glGetUniformLocation(shader->getProg(), "V");
    glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));

    // //send pos data to gpu
    // glBindBuffer(GL_ARRAY_BUFFER, posBufID);
    // glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_STREAM_DRAW);

    // glBindBuffer(GL_ARRAY_BUFFER, alpBufID);
    // glBufferData(GL_ARRAY_BUFFER, alpBuf.size()*sizeof(float), &alpBuf[0], GL_STREAM_DRAW);

    // //send color data to gpu
    // glBindBuffer(GL_ARRAY_BUFFER, colBufID);
    // glBufferData(GL_ARRAY_BUFFER, colBuf.size()*sizeof(float), &colBuf[0], GL_DYNAMIC_DRAW);
    
    glDrawArrays(GL_POINTS, 0, nParticles);

    glBindVertexArray(0);
}

glm::vec3 genPos() {
    static float radius = 1;
    static float x = -radius;
    

    float y = sqrt(radius * radius - x * x);
    std::cout << y << " " << radius << " " << x << std::endl;
    glm::vec3 pos;
    pos.x = x;
    pos.y = -y;

    // radius++;
    x += .1;
    return pos;
}

void ParticleSystem::createParticles() {
    static float color = 0;
    float inc = 1.0f / nParticles;
    for (int i = 1; i <= nParticles; i++) {
        Particle *p = new Particle(i, origin, particleType);
        color += inc;
        // glm::vec3 col = glm::vec3(.3, .2, color);
        glm::vec3 col = rgbVal(color);
        // glm::vec3 pos = glm::vec3(randNum(-1.0f, 1.0f), randNum(-355.0f,-205.0f), randNum(0.5f, 1.0f));
        glm::vec3 pos = glm::vec3(randNum(-5.0f, 5.0f), HEIGHT*color, randNum(-5.0f, 5.0f));
        float scale = 50 * color;
        float t = 20*color;

        float radius = color;
        // pos.x = radius * sin(t);
        // pos.y += radius * (t);
        // pos.z = radius * cos(t);


        radBuf[i] = radius;
        initTimeBuf[i] = t;

        float alpha = 1;

        // col = glm::vec3(1.0f);
        scale = 5;
        // alpha = 1;

        p->setPosition(pos);
        p->setColor(col);
        p->setScale(scale);
        p->setAlpha(alpha);
        p->setT(t);

        particles.push_back(p);
        updateGLData(p);
    }
}

void ParticleSystem::createShader() {
    shader = new Shader(SHADER_PATH "particles.vert", SHADER_PATH "general.frag");
}

void ParticleSystem::createGLBuffers() {
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);
    // glEnable(GL_POINT_SMOOTH);
    GLuint bufs[6];
    glGenBuffers(6, bufs);
    std::cout << "createGLBuffers" << std::endl;
    posBufID = bufs[0];
    colBufID = bufs[1];
    alpBufID = bufs[2];
    scaBufID = bufs[3];
    radBufID = bufs[4];
    initTimeBufID = bufs[5];

    glGenVertexArrays(1, &VAO);
    // Bind the Vertex Array Object first, then bind and set vertex buffer(s) and attribute pointer(s).
    glBindVertexArray(VAO);

    // Send pos data to the GPU
    glBindBuffer(GL_ARRAY_BUFFER, posBufID);
    glBufferData(GL_ARRAY_BUFFER, posBuf.size()*sizeof(float), &posBuf[0], GL_STREAM_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);


    // Send color buffer to GPU
    glBindBuffer(GL_ARRAY_BUFFER, colBufID);
    glBufferData(GL_ARRAY_BUFFER, colBuf.size()*sizeof(float), &colBuf[0], GL_STREAM_DRAW);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (GLvoid*)0);
    
    // Send scale buffer to GPU
    glBindBuffer(GL_ARRAY_BUFFER, scaBufID);
    glBufferData(GL_ARRAY_BUFFER, scaBuf.size()*sizeof(float), &scaBuf[0], GL_STREAM_DRAW);
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(GLfloat), (GLvoid*)0);

    // Send alpha buffer to GPU
    glBindBuffer(GL_ARRAY_BUFFER, alpBufID);
    glBufferData(GL_ARRAY_BUFFER, alpBuf.size()*sizeof(float), &alpBuf[0], GL_STREAM_DRAW);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(GLfloat), (GLvoid*)0);

    //radius
    glBindBuffer(GL_ARRAY_BUFFER, radBufID);
    glBufferData(GL_ARRAY_BUFFER, radBuf.size()*sizeof(float), &radBuf[0], GL_STREAM_DRAW);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(GLfloat), (GLvoid*)0);

    //initial time
    glBindBuffer(GL_ARRAY_BUFFER, initTimeBufID);
    glBufferData(GL_ARRAY_BUFFER, initTimeBuf.size()*sizeof(float), &initTimeBuf[0], GL_STREAM_DRAW);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, sizeof(GLfloat), (GLvoid*)0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    std::cout << "endBuf" << std::endl;
}

void ParticleSystem::updateGLData(Particle *p) {
    int index = p->getIndex();
    glm::vec3 col = p->getColor();
    glm::vec3 pos = p->getPosition();

    posBuf[3*index + 0] = pos.x;
    posBuf[3*index + 1] = pos.y;
    posBuf[3*index + 2] = pos.z;

    colBuf[3*index + 0] = col.x;
    colBuf[3*index + 1] = col.y;
    colBuf[3*index + 2] = col.z;

    scaBuf[index] = p->getScale();
    alpBuf[index] = p->getAlpha();

    glBindBuffer(GL_ARRAY_BUFFER, colBufID);
    glBufferSubData(GL_ARRAY_BUFFER, 3*index*sizeof(float), 3*sizeof(float), &colBuf[3*index]);

    glBindBuffer(GL_ARRAY_BUFFER, posBufID);
    glBufferSubData(GL_ARRAY_BUFFER, 3*index*sizeof(float), 3*sizeof(float), &posBuf[3*index]);

    glBindBuffer(GL_ARRAY_BUFFER, scaBufID);
    glBufferSubData(GL_ARRAY_BUFFER, index*sizeof(float), sizeof(float), &scaBuf[index]);

    glBindBuffer(GL_ARRAY_BUFFER, alpBufID);
    glBufferSubData(GL_ARRAY_BUFFER, index*sizeof(float), sizeof(float), &alpBuf[index]);

    glBindBuffer(GL_ARRAY_BUFFER, radBufID);
    glBufferSubData(GL_ARRAY_BUFFER, index*sizeof(float), sizeof(float), &radBuf[index]);

    glBindBuffer(GL_ARRAY_BUFFER, initTimeBufID);
    glBufferSubData(GL_ARRAY_BUFFER, index*sizeof(float), sizeof(float), &initTimeBuf[index]);
}

float ParticleSystem::randNum(float min, float max) {
    float r = rand() / (float)RAND_MAX;
    return (1.0f - r) * min + r * max;
}

void ParticleSystem::setOrigin(glm::vec3 newPos) {

}

void ParticleSystem::debug() {
    for (Particle *p : particles) {
        std::cout << p->getScale();
        std::cout << "   ";
        std::cout << p->getAlpha();
        std::cout << std::endl;
    }
}
