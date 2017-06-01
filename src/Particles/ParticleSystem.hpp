#include "Particle.hpp"
#include "../Classes/Shader.hpp"
#include <glm/glm.hpp>

#include <GL/glew.h>
#include <GLFW/glfw3.h>


#include <iostream>
#include <vector>

class ParticleSystem {
public:
    ParticleSystem();
    ParticleSystem(glm::vec3 origin, std::string type, int numParticles);
    ~ParticleSystem();

    void init();
    void update(float dt);
    void render(const glm::mat4 &projection, const glm::mat4 &view);
    void createParticles();
    void createShader();
    void createGLBuffers();
    void updateGLData(Particle *p);
    float randNum(float min, float max);
    void setOrigin(glm::vec3 newPos);
    // bool inFrustrum();
    // void collidedWithOrb(const glm::vec3 &playerPos);
    // void setForce(const glm::vec3 &force);
    void debug();
    float getR(float percent);
    float getG(float percent);
    float getB(float percent);
    glm::vec3 rgbVal(float percent);
    float getPercentage(float y);
    void setRainbow();

private:
	int id;
    int nParticles;
    std::vector<Particle *> particles;
    std::string particleType;
    glm::vec3 origin;

    GLuint VBO, VAO, EBO;

    std::vector<float> posBuf;
    std::vector<float> colBuf;
    std::vector<float> alpBuf;
    std::vector<float> scaBuf;
    std::vector<float> radBuf;
    std::vector<float> initTimeBuf;
    GLuint posBufID;
    GLuint colBufID;
    GLuint alpBufID;
    GLuint scaBufID;
    GLuint radBufID;
    GLuint initTimeBufID;

    Shader *shader;

    // GLShaderProgram particleShaderProg;
    // ViewFrustrum frustrum;
};