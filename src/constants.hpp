#ifndef Constants_HPP
#define Constants_HPP

#define SHADER_PATH "/Users/ashleydattalo/graphics/gradGraphics/finalProject/openGL/src/Shaders/"
#define RESOURCE_PATH "/Users/ashleydattalo/graphics/LearnOpenGL/resources/"

float randFloat(float min, float max) {
    float r = rand() / (float)RAND_MAX;
    return (1.0f - r) * min + r * max;
}

void printVec(glm::vec3 toPrint) {
    std::cout << toPrint.x << " ";
    std::cout << toPrint.y << " ";
    std::cout << toPrint.z << " ";
    std::cout << std::endl;
}

#endif