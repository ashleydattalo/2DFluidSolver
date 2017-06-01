#include "Rainbow.hpp"


Rainbow::Rainbow(float min, float max) {
    this->min = min;
    this->max = max;
}

Rainbow::~Rainbow() {}

glm::vec3 Rainbow::getColor(float value) {
    float percent = value / (max - min);
    return glm::vec3(getR(percent), getG(percent), getB(percent));
}

glm::vec3 Rainbow::getColor(glm::vec3 value) {
    glm::vec3 percent = value / (max - min);
    return glm::vec3(getR(percent.x), getG(percent.y), getB(percent.z));
}

float Rainbow::getR(float percent) {
    if (percent <= .5) {
        return 1 - (2 * percent);
    }
    else {
        return (percent-.5) * 2;
    }
}
float Rainbow::getG(float percent) {
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
float Rainbow::getB(float percent) {
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