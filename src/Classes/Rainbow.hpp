#ifndef Rainbow_HPP
#define Rainbow_HPP

#include <glm/glm.hpp>

class  Rainbow{
public:
	Rainbow(float min, float max);
	~Rainbow();
	glm::vec3 getColor(float value);
	glm::vec3 getColor(glm::vec3 value);
	float getR(float percent);
	float getG(float percent);
	float getB(float percent);

private:
	float max;
	float min;
};
  
#endif