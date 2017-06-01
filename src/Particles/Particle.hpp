#pragma once
#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include <glm/glm.hpp>
#include <string>

class Particle
{
public:
	
	Particle();
	Particle(int index, glm::vec3 origin, std::string type);
	Particle(glm::vec3 position, glm::vec3 color, float scale, float alpha);
	virtual ~Particle();
	void init();
	void rebirth(float dt);
	void update(float dt);
	
	float randFloat(float min, float max);

	void setPosition(glm::vec3 pos);
	void setColor(glm::vec3 col);
	void setScale(float sca);
	void setAlpha(float alp);
	void setT(float t);

	void updatePosition(glm::vec3 update);
	int getIndex();
	glm::vec3 getPosition();
	glm::vec3 getColor();
	float getScale();
	float getAlpha();
	void setOrigin(glm::vec3 origin);

private:
	int index;
	glm::vec3 origin;
	std::string type;

	glm::vec3 color;
	glm::vec3 position;
	float scale;
	float alpha;
	float t;

	glm::vec3 velocity;
	glm::vec3 acc;
	glm::vec3 force;

	float mass;
	float lifespan;
	float tEnd;

	bool alphaToggle;
};

#endif
