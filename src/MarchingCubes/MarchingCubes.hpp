#ifndef MatchingCubes_HPP
#define MatchingCubes_HPP

#define NUM_EDGES 12

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <glm/glm.hpp>

class MarchingCubes {
public:
	MarchingCubes();
	void init();

	void setVertList(int cubeIndex);
	glm::vec3 VertexInterp(glm::vec3 p1, glm::vec3 p2, float v1, float v2);
	int getCubeIndex();
	void createGrid(float x, float y, float z);
	float implicitFunction(glm::vec3 point);
	float implicitFunction(float x, float y, float z);

	glm::vec3 gradInterp(int vert1, int vert2);
	glm::vec3 getGradient(glm::vec3 point);
	
	std::vector<glm::vec3> getVertices();
	std::vector<glm::vec3> getNormals();

	int getNumVertices();

	glm::vec3 getCenter();
	glm::vec3 getGridSize();

private:
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> normals;
	float X_SIZE;
	float Y_SIZE;
	float Z_SIZE;
	
	float cx, cy, cz;
	float inc;
	float isolevel;
	float radius;

	int numVertices;
	
	std::vector<glm::vec3> vertList;
	std::vector<glm::vec3> gradList;
	glm::vec3 currPoint;
};
  
#endif