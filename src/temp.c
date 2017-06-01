#include <iostream>
#include <stdlib.h>
#include <glm/glm.hpp>

#define SIZE 10
#define START 1
#define END SIZE + 1
#define at(i,j) ((i) + (SIZE + 2)*j)

void vel_step(glm::vec3 *vel, glm::vec3 *velPrev, float visc, float dt);
void dens_step();

void add_source(int gridSize, float *result, float *source, float dt);
void diffuse(int gridSize, int b, float *x, float *x0, float diff, float dt);
void advect(int gridSize, int b, float *d, float *d0, float *u, float *v, float dt);

void printArr(float *arr, std::string arrName);
void printVel(glm::vec3 *arr, std::string arrName);

void initDensity(float *density);
void initVelocity(glm::vec3 *velocity);

int main() {

	float *density = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
	glm::vec3 *velocity = (glm::vec3 *) malloc(sizeof(glm::vec3) * (SIZE + 2) * (SIZE + 2));

	float visc = 0.5f;
	float dt = 0.1f;

	// game loop
	for (int i = 0; i < 5; i++) {
		std::cout << "i: " << i << std::endl;

		// dens_step();
		// printArr(density, "density");
		printVel(velocity, "velocity");
		std::cout << std::endl << std::endl << std::endl;
	}
}

void vel_step(glm::vec3 *vel, glm::vec3 *velPrev, float visc, float dt) {
	addSourceVel(vel, velPrev, dt);
	swapVel(velPrev, vel);
	diffuseVel(vel, velPrev, visc, dt);
	projectVel(vel, velPrev);
	advectVel(vel, velPrev, dt);
	projectVel(vel, velPrev);
}

void dens_step() {

}

void add_source(int gridSize, float *result, float *source, float dt) {
	int size = (gridSize + 2) * (gridSize + 2);
	for (int i = 0; i < size; i++) {
		result[i] += dt * source[i];
	}
}

void diffuse(float *x, float *x0, float diff, float dt) {
	float a = dt * diff * SIZE * SIZE;

	for (int k = 0; k < 20; k++) {
		for (int i = START; i < END; i++) {
			for (int j = START; j < END; j++) {
				float surrDensity = x[at(i-1,j)] + x[at(i+1, j)] + x[at(i,j-1)] + x[at(i,j+1)];
				x[at(i,j)] = (x0[at(i,k)] + a * surrDensity) / (1 + 4 * a);
			}
		}
		// setBoundaries(gridSize, b, x);
	}
}

void diffuse(glm::vec3 *vel, glm::vec3 *velPrev, float diff, float dt) {
	float a = dt * diff * SIZE * SIZE;

	for (int k = 0; k < 20; k++) {
		for (int i = START; i < END; i++) {
			for (int j = START; j < END; j++) {
				float surrDensity = vel[at(i-1,j)] + vel[at(i+1, j)] + vel[at(i,j-1)] + vel[at(i,j+1)];
				vel[at(i,j)] = (velPrev[at(i,k)] + a * surrDensity) / (1 + 4 * a);
			}
		}
		// setBoundaries(gridSize, b, x);
	}
}

void advect(int gridSize, int b, float *d, float *d0, float *u, float *v, float dt) {
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt * gridSize;

	for (i = 1; i <= gridSize; i++) {
		for (j = 1; j <= gridSize; j++) {
			x = i - dt * u[at(i,j)];
			y = j - dt * v[at(i,j)];

			//set x vals
			if (x < 0.5) {
				x = 0.5;
			}
			if (x > gridSize + 0.5) {
				x = gridSize + 0.5;
			}
			i0 = (int) x;
			i1 = i0 + 1;

			s1 = x - i0;
			s0 = 1 - s1;

			//set y vals
			if (y < 0.5) {
				y = 0.5;
			}
			if (y > gridSize + 0.5) {
				y = gridSize + 0.5;
			}
			j0 = (int) y;
			j1 = j0 + 1;

			t1 = y - j0;
			t0 = 1 - t1;

			d[at(i,j)] = s0 * (t0 * d0[at(i0, j0)] + t1 * d0[at(i0,j1)]) + 
						 s1 * (t0 * d0[at(i1, j0)] + t1 * d0[at(i1,j1)]);
		}
	}
	// setBoundaries(gridSize, b, d);
}






void initDensity(float *density) {
	for (int i = START; i < END; i++) {
		for (int j = START; j < END; j++) {
			float *num = (float *) malloc(sizeof(float));
			*num = 1.0f;
			density[at(i, j)] = j;
		}
	}
}	

void initVelocity(glm::vec3 *velocity) {
	for (int i = START; i < END; i++) {
		for (int j = START; j < END; j++) {
			velocity[at(i, j)] = glm::vec3(i, j, 0.0f);
		}
	}
}

void printArr(float *arr, std::string arrName) {
	std::cout << arrName << ": " << std::endl;
	for (int j = START; j < END; j++) {
		std::cout << "j: ";
		for (int i = START; i < END; i++) {
			std::cout << arr[at(i, j)] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void printVel(glm::vec3 *arr, std::string arrName) {
	std::cout << arrName << ": " << std::endl;
	for (int j = START; j < END; j++) {
		std::cout << "j: ";
		for (int i = START; i < END; i++) {
			glm::vec3 val = arr[at(i, j)];
			std::cout << "(";
			std::cout << val.x;
			std::cout << ",";
			std::cout << val.y;
			std::cout << ",";
			std::cout << val.z;
			std::cout << ") ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}