#include <unistd.h>
#include <iostream>
#include <stdlib.h>
#include <glm/glm.hpp>

#include "Image.h"
#include "types.h"

#define TIME 2000000
#define SIZE 300
#define START 1
#define END SIZE + 1
#define at(i,j) ((i) + (SIZE + 2)*j)
#define SWAP(x0,x) {float *tmp=x0;x0=x;x=tmp;}

void vel_step(float *velX, float *velY, float *velPrevX, float *velPrevY, float visc, float dt);
void dens_step(float *density, float *densityPrev, float *velX, float *velY, float diff, float dt);

void add_source(float *result, float *source, float dt);
void diffuse(int b, float *x, float *x0, float diff, float dt);
void advect(int b, float *d, float *d0, float *u, float *v, float dt);
void project(float *u, float *v, float *p, float *divD);

void set_bnd (int b, float *x);


void printArr(float *arr, std::string arrName);
void printVel(float *velX, float *velY);
void printVel(float *velX, float *velY, std::string name);
void printVel(glm::vec3 *arr, std::string arrName);
void initDensity(float *density);
void initVelocity(float *velocity);
void initVelocity(glm::vec3 *velocity);
void drawDensity(float *density);

void init(float *densityPrev, float *velocityPrevX, float *velocityPrevY);
void addDensitySource(float *density);
void addForce(float *velocityX, float *velocityY, float forceX, float forceY);

int main() {

	float *density = (float *) calloc(1.0f, sizeof(float) * (SIZE + 2) * (SIZE + 2));
	float *densityPrev = (float *) calloc(1.0f, sizeof(float) * (SIZE + 2) * (SIZE + 2));

	float *velocityX = (float *) calloc(1.0f, sizeof(float) * (SIZE + 2) * (SIZE + 2));
	float *velocityY = (float *) calloc(1.0f, sizeof(float) * (SIZE + 2) * (SIZE + 2));

	float *velocityPrevX = (float *) calloc(1.0f, sizeof(float) * (SIZE + 2) * (SIZE + 2));
	float *velocityPrevY = (float *) calloc(1.0f, sizeof(float) * (SIZE + 2) * (SIZE + 2));

	//init(densityPrev, velocityPrevX, velocityPrevY);
	//printVel(velocityPrevX, velocityPrevY, "inital Velocity");

	float visc = 100.0f;
	float diff = 1.0f;
	float dt = 0.1f;

	addDensitySource(densityPrev);
	printArr(densityPrev, "densitySOURCE");
	drawDensity(densityPrev);


	addForce(velocityPrevX, velocityPrevY, 10.0f, 10.0f);
	for (int i = 1; i < 10; i++) {
		std::cout << "i: " << i << std::endl;
		vel_step(velocityX, velocityY, velocityPrevX, velocityPrevY, visc, dt);
		dens_step(density, densityPrev, velocityX, velocityY, diff, dt);

		// printArr(density, "density");
		drawDensity(density);
		// addDensitySource(densityPrev);
		SWAP(densityPrev, density);
		// SWAP(velocityPrevX, velocityX);
		// SWAP(velocityPrevY, velocityY);
		dt++;
	}

	
	std::cout << std::endl << std::endl << std::endl;
}

void addDensitySource(float *density) {
	for (int i = 100; i < 250; i++) {
		for (int j = 100; j < 250; j++) {
			density[at(i, j)] = 1.0f;
		}	
	}
	for (int i = 50; i < 60; i++) {
		for (int j = 0; j < 30; j++) {
			density[at(2+i, 70 + j)] = 1.0f;
			density[at(3+i, 70 + j)] = 1.0f;
			density[at(4+i, 70 + j)] = 1.0f;
			density[at(3+i, 70 + j)] = 1.0f;
			density[at(3+i, 70 + j)] = 1.0f;		
		}
	}
}
void addForce(float *velocityX, float *velocityY, float forceX, float forceY) {
	for (int i = START; i < END; i++) {
		for (int j = START; j < END; j++) {
			velocityX[at(i,j)] = forceX;
			velocityY[at(i,j)] = forceY;
		}
	}
}

void init(float *density, float *velocityX, float *velocityY) {
	for (int i = START; i < END; i++) {
		for (int j = START; j < END; j++) {
			velocityX[at(i,j)] = 0.0f;
			velocityY[at(i,j)] = -2.0f;
		}
	}
}

void vel_step(float *velX, float *velY, float *velPrevX, float *velPrevY, float visc, float dt) {
	// printVel(velX, velY, "initial");
	add_source(velX, velPrevX, dt);
	add_source(velY, velPrevY, dt);
	// printVel(velX, velY, "after addSources");
	
	SWAP(velPrevX, velX);
	diffuse(1, velX, velPrevX, visc, dt);

	SWAP(velPrevY, velY);
	diffuse(2, velY, velPrevY, visc, dt);
	
	project(velX, velY, velPrevX, velPrevY);
	// printVel(velX, velY, "after diffuse");

	SWAP(velPrevX, velX);
	SWAP(velPrevY, velY);
	
	advect(1, velX, velPrevX, velPrevX, velPrevY, dt);
	advect(2, velY, velPrevY, velPrevX, velPrevY, dt);
	project(velX, velY, velPrevX, velPrevY);
	// printVel(velX, velY, "after advect");
}

void dens_step(float *density, float *densityPrev, float *velX, float *velY, float diff, float dt) {
	// printArr(densityPrev, "initial Density");
	
	add_source(density, densityPrev, dt);
	SWAP(densityPrev, density);
	// printArr(densityPrev, "density after add sources");
	
	diffuse(0, density, densityPrev, diff, dt);
	SWAP(densityPrev, density);
	// printArr(densityPrev, "density after diffuse");

	advect(0, density, densityPrev, velX, velY, dt);
	// printArr(density, "density final");
}

void add_source(float *result, float *source, float dt) {
	int size = (SIZE + 2) * (SIZE + 2);
	for (int i = 0; i < size; i++) {
		result[i] += dt * source[i];
	}
}

void diffuse(int b, float *x, float *x0, float diff, float dt) {
	float a = dt * diff * SIZE * SIZE;

	for (int k = 0; k < 20; k++) {
		for (int i = START; i < END; i++) {
			for (int j = START; j < END; j++) {
				float surrDensity = x[at(i-1,j)] + x[at(i+1, j)] + x[at(i,j-1)] + x[at(i,j+1)];
				x[at(i,j)] = (x0[at(i,j)] + a * surrDensity) / (1 + 4 * a);
			}
		}
		set_bnd(b, x);
	}
}

void advect(int b, float *d, float *d0, float *u, float *v, float dt) {
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt * SIZE;

	for (i = START; i < END; i++) {
		for (j = START; j < END; j++) {
			x = i - dt0 * u[at(i,j)];
			y = j - dt0 * v[at(i,j)];

			//set x vals
			if (x < 0.5f) {
				x = 0.5f;
			}
			if (x > SIZE + 0.5f) {
				x = SIZE + 0.5f;
			}
			i0 = (int) x;
			i1 = i0 + 1.0f;

			s1 = x - i0;
			s0 = 1.0f - s1;

			//set y vals
			if (y < 0.5f) {
				y = 0.5f;
			}
			if (y > SIZE + 0.5f) {
				y = SIZE + 0.5f;
			}
			j0 = (int) y;
			j1 = j0 + 1.0f;

			t1 = y - j0;
			t0 = 1.0f - t1;

			d[at(i,j)] = s0 * (t0 * d0[at(i0, j0)] + t1 * d0[at(i0,j1)]) + 
						 s1 * (t0 * d0[at(i1, j0)] + t1 * d0[at(i1,j1)]);
		}
	}
	set_bnd(b, d);
}

void project(float *u, float *v, float *p, float *divD) {
	int i, j, k;
	float h;

	h = 1.0/SIZE;

	for (i = START; i < END; i++) {
		for (j = START; j < END; j++) {
			float term = u[at(i+1,j)] - u[at(i-1, j)] + v[at(i,j+1)] - v[at(i,j-1)];
			divD[at(i,j)] = -0.5f * h * term;
			p[at(i,j)] = 0.0f;
		}
	}
	set_bnd(0, divD);
	set_bnd(0, p);

	for (k = 0; k < 20; k++) {
		for (i = START; i < END; i++) {
			for (j = START; j < END; j++) {
				float term = p[at(i+1,j)] + p[at(i-1, j)] + p[at(i,j+1)] + p[at(i,j-1)];
				p[at(i,j)] = (divD[at(i,j)] + term) / 4;
			}
		}
		set_bnd(0, p);
	}

	for (i = START; i < END; i++) {
		for (j = START; j < END; j++) {
			u[at(i,j)] -= 0.5f * (p[at(i+1,j)] - p[at(i-1,j)])/h;
			v[at(i,j)] -= 0.5f * (p[at(i,j+1)] - p[at(i,j-1)])/h;
		}
	}	
	set_bnd(1, u);
	set_bnd(2, v);
}

void set_bnd (int b, float *x) {
	for (int i = 1; i <= SIZE; i++) {
		x[at(0, i)] 	 = b==1 ? -x[at(1, i)] 	  : x[at(1,i)];
		x[at(SIZE+1, i)] = b==1 ? -x[at(SIZE, i)] : x[at(SIZE,i)];
		x[at(i, 0)]      = b==2 ? -x[at(i, 1)]    : x[at(i, 1)];
		x[at(i, SIZE+1)] = b==2 ? -x[at(i, SIZE)] : x[at(i,SIZE)];
	}
	x[at(0,0)] 			 = 0.5f * (x[at(1,0)] + x[at(0,1)]);
	x[at(0,SIZE+1)]      = 0.5f * (x[at(1,SIZE+1)] + x[at(0,SIZE)]);
	x[at(SIZE+1,0)]      = 0.5f * (x[at(SIZE,0)] + x[at(SIZE+1,1)]);
	x[at(SIZE+1,SIZE+1)] = 0.5f * (x[at(SIZE,SIZE+1)] + x[at(SIZE+1,SIZE)]);
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
		std::cout << "j = "<< j << ": ";
		for (int i = START; i < END; i++) {
			std::cout << arr[at(i, j)] << " ";
		}
		std::cout << std::endl;
	}

	// for(int i = 0; i < (SIZE +2) * (SIZE +2); i++) {
	// 	std::cout << arr[i] << " ";
	// }
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
void printVel(float *velX, float *velY) {
	std::cout << "Velocity: " << std::endl;
	for (int j = START; j < END; j++) {
		std::cout << "j: ";
		for (int i = START; i < END; i++) {
			std::cout << "(" << velX[at(i, j)] << ", " << velY[at(i, j)] << ") ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void printVel(float *velX, float *velY, std::string name) {
	std::cout << "Velocity " << name << ": " << std::endl;
	for (int j = START; j < END; j++) {
		std::cout << "j: ";
		for (int i = START; i < END; i++) {
			std::cout << "(" << velX[at(i, j)] << ", " << velY[at(i, j)] << ") ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

Image img(SIZE + 2, SIZE + 2);

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
			// if (val >= .0001) {
				color.r = val;
				color.g = val;
				color.b = val;
				img.pixel(i, j, color);
			// }
			// else {
			// 	color.r = 1.0;
			// 	img.pixel(i, j, color);
			// }
		}
	}
	img.WriteTga((char *) "density.tga", true);
	system("open density.tga");
	usleep(TIME);
	system("rm density.tga");
}