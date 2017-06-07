#include <unistd.h>
#include <iostream>
#include <stdlib.h>
#include <glm/glm.hpp>

#include "Image.h"
#include "types.h"

#define TIME 2000000
#define SIZE 200
#define START 1
#define END SIZE + 1
#define at(i,j) ((i) + (SIZE + 2)*j)
// #define SWAP(x0,x) {float *tmp=x0;x0=x;x=tmp;}

#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define FOR_EACH_CELL for ( i=1 ; i<=N ; i++ ) { for ( j=1 ; j<=N ; j++ ) {
#define END_FOR }}


// void vel_step(float *velX, float *velY, float *velPrevX, float *velPrevY, float visc, float dt);
// void dens_step(float *density, float *densityPrev, float *velX, float *velY, float diff, float dt);

// void add_source(float *result, float *source, float dt);
// void diffuse(int b, float *x, float *x0, float diff, float dt);
// void advect(int b, float *d, float *d0, float *u, float *v, float dt);
// void project(float *u, float *v, float *p, float *divD);

// void set_bnd (int b, float *x);

void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt );
void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt );


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
void addDensity(float *density, float dt);

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
// int main() {

// 	float *density = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
// 	float *densityPrev = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

// 	float *velocityX = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
// 	float *velocityY = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

// 	float *velocityPrevX = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
// 	float *velocityPrevY = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

// 	//init(densityPrev, velocityPrevX, velocityPrevY);
// 	//printVel(velocityPrevX, velocityPrevY, "inital Velocity");

// 	float visc = 100.0f;
// 	float diff = 1.0f;
// 	float dt = 0.1f;

// 	addDensitySource(densityPrev);
// 	printArr(densityPrev, "densitySOURCE");
// 	drawDensity(densityPrev);

// 	addForce(velocityPrevX, velocityPrevY, 10.0f, 0.0f);
// 	for (int i = 1; i < 200; i++) {
// 		std::cout << "i: " << i << std::endl;
// 		addDensity(densityPrev, dt);
// 		vel_step(SIZE, velocityX, velocityY, velocityPrevX, velocityPrevY, visc, dt);
// 		dens_step(SIZE, density, densityPrev, velocityX, velocityY, diff, dt);

// 		// printArr(density, "density");
// 		// addDensitySource(densityPrev);
// 		SWAP(densityPrev, density);
// 		// SWAP(velocityPrevX, velocityX);
// 		// SWAP(velocityPrevY, velocityY);
// 		dt++;
// 		if (i % 10 == 0) {
// 			// drawDensity(density);
// 		}
// 	}
// 	drawDensity(density);

	
// 	std::cout << std::endl << std::endl << std::endl;
// }

int main() {
	//allocated data
	float *density = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
	float *velocityX = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
	float *velocityY = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

	float *densityPrev = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
	float *velocityPrevX = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
	float *velocityPrevY = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));


	float visc = 100.0f;
	float diff = 1.0f;
	float dt = 0.1f;

	addDensitySource(densityPrev);
	drawDensity(densityPrev);

	for (int i = 1; i < 100; i++) {
		// std::cout << "i: " << i << std::endl;

		// addForce(velocityPrevX, velocityPrevY, 10.0f, 0.0f);
		addForceAwayFromCenter(velocityPrevX, velocityPrevY);
		// addForce(velocityPrevX, velocityPrevY, 0, 2);
		vel_step(SIZE, velocityX, velocityY, velocityPrevX, velocityPrevY, visc, dt);
		dens_step(SIZE, density, densityPrev, velocityX, velocityY, diff, dt);

		if (i %10 == 0) {
			drawDensity(density);

		}
		// printArr(density, "density");

		
		dt++;
		if (i % 10 == 0) {
			drawDensity(density);
		}
	}
	drawDensity(density);
	printArr(density, "final density");
}

void addForce(float *velocityX, float *velocityY, float forceX, float forceY) {
    for (int i = 50; i < 150; i++) {
        for (int j = 50; j < 150; j++) {
            velocityX[at(i,j)] = forceX;
            velocityY[at(i,j)] = forceY;
        }
    }
    // int i, j;
    // int N = SIZE;
    // FOR_EACH_CELL
    //     velocityX[at(i,j)] = forceX;
    //     velocityY[at(i,j)] = forceY;
    // END_FOR
}

void addDensity(float *density, float dt) {
	float start = 0.0f;
	if (dt-20 > 0.0f) {
		start = dt - 20;
	}
	for (int i = start; i < dt + 20; i++) {
		density[at(i, i)] = 1.0f;
	}
}

void add_source ( int N, float * x, float * s, float dt )
{
	int i, size=(N+2)*(N+2);
	for ( i=0 ; i<size ; i++ ) x[i] += dt*s[i];
}

void set_bnd ( int N, int b, float * x )
{
	int i;

	for ( i=1 ; i<=N ; i++ ) {
		x[IX(0  ,i)] = b==1 ? -x[IX(1,i)] : x[IX(1,i)];
		x[IX(N+1,i)] = b==1 ? -x[IX(N,i)] : x[IX(N,i)];
		x[IX(i,0  )] = b==2 ? -x[IX(i,1)] : x[IX(i,1)];
		x[IX(i,N+1)] = b==2 ? -x[IX(i,N)] : x[IX(i,N)];
	}
	x[IX(0  ,0  )] = 0.5f*(x[IX(1,0  )]+x[IX(0  ,1)]);
	x[IX(0  ,N+1)] = 0.5f*(x[IX(1,N+1)]+x[IX(0  ,N)]);
	x[IX(N+1,0  )] = 0.5f*(x[IX(N,0  )]+x[IX(N+1,1)]);
	x[IX(N+1,N+1)] = 0.5f*(x[IX(N,N+1)]+x[IX(N+1,N)]);
}

void lin_solve ( int N, int b, float * x, float * x0, float a, float c )
{
	int i, j, k;

	for ( k=0 ; k<20 ; k++ ) {
		FOR_EACH_CELL
			x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;
		END_FOR
		set_bnd ( N, b, x );
	}
}

void diffuse ( int N, int b, float * x, float * x0, float diff, float dt )
{
	float a=dt*diff*N*N;
	lin_solve ( N, b, x, x0, a, 1+4*a );
}

void advect ( int N, int b, float * d, float * d0, float * u, float * v, float dt )
{
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt*N;
	FOR_EACH_CELL
		x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];
		if (x<0.5f) x=0.5f; if (x>N+0.5f) x=N+0.5f; i0=(int)x; i1=i0+1;
		if (y<0.5f) y=0.5f; if (y>N+0.5f) y=N+0.5f; j0=(int)y; j1=j0+1;
		s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
		d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+
					 s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)]);
	END_FOR
	set_bnd ( N, b, d );
}

void project ( int N, float * u, float * v, float * p, float * div )
{
	int i, j;

	FOR_EACH_CELL
		div[IX(i,j)] = -0.5f*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
		p[IX(i,j)] = 0;
	END_FOR	
	set_bnd ( N, 0, div ); set_bnd ( N, 0, p );

	lin_solve ( N, 0, p, div, 1, 4 );

	FOR_EACH_CELL
		u[IX(i,j)] -= 0.5f*N*(p[IX(i+1,j)]-p[IX(i-1,j)]);
		v[IX(i,j)] -= 0.5f*N*(p[IX(i,j+1)]-p[IX(i,j-1)]);
	END_FOR
	set_bnd ( N, 1, u ); set_bnd ( N, 2, v );
}

void dens_step ( int N, float * x, float * x0, float * u, float * v, float diff, float dt )
{
	add_source ( N, x, x0, dt );
	SWAP ( x0, x ); diffuse ( N, 0, x, x0, diff, dt );
	SWAP ( x0, x ); advect ( N, 0, x, x0, u, v, dt );
}

void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt )
{
	add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
	SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt );
	SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt );
	project ( N, u, v, u0, v0 );
	SWAP ( u0, u ); SWAP ( v0, v );
	advect ( N, 1, u, u0, u0, v0, dt ); advect ( N, 2, v, v0, u0, v0, dt );
	project ( N, u, v, u0, v0 );
}









// MY FUNCTIONS

void addDensitySource(float *density) {
	// for (int i = 100; i < 250; i++) {
	// 	for (int j = 100; j < 250; j++) {
	// 		density[at(i, j)] = 1.0f;
	// 	}	
	// }
	// for (int i = 50; i < 60; i++) {
	// 	for (int j = 0; j < 30; j++) {
	// 		density[at(2+i, 70 + j)] = 1.0f;
	// 		density[at(3+i, 70 + j)] = 1.0f;
	// 		density[at(4+i, 70 + j)] = 1.0f;
	// 		density[at(3+i, 70 + j)] = 1.0f;
	// 		density[at(3+i, 70 + j)] = 1.0f;		
	// 	}
	// }
	for (int i = 0; i < 50; i++) {
		for (int j = 0; j < 50; j++) {
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
// void addForce(float *velocityX, float *velocityY, float forceX, float forceY) {
// 	for (int i = START; i < END; i++) {
// 		for (int j = START; j < END; j++) {
// 			velocityX[at(i,j)] = forceX;
// 			velocityY[at(i,j)] = forceY;
// 		}
// 	}
// }

void init(float *density, float *velocityX, float *velocityY) {
	for (int i = START; i < END; i++) {
		for (int j = START; j < END; j++) {
			velocityX[at(i,j)] = 0.0f;
			velocityY[at(i,j)] = -2.0f;
		}
	}
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