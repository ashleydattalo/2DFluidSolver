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
#define at(i,j) ((i) + (SIZE + 2)*(j))
#define SWAP(x0,x) {float *tmp=x0;x0=x;x=tmp;}

void dens_step(int N, float * x, float * x0, float * u, float * v, float diff, float dt);
void vel_step(int N, float * u, float * v, float * u0, float * v0, float visc, float dt);

void drawDensity(float *density);
void addDensity(float *density);
void addForce(float *velocityX, float *velocityY, float forceX, float forceY);
void addForceAwayFromCenter(float *velocityX, float *velocityY);

void printArr(float *arr, std::string arrName);

Image img(SIZE + 2, SIZE + 2);

int main() {
	//allocated data
	float *density = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
	float *densityPrev = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

	float *velocityX = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
	float *velocityY = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));

	float *velocityPrevX = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));
	float *velocityPrevY = (float *) malloc(sizeof(float) * (SIZE + 2) * (SIZE + 2));


	float visc = 100.0f;
	float diff = 1.0f;
	float dt = 0.1f;

	addDensity(densityPrev);
	// drawDensity(densityPrev);

	for (int i = 1; i < 30; i++) {
		std::cout << "i: " << i << std::endl;

		// addForce(velocityPrevX, velocityPrevY, 10.0f, 0.0f);
		// addForceAwayFromCenter(velocityPrevX, velocityPrevY);

		vel_step(SIZE, velocityX, velocityY, velocityPrevX, velocityPrevY, visc, dt);
		dens_step(SIZE, density, densityPrev, velocityX, velocityY, diff, dt);

		if (i % 10 == 0) {
			drawDensity(density);
		}
		dt++;
	}
	drawDensity(density);
	printArr(density, "final density");
}

void add_source(int N, float *x, float *s, float dt) {
	int size = (N + 2) * (N + 2);
	for (int i = 0; i < size; i++) {
		x[i] += dt * s[i];
	}
}

void set_bnd(int N, int b, float * x ) {
	int i;
	for ( i=1 ; i<=N ; i++ ) {
		x[at(0, i)] 	 = b==1 ? -x[at(1, i)] : x[at(1,i)];
		x[at(N+1, i)]    = b==1 ? -x[at(N, i)] : x[at(N,i)];
		x[at(i, 0)]      = b==2 ? -x[at(i, 1)] : x[at(i, 1)];
		x[at(i, N+1)]    = b==2 ? -x[at(i, N)] : x[at(i,N)];
	}
	x[at(0 ,0 )] = 0.5*(x[at(1,0 )]+x[at(0 ,1)]);
	x[at(0 ,N+1)] = 0.5*(x[at(1,N+1)]+x[at(0 ,N )]);
	x[at(N+1,0 )] = 0.5*(x[at(N,0 )]+x[at(N+1,1)]);
	x[at(N+1,N+1)] = 0.5*(x[at(N,N+1)]+x[at(N+1,N)]);
}

void diffuse(int N, int b, float *x, float *x0, float diff, float dt) {
	int i, j, k;
	float a = dt * diff * N * N;

	for (k = 0; k < 20; k++) {
		for(int i = 1; i <= N; i++) {
			for (int j = 0; j <= N; j++) {
				float surrDensity = x[at(i-1,j)] + x[at(i+1,j)] + x[at(i,j-1)] + x[at(i,j+1)];
				x[at(i,j)] = (x0[at(i,j)] + a * surrDensity) / (1 + 4*a);
			}
		}
		set_bnd(N,b,x);
	}
}

void advect(int N, int b, float *d, float *d0, float *u, float *v, float dt) {
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt * N;
	for(int i = 1; i <= N; i++) {
		for (int j = 0; j <= N; j++) {
			x = i-dt0*u[at(i,j)]; y = j-dt0*v[at(i,j)];
			if (x<0.5) x=0.5; if (x>N+0.5) x=N+0.5; i0=(int)x; i1=i0+1;
			if (y<0.5) y=0.5; if (y>N+0.5) y=N+0.5; j0=(int)y; j1=j0+1; 
			s1 = x-i0; s0 = 1-s1; 
			t1 = y-j0; t0 = 1-t1;
			d[at(i,j)] = s0*(t0*d0[at(i0,j0)]+t1*d0[at(i0,j1)])+s1*(t0*d0[at(i1,j0)]+t1*d0[at(i1,j1)]);
		}
	}
	set_bnd(N, b, d);
}

void project ( int N, float * u, float * v, float * p, float * div ) {
	int i, j, k;
	float h;

	h = 1.0 / N;
	for(int i = 1; i <= N; i++) {
		for (int j = 0; j <= N; j++) {
			div[at(i,j)] = -0.5*h*(u[at(i+1,j)]-u[at(i-1,j)]+v[at(i,j+1)]-v[at(i,j-1)]);
			p[at(i,j)] = 0;
		}
	}

	set_bnd(N, 0, div);
	set_bnd(N, 0, p);

	for (k = 0; k < 20; k++) {
		for(int i = 1; i <= N; i++) {
			for (int j = 0; j <= N; j++) {
				p[at(i,j)] = (div[at(i,j)]+p[at(i-1,j)]+p[at(i+1,j)]+p[at(i,j-1)]+p[at(i,j+1)])/4;
			}
		}
		set_bnd(N, 0, p);
	}

	for(int i = 1; i <= N; i++) {
		for (int j = 0; j <= N; j++) {
		u[at(i,j)] -= 0.5*(p[at(i+1,j)]-p[at(i-1,j)])/h;
	    v[at(i,j)] -= 0.5*(p[at(i,j+1)]-p[at(i,j-1)])/h;
	}
}
	set_bnd(N, 1, u); 
	set_bnd(N, 2, v);
}

void dens_step(int N, float * x, float * x0, float * u, float * v, float diff, float dt) {
	add_source(N, x, x0, dt);
	SWAP(x0, x); 
	diffuse(N, 0, x, x0, diff, dt);
	SWAP(x0, x);
	advect(N, 0, x, x0, u, v, dt);
}

void vel_step ( int N, float * u, float * v, float * u0, float * v0, float visc, float dt ) {
	add_source ( N, u, u0, dt ); add_source ( N, v, v0, dt );
	SWAP ( u0, u ); diffuse ( N, 1, u, u0, visc, dt );
	SWAP ( v0, v ); diffuse ( N, 2, v, v0, visc, dt );
	project ( N, u, v, u0, v0 );
	SWAP ( u0, u ); SWAP ( v0, v );
	advect ( N, 1, u, u0, u0, v0, dt ); advect ( N, 2, v, v0, u0, v0, dt ); project ( N, u, v, u0, v0 );
}


//adding sources to scene
void addDensity(float *density) {
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
	// density[at(5,5)] = 1.0f;
	// density[at(3,5)] = 1.0f;
	// density[at(4,5)] = 1.0f;
	// density[at(5,3)] = 1.0f;
	// density[at(5,4)] = 1.0f;
	// density[at(5,5)] = 1.0f;
}
void addForce(float *velocityX, float *velocityY, float forceX, float forceY) {
	for (int i = START; i < END; i++) {
		for (int j = START; j < END; j++) {
			velocityX[at(i,j)] = forceX;
			velocityY[at(i,j)] = forceY;
		}
	}
}

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

// Drawing the density
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
			color.r = val;
			color.g = val;
			color.b = val;
			img.pixel(i, j, color);
		}
	}
	img.WriteTga((char *) "density.tga", true);
	system("open density.tga");
	usleep(TIME);
	system("rm density.tga");
}

//debugging
void printArr(float *arr, std::string arrName) {
	std::cout << arrName << ": " << std::endl;
	for (int j = START; j < END; j++) {
		std::cout << "j = "<< j << ": ";
		for (int i = START; i < END; i++) {
			std::cout << arr[at(i, j)] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}