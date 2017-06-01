#version 330 core
out vec4 color;
uniform vec4 ourColor;
uniform float time;

in vec4 vecCol;
in vec4 vecPos;

float getCol(float x, float y);

void main()
{
    vec2 cord = 2.0 * gl_PointCoord - 1.0;
    if (dot(cord, cord) <= 1.0f) {
    	float x = cord.x;
    	float y = cord.y;

    	color = vec4(vecPos.xyz, getCol(x,y));
    	color = vecCol;
    }
    else {
        color = vec4(1.0f, 1.0f, 1.0f, 0.0f);
    }
}

float getCol(float x, float y) {
	return x*x + y*y;
}