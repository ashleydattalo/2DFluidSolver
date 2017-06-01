#version 330 core
layout (location = 0) in vec4 aPos;
layout (location = 1) in vec3 aCol;
layout (location = 2) in float aAlp;
layout (location = 3) in float aSca;
layout (location = 4) in float radius;
layout (location = 5) in float initialTime;
uniform mat4 P;
uniform mat4 V;
uniform float time;

out vec4 vecCol;
out vec4 vecPos;

vec4 getPos();

void main()
{
	//aPos.x + abs(aPos.x);
	//vec4 tempPos = vec4(aPos.x, aPos.y, aPos.z, aPos.w);
    //gl_Position = P * V * tempPos;
    //vecPos = P * V * tempPos;


    vec4 newPos;
    float newTime = initialTime + time;
	newPos.x = radius * sin(time);
	newPos.y = radius * (time);
    newPos.z = radius * cos(time);

    newPos += aPos;
    
    gl_Position = P * V * newPos;
    vecPos = P * V * newPos;
    

    vecCol = vec4(aCol, aAlp);
    gl_PointSize = aSca;
}

vec4 getPos() {
	vec4 position = aPos;
	position.x += radius * cos(initialTime + time);
	position.y += radius * (initialTime + time);
    position.z += radius * sin(initialTime + time);
    return position;
}