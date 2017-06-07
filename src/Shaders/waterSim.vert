#version 330 core
layout (location = 0) in vec3 position;
//layout (location = 1) in float density;

out vec3 densityColor;

void main()
{
    gl_Position = vec4(position.x, position.y, 0.0, 1.0);
    densityColor = vec3(position.x, position.y, position.z);
}