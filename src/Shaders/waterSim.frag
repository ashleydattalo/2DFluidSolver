
#version 330 core

in float densityColor;
out vec4 color;

void main()
{
    vec3 density = vec3(densityColor);
    color = vec4(densityColor*.2, densityColor * .4, densityColor *.4, 1.0f);
}