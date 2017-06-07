
#version 330 core

in float densityColor;
out vec4 color;

void main()
{
    vec3 density = vec3(densityColor);
    color = vec4(densityColor*.2, densityColor * .9, densityColor *.9, 1.0f);
}