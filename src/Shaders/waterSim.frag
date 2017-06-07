
#version 330 core

in float densityColor;
out vec4 color;

void main()
{
    vec3 density = vec3(densityColor);
    color = vec4(densityColor*.7, densityColor * .3, densityColor *.2, 1.0f);
}