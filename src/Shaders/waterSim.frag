
#version 330 core

in vec3 densityColor;
out vec4 color;

void main()
{
    
    //rainbow
    //color = vec4(densityColor, 1.0f);
    

    //top stays colored
    //color = vec4( (densityColor.y+1)/2 * densityColor.z, (densityColor.y+1)/2 * densityColor.z * .2, (densityColor.y+1)/2 * densityColor.z* .2, 1.0f);

    //sin ~ looks really cool!
    //color = vec4(sin(densityColor.z), cos(densityColor.z), densityColor.z * .3, 1.0f);


    // turquiose blue
    //color = vec4(densityColor.x*.6, densityColor.y * .9, densityColor.x *.9, 1.0f);


    //on fire!
    color = vec4(densityColor.z*.7, densityColor.z * .3, densityColor.z *.2, 1.0f);

    //normal
    //color = vec4(vec3(densityColor.z), 1.0f);
}

