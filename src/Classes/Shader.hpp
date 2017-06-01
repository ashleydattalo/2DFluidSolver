#ifndef SHADER_HPP
#define SHADER_HPP

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
  
#include <GL/glew.h>

class Shader
{
public:
	// Constructor reads and builds the shader
	Shader(const GLchar* vertexPath, const GLchar* fragmentPath);
  	// Use the program
  	void use();
  	GLuint getProg();

private:
	GLuint program;
};
  
#endif