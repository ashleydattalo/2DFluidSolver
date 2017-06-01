#ifndef TEXTURE_HPP
#define TEXTURE_HPP

#include <GL/glew.h>
#include <string>

class Texture {
public:
    Texture(std::string texname);
    ~Texture();

    void loadTexture();
    void setTarget(GLenum target);
    void bind();
    void unbind();
    
private:
    GLuint handle, target;

    std::string texname;

    static GLuint createTexture2D(std::string texname);
};

#endif