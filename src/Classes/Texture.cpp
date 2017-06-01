#include <GL/glew.h>
#include <string>
#include <iostream>
#include "Texture.hpp"

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

Texture::Texture(std::string texname) {
    this->target = GL_TEXTURE0;
    this->texname = texname;
}

Texture::~Texture() {
    glDeleteTextures(1, &handle);
}

void Texture::loadTexture() {
    this->handle = createTexture2D(texname);
}

void Texture::setTarget(GLenum target) {
    this->target =  target;
}

void Texture::bind() {
    glActiveTexture(target);
    glBindTexture(GL_TEXTURE_2D, handle);
}

void Texture::unbind() {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

GLuint Texture::createTexture2D(std::string texname) {
    int width, height, channels;
    unsigned char *image = stbi_load(texname.c_str(), &width, &height, &channels, 0);
    if (image == NULL) {
        std::cout
            << "ERROR::TEXTURE::LOAD_FAILED::"
            << texname
            << std::endl;
    }

    GLuint texture;
    glGenTextures(1, &texture);

    glBindTexture(GL_TEXTURE_2D, texture);
    if (channels == 4)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
    else
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);

    glGenerateMipmap(GL_TEXTURE_2D);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    stbi_image_free(image);
    glBindTexture(GL_TEXTURE_2D, 0);

    return texture;
}