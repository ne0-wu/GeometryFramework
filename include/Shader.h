#pragma once

#include <string>
#include <glad/glad.h>

class Shader
{
public:
	Shader(const std::string &vertexShaderSource, const std::string &fragmentShaderSource, bool fromFile = true);
	void use();
	static void checkCompileErrors(GLuint shader, const std::string &type);
	GLuint getID() const;
	~Shader();

private:
	GLuint ID;
	static std::string readShaderCodeFromFile(const std::string &shaderPath);
	static GLuint compileShader(const GLchar *shaderCode, GLenum shaderType);
};