#pragma once

#include <string>

#include <glad/gl.h>

class Shader
{
public:
	Shader() {}
	Shader(const std::string &vertexShaderSource, const std::string &fragmentShaderSource, bool fromFile = true);
	void use() { glUseProgram(ID); }

	static void checkCompileErrors(GLuint shader, const std::string &type);
	GLuint getID() const { return ID; }
	~Shader() {}

private:
	GLuint ID;
	static std::string readShaderCodeFromFile(const std::string &shaderPath);
	static GLuint compileShader(const GLchar *shaderCode, GLenum shaderType);
};