#include <iostream>
#include <fstream>
#include <sstream>

#include "Shader.h"

Shader::Shader(const std::string &vertexShaderCode, const std::string &fragmentShaderCode, bool fromFile)
{
	GLuint vertexShader, fragmentShader;
	if (fromFile)
	{
		vertexShader = compileShader(readShaderCodeFromFile(vertexShaderCode).c_str(), GL_VERTEX_SHADER);
		checkCompileErrors(vertexShader, "VERTEX");
		fragmentShader = compileShader(readShaderCodeFromFile(fragmentShaderCode).c_str(), GL_FRAGMENT_SHADER);
		checkCompileErrors(fragmentShader, "FRAGMENT");
	}
	else
	{
		vertexShader = compileShader(vertexShaderCode.c_str(), GL_VERTEX_SHADER);
		fragmentShader = compileShader(fragmentShaderCode.c_str(), GL_FRAGMENT_SHADER);
	}

	ID = glCreateProgram();
	glAttachShader(ID, vertexShader);
	glAttachShader(ID, fragmentShader);
	glLinkProgram(ID);

	checkCompileErrors(ID, "PROGRAM");

	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);
}

void Shader::use()
{
	glUseProgram(ID);
}

GLuint Shader::getID() const
{
	return ID;
}

void Shader::checkCompileErrors(GLuint shader, const std::string &type)
{
	GLint success;
	GLchar infoLog[1024];
	if (type != "PROGRAM")
	{
		glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
		if (!success)
		{
			glGetShaderInfoLog(shader, 1024, NULL, infoLog);
			std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << "\n"
					  << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
		}
	}
	else
	{
		glGetProgramiv(shader, GL_LINK_STATUS, &success);
		if (!success)
		{
			glGetProgramInfoLog(shader, 1024, NULL, infoLog);
			std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << "\n"
					  << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
		}
	}
}

std::string Shader::readShaderCodeFromFile(const std::string &shaderPath)
{
	std::ifstream shaderFile;
	shaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	std::stringstream shaderStream;
	try
	{
		shaderFile.open(shaderPath);
		shaderStream << shaderFile.rdbuf();
		shaderFile.close();
	}
	catch (std::ifstream::failure &e)
	{
		std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
	}
	return shaderStream.str();
}

GLuint Shader::compileShader(const GLchar *shaderCode, GLenum shaderType)
{
	GLuint shader;
	shader = glCreateShader(shaderType);
	glShaderSource(shader, 1, &shaderCode, NULL);
	glCompileShader(shader);
	return shader;
}