#pragma once

#include <vector>
#include <string>

#include <glad/glad.h>

class Mesh;

class GLMesh
{
private:
	GLuint VAO, VBO, EBO;

	void initializeMeshBuffers(const std::vector<GLfloat> &vertexList, const std::vector<GLuint> &faceList);

public:
	unsigned int numVertices, numFaces;

	GLMesh();
	GLMesh(const std::vector<GLfloat> &vertexList, const std::vector<GLuint> &faceList);
	GLMesh(const std::string &filename);
	GLMesh(Mesh &mesh);

	void draw();
	void update();
	void setMesh(const std::vector<GLfloat> &vertexList, const std::vector<GLuint> &faceList);
};