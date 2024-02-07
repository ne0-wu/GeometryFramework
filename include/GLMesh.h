#pragma once

#include <glad/glad.h>
#include <Eigen/Core>

#include <vector>

#include "Mesh.h"

class GLMesh
{
private:
	GLuint VAO, VBO, EBO;

	void initializeMeshBuffers(const std::vector<GLfloat> &vertexList, const std::vector<GLuint> &faceList);

public:
	unsigned int numVertices, numFaces;

	GLMesh(const std::vector<GLfloat> &vertexList, const std::vector<GLuint> &faceList);
	GLMesh(Mesh &mesh);
	void draw();
};