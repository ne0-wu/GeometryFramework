#include "GLMesh.h"

#include "Mesh.h"

GLMesh::GLMesh()
{
	std::vector<GLfloat> vertexList = {
		0.5f, 0.5f, 0.0f,	// Top right
		0.5f, -0.5f, 0.0f,	// Bottom right
		-0.5f, -0.5f, 0.0f, // Bottom left
		-0.5f, 0.5f, 0.0f	// Top left
	};
	std::vector<GLuint> faceList = {
		0, 1, 3, // First triangle
		1, 2, 3	 // Second triangle
	};
	numVertices = vertexList.size() / 3;
	numFaces = faceList.size() / 3;
	initializeMeshBuffers(vertexList, faceList);
}

GLMesh::GLMesh(const std::vector<GLfloat> &vertexList, const std::vector<GLuint> &faceList)
{
	numVertices = vertexList.size() / 3;
	numFaces = faceList.size() / 3;
	initializeMeshBuffers(vertexList, faceList);
}

GLMesh::GLMesh(const std::string &filename)
{
	Mesh mesh(filename);

	std::vector<GLfloat> vertexList = mesh.vertexListFloat();
	std::vector<GLuint> faceList = mesh.faceList();

	numVertices = vertexList.size() / 3;
	numFaces = faceList.size() / 3;
	initializeMeshBuffers(vertexList, faceList);
}

GLMesh::GLMesh(Mesh &mesh)
{

	std::vector<GLfloat> vertexList = mesh.vertexListFloat();
	std::vector<GLuint> faceList = mesh.faceList();

	numVertices = vertexList.size() / 3;
	numFaces = faceList.size() / 3;
	initializeMeshBuffers(vertexList, faceList);
}

void GLMesh::initializeMeshBuffers(const std::vector<GLfloat> &vertexList, const std::vector<GLuint> &faceList)
{
	// 1. Generate vertex array object
	glGenVertexArrays(1, &VAO);
	// 2. Bind vertex array object
	glBindVertexArray(VAO);

	// 3. Generate vertex buffer object
	glGenBuffers(1, &VBO);
	// 4. Bind vertex buffer object
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	// 5. Copy vertex data to vertex buffer object
	glBufferData(GL_ARRAY_BUFFER, vertexList.size() * sizeof(GLfloat), vertexList.data(), GL_STATIC_DRAW);

	// 6. Generate index buffer object
	glGenBuffers(1, &EBO);
	// 7. Bind index buffer object
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	// 8. Copy index data to index buffer object
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, faceList.size() * sizeof(GLuint), faceList.data(), GL_STATIC_DRAW);

	// 9. Set vertex attribute pointers
	// Vertex position
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void *)0);
	glEnableVertexAttribArray(0);

	// 10. Unbind vertex array object
	glBindVertexArray(0);
}

void GLMesh::draw()
{
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, numFaces * 3, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}