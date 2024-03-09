#include <iostream>
#include <filesystem>
#include <vector>

#include "Scene.h"
#include "Mesh.h"
#include "GLMesh.h"
#include "Shader.h"

int main()
{
	std::filesystem::path currentPath = std::filesystem::current_path();
	std::cout << "Current working directory: " << currentPath << std::endl;

	Scene scene;

	Mesh mesh("meshes/ball.obj");
	mesh.fitIntoUnitBall();
	mesh.scale(0.9);

	std::cout << "Number of vertices: " << mesh.numVertices() << std::endl;

	bool result;

	result = mesh.simplifyQEM((int)(mesh.numVertices() * 0.8));
	mesh.save("0.8.obj");

	if (result)
	{
		result = mesh.simplifyQEM((int)(mesh.numVertices() * 0.6));
		mesh.save("0.6.obj");
	}

	if (result)
	{
		result = mesh.simplifyQEM((int)(mesh.numVertices() * 0.4));
		mesh.save("0.4.obj");
	}

	std::cout << "Number of vertices: " << mesh.numVertices() << std::endl;

	GLMesh glMesh(mesh);

	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glEnable(GL_DEPTH_TEST);

	int vertexColorLocation = glGetUniformLocation(scene.shaderProgram.getID(), "color");

	// render loop
	while (!scene.shouldClose())
	{
		// input
		// -----
		scene.processInput();

		// render
		// ------
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glUniform4f(vertexColorLocation, 1.0f, 1.0f, 1.0f, 1.0f);
		glMesh.draw();

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(2.0f);
		glUniform4f(vertexColorLocation, 0.0f, 0.0f, 0.0f, 1.0f);
		glMesh.draw();

		scene.swapBuffers();
		scene.pollEvents();
	}

	return 0;
}