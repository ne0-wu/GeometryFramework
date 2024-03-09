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

	Shader shaderProgram("shaders/basic.vert", "shaders/basic.frag");
	scene.shaderProgram = &shaderProgram;

	Mesh mesh("meshes/cathead.obj");
	mesh.fitIntoUnitBall();
	mesh.scale(0.9);

	std::cout << "Number of vertices: " << mesh.numVertices() << std::endl;

	// mesh.simplifyQEM(120);

	std::cout << "Number of vertices: " << mesh.numVertices() << std::endl;

	GLMesh glMesh(mesh);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	// glEnable(GL_DEPTH_TEST);

	// render loop
	while (!scene.shouldClose())
	{
		// scene.shaderProgram.use();

		// input
		// -----
		scene.processInput();

		// shaderProgram.use();
		scene.shaderProgram->use();

		// render
		// ------
		glClear(GL_COLOR_BUFFER_BIT);

		glMesh.draw();

		scene.swapBuffers();
		scene.pollEvents();
	}

	return 0;
}