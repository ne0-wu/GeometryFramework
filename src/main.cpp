#include <iostream>
#include <filesystem>
#include <vector>

#include <glad/glad.h>

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

	Mesh mesh("meshes/cathead.obj");
	mesh.fitIntoUnitBall();
	mesh.scale(0.9);
	GLMesh glMesh(mesh);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

	// render loop
	while (!scene.shouldClose())
	{
		// input
		// -----
		scene.processInput();

		// render
		// ------
		glClear(GL_COLOR_BUFFER_BIT);

		shaderProgram.use();

		glMesh.draw();

		scene.swapBuffers();
	}

	return 0;
}