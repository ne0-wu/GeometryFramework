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

	Mesh mesh("meshes/spot_quadrangulated.obj");
	mesh.fitIntoUnitBall();
	mesh.scale(0.9);

	int numInitialVertices = mesh.numVertices();
	std::cout << "Number of vertices (initial): " << numInitialVertices << std::endl;

	bool result;

	for (double remainingRatio = 0.9; remainingRatio > 0.8; remainingRatio -= 0.1)
	{
		int targetNumVertices = (int)(numInitialVertices * remainingRatio);
		std::cout << "Number of vertices (" << remainingRatio << "): " << targetNumVertices << std::endl;
		mesh.simplifyQEM(targetNumVertices);
		mesh.save("output_initial=" + std::to_string(numInitialVertices) +
				  "_target=" + std::to_string(targetNumVertices) +
				  "_result=" + std::to_string(mesh.numVertices()) + ".obj");
	}

	// GLMesh glMesh(mesh);

	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glEnable(GL_DEPTH_TEST);

	int vertexColorLocation = glGetUniformLocation(scene.shaderProgram.getID(), "color");

	auto initialState = glfwGetKey(scene.window, GLFW_KEY_Q);

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
		mesh.render();
		// glMesh.draw();

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(2.0f);
		glUniform4f(vertexColorLocation, 0.0f, 0.0f, 0.0f, 1.0f);
		// glMesh.draw();
		mesh.render();

		if (glfwGetKey(scene.window, GLFW_KEY_Q) == GLFW_RELEASE && initialState == GLFW_PRESS)
		{
			mesh.simplifyQEM((int)(mesh.numVertices() * 0.9));
			initialState = GLFW_RELEASE;
		}
		else if (glfwGetKey(scene.window, GLFW_KEY_Q) == GLFW_PRESS)
		{
			initialState = GLFW_PRESS;
		}

		scene.swapBuffers();
		scene.pollEvents();
	}

	return 0;
}