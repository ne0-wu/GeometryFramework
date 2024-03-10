#include <iostream>
#include <filesystem>
#include <vector>

#include "Scene.h"
#include "Mesh.h"

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

	// for (double remainingRatio = 0.9; remainingRatio > 0.8; remainingRatio -= 0.1)
	// {
	// 	int targetNumVertices = (int)(numInitialVertices * remainingRatio);
	// 	std::cout << "Number of vertices (" << remainingRatio << "): " << targetNumVertices << std::endl;
	// 	mesh.simplifyQEM(targetNumVertices);
	// 	mesh.save("output_initial=" + std::to_string(numInitialVertices) +
	// 			  "_target=" + std::to_string(targetNumVertices) +
	// 			  "_result=" + std::to_string(mesh.numVertices()) + ".obj");
	// }

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glEnable(GL_DEPTH_TEST);

	auto initialState = glfwGetKey(scene.window, GLFW_KEY_Q);

	scene.addMesh(mesh);

	// render loop
	while (!scene.shouldClose())
	{
		// input
		// -----
		scene.processInput();

		if (glfwGetKey(scene.window, GLFW_KEY_Q) == GLFW_RELEASE && initialState == GLFW_PRESS)
		{
			scene.meshes[0].simplifyQEM((int)(scene.meshes[0].numVertices() * 0.9));
			initialState = GLFW_RELEASE;
		}
		initialState = glfwGetKey(scene.window, GLFW_KEY_Q);

		// render
		// ------
		scene.draw();

		scene.swapBuffers();
		scene.pollEvents();
	}

	return 0;
}