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

	// --------------------------------------------------

	// Load mesh and fit into unit ball
	Mesh mesh("meshes/cathead.obj");
	mesh.fitIntoUnitBall();
	mesh.scale(0.9);

	std::cout << mesh.numVertices() << " vertices" << std::endl;

	// Point cloud test

	PointCloud pointCloud = mesh.pointCloud(200, true);

	std::cout << "Point cloud size: " << pointCloud.points.size() << std::endl;
	// for (auto p : pointCloud.points)
	// {
	// 	std::cout << p << std::endl;
	// }
	// for (auto n : pointCloud.normals)
	// {
	// 	std::cout << n << std::endl;
	// }

	testOctree(pointCloud);

	// --------------------------------------------------

	// Add the mesh to the scene
	scene.addMesh(mesh);

	// --------------------------------------------------

	auto initialState = glfwGetKey(scene.window, GLFW_KEY_Q);

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