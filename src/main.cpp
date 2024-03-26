#include <iostream>
#include <filesystem>
#include <vector>

#include "Scene.h"
#include "Mesh.h"

#include "GeometryProcessing.h"

int main()
{
	std::filesystem::path currentPath = std::filesystem::current_path();
	std::cout << "Current working directory: " << currentPath << std::endl;

	Scene scene;

	// --------------------------------------------------

	// Load mesh and fit into unit ball
	Mesh mesh("meshes/spot_simplified.obj");
	mesh.fitIntoUnitBall();

	std::cout << mesh.numVertices() << " vertices"
			  << ", " << mesh.numEdges() << " edges"
			  << ", " << mesh.numFaces() << " faces"
			  << std::endl;

	PointCloud pointCloud(mesh);

	std::cout << "Point cloud size: " << pointCloud.data.size() << std::endl;

	// poissonSurfaceReconstruction(pointCloud);

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