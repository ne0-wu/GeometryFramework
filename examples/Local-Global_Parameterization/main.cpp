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
	Mesh mesh("meshes/camelhead.obj");
	mesh.fitIntoUnitBall();

	std::cout << mesh.numVertices() << " vertices"
			  << ", " << mesh.numEdges() << " edges"
			  << ", " << mesh.numFaces() << " faces"
			  << std::endl;

	Eigen::MatrixX2d result = localGlobalParameterization(mesh, 0, LocalGlobalTarget::ARAP);

	auto meshFrom2D = [](const Eigen::MatrixX2d &vertexList2D, const std::vector<unsigned int> &faceList) -> Mesh
	{
		std::vector<double> vertexList;
		for (int i = 0; i < vertexList2D.rows(); i++)
		{
			vertexList.push_back(vertexList2D(i, 0));
			vertexList.push_back(vertexList2D(i, 1));
			vertexList.push_back(0);
		}
		Mesh result(vertexList, faceList);
		result.fitIntoUnitBall();
		return result;
	};

	// std::vector<double> paraResult(mesh.numVertices() * 3);
	// for (int i = 0; i < mesh.numVertices(); i++)
	// {
	// 	paraResult[3 * i] = result(i, 0);
	// 	paraResult[3 * i + 1] = result(i, 1);
	// 	paraResult[3 * i + 2] = 0;
	// }

	// Mesh tutteTest(paraResult, mesh.faceList());

	// tutteTest.fitIntoUnitBall();

	// tutteTest.move(Eigen::Vector3d(0, 0, -1));

	// PointCloud pointCloud(mesh);
	// std::cout << "Point cloud size: " << pointCloud.data.size() << std::endl;
	// poissonSurfaceReconstruction(pointCloud);

	// --------------------------------------------------

	// Add the mesh to the scene
	// scene.addMesh(mesh);
	// scene.addMesh(tutteTest);
	// scene.addMesh(meshFrom2D(localGlobalParameterization(mesh, 0, LocalGlobalTarget::ARAP), mesh.faceList()));
	scene.addMesh(mesh);

	// --------------------------------------------------

	// auto initialStateKeyQ = glfwGetKey(scene.window, GLFW_KEY_Q);
	auto initialStateKeyP = glfwGetKey(scene.window, GLFW_KEY_P);
	int numIter = -2;

	// render loop
	while (!scene.shouldClose())
	{
		// input
		// -----
		scene.processInput();

		// if (glfwGetKey(scene.window, GLFW_KEY_Q) == GLFW_RELEASE && initialStateKeyQ == GLFW_PRESS)
		// {
		// 	scene.meshes[0].simplifyQEM((int)(scene.meshes[0].numVertices() * 0.9));
		// 	initialStateKeyQ = GLFW_RELEASE;
		// }
		// initialStateKeyQ = glfwGetKey(scene.window, GLFW_KEY_Q);

		if (glfwGetKey(scene.window, GLFW_KEY_P) == GLFW_RELEASE && initialStateKeyP == GLFW_PRESS)
		{
			scene.meshes.clear();
			numIter += 2;
			scene.addMesh(meshFrom2D(localGlobalParameterization(mesh, numIter, LocalGlobalTarget::ASAP), mesh.faceList()));
			scene.setCamera(defaultCamera);
		}
		initialStateKeyP = glfwGetKey(scene.window, GLFW_KEY_P);

		// render
		// ------
		scene.draw();

		scene.swapBuffers();
		scene.pollEvents();
	}

	return 0;
}