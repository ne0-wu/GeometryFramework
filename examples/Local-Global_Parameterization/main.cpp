#include <iostream>
#include <filesystem>
#include <vector>

#include "Scene.h"
#include "Mesh.h"
#include "Window.h"

#include "GeometryProcessing.h"

// int main2()
// {
// 	Scene scene;

// 	// --------------------------------------------------

// 	// Load mesh and fit into unit ball
// 	Mesh mesh("meshes/camelhead.obj");
// 	mesh.fitIntoUnitBall();

// 	std::cout << mesh.numVertices() << " vertices"
// 			  << ", " << mesh.numEdges() << " edges"
// 			  << ", " << mesh.numFaces() << " faces"
// 			  << std::endl;

// 	Eigen::MatrixX2d result = localGlobalParameterization(mesh, 0, LocalGlobalTarget::ARAP);

// 	auto meshFrom2D = [](const Eigen::MatrixX2d &vertexList2D, const std::vector<unsigned int> &faceList) -> Mesh
// 	{
// 		std::vector<double> vertexList;
// 		for (int i = 0; i < vertexList2D.rows(); i++)
// 		{
// 			vertexList.push_back(vertexList2D(i, 0));
// 			vertexList.push_back(vertexList2D(i, 1));
// 			vertexList.push_back(0);
// 		}
// 		Mesh result(vertexList, faceList);
// 		result.fitIntoUnitBall();
// 		return result;
// 	};

// 	// --------------------------------------------------

// 	// Add the mesh to the scene
// 	scene.addMesh(mesh);

// 	// --------------------------------------------------

// 	auto initialStateKeyP = glfwGetKey(scene.window, GLFW_KEY_P);
// 	int numIter = -2;

// 	// render loop
// 	while (!scene.shouldClose())
// 	{
// 		// input
// 		// -----
// 		scene.processInput();

// 		if (glfwGetKey(scene.window, GLFW_KEY_P) == GLFW_RELEASE && initialStateKeyP == GLFW_PRESS)
// 		{
// 			scene.meshes.clear();
// 			numIter += 2;
// 			scene.addMesh(meshFrom2D(localGlobalParameterization(mesh, numIter, LocalGlobalTarget::ASAP), mesh.faceList()));
// 			scene.setCamera(defaultCamera);
// 		}
// 		initialStateKeyP = glfwGetKey(scene.window, GLFW_KEY_P);

// 		// render
// 		// ------
// 		scene.draw();

// 		scene.swapBuffers();
// 		scene.pollEvents();
// 	}

// 	return 0;
// }

int main()
{
	Window window(1280, 720, "Local-Global Parameterization");

	// Setup Dear ImGui style
	ImGui::StyleColorsDark();
	// ImGui::StyleColorsLight();

	//
	Mesh mesh("meshes/camelhead.obj");
	mesh.fitIntoUnitBall();

	std::cout << mesh.numVertices() << " vertices"
			  << ", " << mesh.numEdges() << " edges"
			  << ", " << mesh.numFaces() << " faces"
			  << std::endl;

	Scene scene;
	// scene.window = window.window;
	scene.addMesh(GLMesh(mesh));

	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glEnable(GL_DEPTH_TEST | GL_MULTISAMPLE);

	while (!window.shouldClose())
	{
		window.pollEvents();

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		{
			ImGui::Begin("Example");

			ImGui::BeginChild("Left Pane", ImVec2(200, 0), true);

			ImGui::Button("Button");
			static float value = 0.5f;
			ImGui::SliderFloat("Slider", &value, 0.0f, 1.0f);

			ImGui::EndChild();

			ImGui::SameLine();

			ImGui::BeginChild("Right Pane");

			ImGui::Text("OpenGL Scene Placeholder");
			ImGui::EndChild();

			ImGui::End();
		}

		ImGui::Render();

		glViewport(0, 0, 500, 500);
		scene.draw();

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		window.swapBuffers();
	}

	return 0;
}