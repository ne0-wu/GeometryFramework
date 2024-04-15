#include <iostream>
#include <filesystem>
#include <vector>

#include "Scene.h"
#include "Mesh.h"
#include "Window.h"

#include "GeometryProcessing.h"

int main()
{
	Window window(1280, 720, "Local-Global Parameterization");

	//
	Mesh mesh("meshes/camelhead.obj");
	mesh.fitIntoUnitBall();

	std::cout << mesh.numVertices() << " vertices"
			  << ", " << mesh.numEdges() << " edges"
			  << ", " << mesh.numFaces() << " faces"
			  << std::endl;

	Scene scene;
	scene.window = window.window;
	scene.addMesh(GLMesh(mesh));

	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glfwWindowHint(GLFW_SAMPLES, 4);
	glEnable(GL_DEPTH_TEST);

	while (!window.shouldClose())
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		window.pollEvents();

		window.processInput();
		scene.processInput();

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		{
			ImGui::Begin("Local-Global Parameterization");

			int iterations = 0;
			ImGui::SliderInt("Iterations", &iterations, 1, 100);

			ImGui::End();
		}

		ImGui::Render();

		scene.draw();

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		window.swapBuffers();
	}

	return 0;
}