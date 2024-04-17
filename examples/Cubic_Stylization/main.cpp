#include <iostream>
#include <filesystem>
#include <vector>

#include "Scene.h"
#include "Window.h"

#include "Mesh.h"
#include "GeometryProcessing.h"

int main()
{
	Window window(1280, 720, "Cubic Stylization");

	Scene scene;
	scene.window = window.window;

	Mesh original_mesh("meshes/spot_quadrangulated.obj");
	original_mesh.fitIntoUnitBall();

	Mesh mesh = original_mesh;
	scene.addMesh(GLMesh(mesh));

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glfwWindowHint(GLFW_SAMPLES, 4);
	glEnable(GL_DEPTH_TEST);

	// GUI states
	float lambda = 0.1, lambda_prev = 0.09;
	int num_iter = 10, num_iter_prev = 10;

	while (!window.shouldClose())
	{
		window.pollEvents();
		window.processInput();
		scene.processInput();

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		{
			ImGui::Begin("Cubic Stylization");

			ImGui::SliderFloat("Lambda", &lambda, 0.0, 1.0);
			ImGui::SliderInt("Iterations", &num_iter, 1, 100);

			ImGui::End();
		}

		ImGui::Render();

		scene.update();
		scene.draw();

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		window.swapBuffers();

		if (lambda != lambda_prev || num_iter != num_iter_prev)
		{
			lambda_prev = lambda;
			num_iter_prev = num_iter;

			CubicStylization cubic(original_mesh, lambda, num_iter);
			cubic.stylize();
			mesh = cubic.get_stylized_mesh();
			mesh.fitIntoUnitBall();

			scene.glMeshes[0].setMesh(std::make_shared<Mesh>(mesh));
		}
	}

	return 0;
}