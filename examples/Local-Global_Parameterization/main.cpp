#include <iostream>
#include <filesystem>
#include <vector>

#include "Scene.h"
#include "Mesh.h"
#include "Window.h"

#include "GeometryProcessing.h"

void flattenMesh(Mesh &mesh, Eigen::MatrixX2d uv)
{
	for (auto v : mesh.vertices())
	{
		auto uv2d = uv.row(v.idx());
		mesh.set_point(v, Mesh::Point(uv2d(0), uv2d(1), 0));
	}
}

int main()
{
	Window window(1280, 720, "Local-Global Parameterization");

	Mesh mesh("meshes/camelhead.obj");
	mesh.fitIntoUnitBall();

	auto localGlobal = LocalGlobalParameterization(mesh);
	localGlobal.flatten();

	std::cout << mesh.numVertices() << " vertices, "
			  << mesh.numEdges() << " edges, "
			  << mesh.numFaces() << " faces"
			  << std::endl;

	Scene scene;
	scene.window = window.window;

	auto mesh2d = mesh;
	flattenMesh(mesh2d, localGlobal.getUV());
	mesh2d.moveCenterToOrigin();

	scene.addMesh(GLMesh(mesh2d));

	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glfwWindowHint(GLFW_SAMPLES, 4);
	glEnable(GL_DEPTH_TEST);

	// GUI states
	int numIter = 0;

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

			ImGui::SliderInt("Iterations", &numIter, 1, 100);

			ImGui::End();
		}

		ImGui::Render();

		scene.update();
		scene.draw();

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		window.swapBuffers();

		if (numIter != localGlobal.getNumIter())
		{
			localGlobal.setNumIter(numIter);
			localGlobal.flatten();
			flattenMesh(mesh2d, localGlobal.getUV());
			mesh2d.moveCenterToOrigin();
			scene.glMeshes[0].shouldUpdate = true;
		}
	}

	return 0;
}