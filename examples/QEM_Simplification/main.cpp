// #include "Mesh.h"
// #include "GeometryProcessing.h"

// int main()
// {
//     Mesh mesh("meshes/spot_quadrangulated.obj");

//     QEMSimplification qem(mesh);

//     std::cout << "Number of vertices: " << qem.getMesh().numVertices() << std::endl;

//     qem.simplify(500);

//     qem.getMesh().save("output.obj");

//     return 0;
// }

#include <iostream>
#include <filesystem>
#include <vector>

#include "Scene.h"
#include "Mesh.h"
#include "Window.h"

#include "GeometryProcessing.h"

int main()
{
	Window window(1280, 720, "QEM Mesh Simplification");

	Scene scene;
	scene.window = window.window;

	Mesh originalMesh("meshes/spot_quadrangulated.obj");
	originalMesh.fitIntoUnitBall();

	Mesh mesh = originalMesh;
	scene.addMesh(GLMesh(mesh));

	glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
	glfwWindowHint(GLFW_SAMPLES, 4);
	glEnable(GL_DEPTH_TEST);

	// GUI states
	int numVert = originalMesh.n_vertices() / 2;

	while (!window.shouldClose())
	{
		window.pollEvents();
		window.processInput();
		scene.processInput();

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		{
			ImGui::Begin("QEM Mesh Simplification");

			ImGui::SliderInt("Num of Vertices", &numVert, originalMesh.n_vertices() / 5, originalMesh.n_vertices());

			ImGui::End();
		}

		ImGui::Render();

		scene.update();
		scene.draw();

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		window.swapBuffers();

		if (numVert < mesh.n_vertices())
		{
			QEMSimplification qem(mesh);
			qem.simplify(numVert);
			mesh = qem.getMesh();
			scene.glMeshes[0].setMesh(std::make_shared<Mesh>(mesh));
		}
		if (numVert > mesh.n_vertices())
		{
			QEMSimplification qem(originalMesh);
			qem.simplify(numVert);
			mesh = qem.getMesh();
			scene.glMeshes[0].setMesh(std::make_shared<Mesh>(mesh));
		}
	}

	return 0;
}