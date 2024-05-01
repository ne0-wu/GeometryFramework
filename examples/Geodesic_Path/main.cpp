#include <iostream>
#include <filesystem>
#include <vector>

#include "Renderer/Scene.h"
#include "Renderer/Window.h"

#include "Mesh.h"
#include "GeometryProcessing.h"

int main()
{
	// Test on a grid mesh
	// -------------------
	int n = 5;

	// Build up the grid mesh
	Mesh grid;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			grid.add_vertex(Mesh::Point(i, j, 0));

	for (int i = 0; i < n - 1; ++i)
		for (int j = 0; j < n - 1; ++j)
		{
			auto v0 = grid.vertex_handle(i * n + j);
			auto v1 = grid.vertex_handle((i + 1) * n + j);
			auto v2 = grid.vertex_handle((i + 1) * n + (j + 1));
			auto v3 = grid.vertex_handle(i * n + (j + 1));

			grid.add_face(v0, v1, v2);
			grid.add_face(v0, v2, v3);
		}

	// Test
	if (false)
	{
		int num_test = 0, num_succ = 0;
		for (auto v : grid.vertices())
			for (auto u : grid.vertices())
			{
				std::cout << "pos of v: " << grid.point(v).transpose() << std::endl;
				std::cout << "pos of u: " << grid.point(u).transpose() << std::endl;

				GeodesicPath geodesic_path(grid, v, u);
				double on_edge = geodesic_path.on_edge_distance();
				double geodeisc = geodesic_path.geodesic_distance();
				std::cout << "On-edge distance:  " << on_edge << std::endl;
				std::cout << "Geodesic distance: " << geodeisc << std::endl;

				std::cout << "Real distance:     "
						  << (grid.point(v) - grid.point(u)).norm() << std::endl;

				if (abs(geodeisc - (grid.point(v) - grid.point(u)).norm()) < 1e-6)
					num_succ++;
				num_test++;
			}

		std::cout << "Tested on the grid mesh: " << num_succ << " / " << num_test << " passed\n";
	}

	Mesh original_mesh("meshes/spot_quadrangulated.obj");
	// Mesh original_mesh = grid;
	original_mesh.fitIntoUnitBall();

	Window window(1280, 720, "Geodesic Path");

	Scene scene;
	scene.window = window.window;

	Mesh mesh = original_mesh;
	scene.addMesh(GLMesh(mesh));

	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glfwWindowHint(GLFW_SAMPLES, 4);
	glEnable(GL_DEPTH_TEST);

	struct GUI_State
	{
		int vertex_i = 0;
		int vertex_j = 1;

		double on_edge_distance = 0.0;
		double geodesic_distance = 0.0;

		bool operator!=(const GUI_State &other) const
		{
			return vertex_i != other.vertex_i || vertex_j != other.vertex_j;
		}
	} gui_state_now, gui_state_prev;

	while (!window.shouldClose())
	{
		window.pollEvents();
		window.processInput();
		scene.processInput();

		ImGui_ImplOpenGL3_NewFrame();
		ImGui_ImplGlfw_NewFrame();
		ImGui::NewFrame();

		{
			ImGui::Begin("Geodesic Path");

			ImGui::SliderInt("Vertex i", &gui_state_now.vertex_i,
							 0, mesh.n_vertices() - 1);
			ImGui::SliderInt("Vertex j", &gui_state_now.vertex_j,
							 0, mesh.n_vertices() - 1);

			ImGui::Text("On-edge distance:  %.6f", gui_state_now.on_edge_distance);
			ImGui::Text("Geodesic distance: %.6f", gui_state_now.geodesic_distance);

			ImGui::End();
		}

		ImGui::Render();

		scene.update();
		scene.draw();

		ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

		window.swapBuffers();

		if (gui_state_now != gui_state_prev)
		{
			gui_state_prev = gui_state_now;

			Mesh::VertexHandle v = original_mesh.vertex_handle(gui_state_now.vertex_i);
			Mesh::VertexHandle u = original_mesh.vertex_handle(gui_state_now.vertex_j);
			GeodesicPath geodesic_path(original_mesh, v, u);

			gui_state_now.on_edge_distance = geodesic_path.on_edge_distance();
			gui_state_now.geodesic_distance = geodesic_path.geodesic_distance();

			mesh = geodesic_path.get_intrinsic_mesh();
			scene.glMeshes[0].setMesh(std::make_shared<Mesh>(mesh));
		}
	}

	return 0;
}