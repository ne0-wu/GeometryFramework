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
	for (auto v : grid.vertices())
		for (auto u : grid.vertices())
		{
			std::cout << "pos of v: " << grid.point(v).transpose() << std::endl;
			std::cout << "pos of u: " << grid.point(u).transpose() << std::endl;

			GeodesicPath geodesic_path(grid, v, u);
			std::cout << "On-edge distance:  " << geodesic_path.on_edge_distance() << std::endl;
			std::cout << "Geodesic distance: " << geodesic_path.geodesic_distance() << std::endl;

			std::cout << "Real distance:     "
					  << (grid.point(v) - grid.point(u)).norm() << std::endl;
		}

	return 0;
}