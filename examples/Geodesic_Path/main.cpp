#include <iostream>
#include <filesystem>
#include <vector>

#include "Renderer/Scene.h"
#include "Renderer/Window.h"

#include "Mesh.h"
#include "GeometryProcessing.h"

int main()
{
	Mesh mesh("meshes/spot_quadrangulated.obj");

	auto v0 = mesh.vertex_handle(0);

	for (auto v : mesh.vv_range(v0))
	{
		std::cout << v.idx() << ", ";
	}

	for (int i = 0; i < 1; i++)
		for (int j : {1165, 1158, 1159, 812, 813, 767, 768, 764})
		{
			std::cout << "i: " << i << " j: " << j << std::endl;

			Mesh::VertexHandle start = mesh.vertex_handle(i);
			Mesh::VertexHandle end = mesh.vertex_handle(j);
			GeodesicPath geodesic_path(mesh, start, end);

			std::cout << "Geodesic distance: " << geodesic_path.geodesic_distance() << std::endl;
		}

	return 0;
}