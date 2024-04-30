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
	Mesh::VertexHandle start = mesh.vertex_handle(1);
	Mesh::VertexHandle end = mesh.vertex_handle(100);
	GeodesicPath geodesic_path(mesh, start, end);

	std::cout << "Geodesic distance: " << geodesic_path.geodesic_distance() << std::endl;

	return 0;
}