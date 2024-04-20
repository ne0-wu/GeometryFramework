#include "GeometryProcessing.h"

GeodesicPath::GeodesicPath(Mesh const &input_mesh)
	: mesh(input_mesh),
	  edge_length(mesh), angle(mesh)
{
	// Initialize edge length and angle properties
	for (auto e : mesh.edges())
		edge_length[e] = mesh.calc_edge_length(e);
}