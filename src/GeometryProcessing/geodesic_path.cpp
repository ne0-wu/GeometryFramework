#include <queue>

#include "GeometryProcessing.h"

GeodesicPath::GeodesicPath(Mesh const &input_mesh, Mesh::VertexHandle source, Mesh::VertexHandle target)
	: mesh(input_mesh), src(source), tgt(target),
	  edge_length(mesh), direction(mesh), angle_sum(mesh)
{
	// Initialize edge length and direction properties
	for (auto e : mesh.edges())
		edge_length[e] = mesh.calc_edge_length(e);
	for (auto v : mesh.vertices())
	{
		angle_sum[v] = 0;
		for (auto h : mesh.voh_range(v))
		{
			direction[h] = angle_sum[v];
			angle_sum[v] += mesh.calc_dihedral_angle(h);
		}
		for (auto h : mesh.voh_range(v))
			direction[h] *= 2 * M_PI / angle_sum[v];
	}

	dijkstra();
}

void GeodesicPath::dijkstra()
{
	// Initialize the distance and previous vertex properties
	OpenMesh::VProp<double> distance(mesh);
	OpenMesh::VProp<Mesh::VertexHandle> previous(mesh);
	for (auto v : mesh.vertices())
	{
		distance[v] = std::numeric_limits<double>::infinity();
		previous[v] = Mesh::VertexHandle();
	}
	distance[src] = 0;

	// Initialize the priority queue
	std::priority_queue<std::pair<double, Mesh::VertexHandle>> pq;
	pq.push({0, src});

	// Dijkstra's algorithm
	while (!pq.empty())
	{
		auto [d, u] = pq.top();
		pq.pop();

		for (auto h : mesh.vih_range(u))
		{
			auto v = mesh.to_vertex_handle(h);
			auto e = mesh.edge_handle(h);
			auto w = edge_length[e];

			if (distance[u] + w < distance[v])
			{
				distance[v] = distance[u] + w;
				previous[v] = u;
				pq.push({distance[v], v});
			}
		}
	}

	// Extract the path
	path.clear();
	auto v = tgt;
	while (v != src)
	{
		path.push_back(v);
		v = previous[v];
	}
	path.push_back(src);
	std::reverse(path.begin(), path.end());
}

void GeodesicPath::intrinsic_flip(Mesh::EdgeHandle e)
{
	auto h0 = mesh.smart_halfedge_handle(e, 0);
	auto h1 = mesh.smart_halfedge_handle(e, 1);

	// Edge v0--v1 will be flipped to v2--v3
	auto v0 = h0.to();
	auto v1 = h1.to();

	auto v2 = h0.next().to();
	auto v3 = h1.next().to();

	// Edge lengths
	double l01 = edge_length[e];

	// Calculate the new edge length, on the 2D plane
	auto p0 = Eigen::Vector2d(0, 0);
	auto p1 = Eigen::Vector2d(edge_length(e), 0);
}