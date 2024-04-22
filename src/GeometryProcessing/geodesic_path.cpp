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
		double angle_sum_temp = 0.0;
		for (auto h : mesh.voh_ccw_range(v)) // out-going halfedges around v in counter-clockwise order
		{
			direction[h] = angle_sum_temp;
			angle_sum_temp += mesh.calc_dihedral_angle(h);
		}
		for (auto h : mesh.voh_ccw_range(v))
			direction[h] *= 2 * M_PI / angle_sum_temp;
		angle_sum[v] = angle_sum_temp;
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

		for (auto h : mesh.voh_range(u))
		{
			auto v = h.to();
			auto l = edge_length[h.edge()];

			if (distance[u] + l < distance[v])
			{
				distance[v] = distance[u] + l;
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
	double l02 = edge_length[h0.next().edge()];
	double l12 = edge_length[h0.prev().edge()];
	double l13 = edge_length[h1.next().edge()];
	double l03 = edge_length[h1.prev().edge()];
	double ang0_12 = acos((l01 * l01 + l02 * l02 - l12 * l12) / (2 * l01 * l02));
	double ang0_13 = acos((l01 * l01 + l03 * l03 - l13 * l13) / (2 * l01 * l03));

	// Flatten the two triangles
	auto p0 = Eigen::Vector2d(0, 0);
	auto p1 = Eigen::Vector2d(edge_length(e), 0);
	auto p2 = l02 * Eigen::Vector2d(cos(ang0_12), -sin(ang0_12));
	auto p3 = l03 * Eigen::Vector2d(cos(ang0_13), sin(ang0_13));

	// Check if the flip is valid, i.e. the new edge crosses v0--v1
	double x0 = (p2[1] * p3[0] - p2[0] * p3[1]) / (p2[1] - p3[0]);
	assert(x0 > 0 && x0 < l01);

	// Compute values for the new edge
	double l23 = (p2 - p3).norm();
	double ang2_13 = acos((l12 * l12 + l23 * l23 - l13 * l13) / (2 * l12 * l23));
	double ang3_02 = acos((l03 * l03 + l23 * l23 - l02 * l02) / (2 * l03 * l23));

	// Update the mesh
	mesh.flip(e);
	edge_length[e] = l23;
}