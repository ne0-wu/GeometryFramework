#include <queue>

#include "GeometryProcessing.h"

GeodesicPath::GeodesicPath(Mesh const &input_mesh, Mesh::VertexHandle source, Mesh::VertexHandle target)
	: original_mesh(input_mesh), mesh(input_mesh), src(source), tgt(target),
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

	dijkstra(); // Initialize the path
}

// Returns the angle at the *from* vertex of the halfedge
double GeodesicPath::angle(Mesh::HalfedgeHandle h) const
{
	double l0 = edge_length[mesh.edge_handle(h)],
		   l1 = edge_length[mesh.edge_handle(mesh.next_halfedge_handle(h))],
		   l2 = edge_length[mesh.edge_handle(mesh.prev_halfedge_handle(h))];
	return acos((l0 * l0 + l2 * l2 - l1 * l1) / (2 * l0 * l2));
}

void GeodesicPath::dijkstra()
{
	OpenMesh::VProp<double> dist(mesh);
	OpenMesh::VProp<Mesh::HalfedgeHandle> prev(mesh);

	for (auto v : mesh.vertices())
		dist[v] = std::numeric_limits<double>::infinity();
	dist[src] = 0;

	std::vector<Mesh::VertexHandle> reserved_memory;
	std::priority_queue<Mesh::VertexHandle, std::vector<Mesh::VertexHandle>,
						std::function<bool(Mesh::VertexHandle, Mesh::VertexHandle)>>
		pq([&dist](Mesh::VertexHandle v0, Mesh::VertexHandle v1)
		   { return dist[v0] > dist[v1]; },
		   reserved_memory);

	for (auto v : mesh.vertices())
		pq.push(v);

	while (!pq.empty())
	{
		auto u = pq.top();
		pq.pop();

		for (auto h : mesh.voh_range(u))
		{
			auto v = h.to();
			auto w = dist[u] + edge_length[h.edge()];
			if (w < dist[v])
			{
				dist[v] = w;
				prev[v] = h;
			}
		}
	}

	// Reconstruct the path
	assert(dist[tgt] < std::numeric_limits<double>::infinity());
	path.clear();
	for (auto v = tgt; v != src; v = mesh.from_vertex_handle(prev[v]))
		path.push_back(v);
	std::reverse(path.begin(), path.end());

	std::cout << "Length of shortest path on edges: " << dist[tgt] << std::endl;
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
	auto h12 = h1.next();
	auto h03 = h0.next();
	mesh.flip(e);
	edge_length[e] = l23;
	auto h23 = h12.next();
	auto h32 = h03.next();
	assert(h23.edge() == e && h32.edge() == e);

	direction[h23] = direction[h12.opp()] + ang2_13 / angle_sum[v2] * 2 * M_PI;
	direction[h32] = direction[h03.opp()] + ang3_02 / angle_sum[v3] * 2 * M_PI;
}

// Flips out the left-hand-side wedge of the two halfedges
void GeodesicPath::flip_out(int i, bool left_hand_side)
{
	OpenMesh::SmartHalfedgeHandle h0, h1;
	if (left_hand_side)
	{
		h0 = mesh.find_halfedge(path[i - 1], path[i]);
		h1 = mesh.find_halfedge(path[i], path[i + 1]);
	}
	else
	{
		h0 = mesh.find_halfedge(path[i + 1], path[i]);
		h1 = mesh.find_halfedge(path[i], path[i - 1]);
	}

	// Find the first edge with beta < pi to flip
	bool should_flip = true;

	while (should_flip)
	{
		should_flip = false;
		auto edge_to_flip = h0.next().edge();

		for (auto h = h0.next().opp(); h.opp() != h1; h = h.next().opp())
		{
			double beta = angle(h) + angle(h.opp().next());
			if (beta < M_PI)
			{
				should_flip = true;
				edge_to_flip = mesh.edge_handle(h);
				break;
			}
		}

		if (should_flip)
			intrinsic_flip(edge_to_flip);
	}

	// Update the path
	std::vector<Mesh::VertexHandle> new_path;
	for (auto h = h0.next().opp(); h.opp() != h1; h = h.next().opp())
		new_path.push_back(h.from());
	path.erase(path.begin() + i, path.begin() + i + 1);
	if (left_hand_side)
		path.insert(path.begin() + i, new_path.begin(), new_path.end());
	else
		path.insert(path.begin() + i, new_path.rbegin(), new_path.rend());
}

void GeodesicPath::find_geodesic_path()
{
	// std::vector<double> wedge_angles(path.size() - 1);
	double min_wedge_angle = std::numeric_limits<double>::infinity();

	while (min_wedge_angle < std::numeric_limits<double>::infinity())
	{
		// Find the smallest wedge angle
		int min_idx = 0;
		bool left_hand_side = true; // which wedge to flip out
		for (int i = 1; i <= path.size() - 1; i++)
		{
			// auto h0 = path[i];
			// auto h1 = path[i + 1];
			auto h0 = mesh.find_halfedge(path[i - 1], path[i]);
			auto h1 = mesh.find_halfedge(path[i], path[i + 1]);
			double diff_rght = direction[h1] - direction[h0],
				   diff_left = direction[h0] - direction[h1];
			while (diff_rght < 0)
				diff_rght += 2 * M_PI;
			while (diff_rght > 2 * M_PI)
				diff_rght -= 2 * M_PI;
			while (diff_left < 0)
				diff_left += 2 * M_PI;
			while (diff_left > 2 * M_PI)
				diff_left -= 2 * M_PI;
			diff_rght = diff_rght * angle_sum[mesh.to_vertex_handle(h0)] / (2 * M_PI);
			diff_left = diff_left * angle_sum[mesh.to_vertex_handle(h1)] / (2 * M_PI);
			if (diff_rght > M_PI && diff_left > M_PI)
				continue;
			if (std::min(diff_rght, diff_left) < min_wedge_angle)
			{
				min_wedge_angle = std::min(diff_rght, diff_left);
				min_idx = i;
				left_hand_side = diff_left < diff_rght;
			}
		}

		// Flip the smallest wedge
		flip_out(min_idx, left_hand_side);
	}
}