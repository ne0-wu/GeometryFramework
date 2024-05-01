#include <queue>

#include "GeometryProcessing.h"

GeodesicPath::GeodesicPath(Mesh const &input_mesh, Mesh::VertexHandle source, Mesh::VertexHandle target)
	: original_mesh(input_mesh), intrinsic_mesh(input_mesh), src(source), tgt(target),
	  edge_length(intrinsic_mesh), direction(intrinsic_mesh), angle_sum(intrinsic_mesh),
	  orig_direction(original_mesh)
{
	// Initialize edge length and direction properties
	for (auto e : intrinsic_mesh.edges())
		edge_length[e] = intrinsic_mesh.calc_edge_length(e);
	for (auto v : intrinsic_mesh.vertices())
	{
		double angle_sum_temp = 0.0;
		for (auto h : intrinsic_mesh.voh_ccw_range(v)) // out-going halfedges around v in counter-clockwise order
		{
			direction[h] = angle_sum_temp;
			angle_sum_temp += intrinsic_mesh.calc_dihedral_angle(h);
		}
		for (auto h : intrinsic_mesh.voh_ccw_range(v))
		{
			direction[h] *= 2 * M_PI / angle_sum_temp;
			orig_direction[h] = direction[h];
		}
		angle_sum[v] = angle_sum_temp;
	}

	dijkstra(); // Initialize the path
}

// Returns the angle at the *from* vertex of the halfedge
double GeodesicPath::angle(Mesh::HalfedgeHandle h) const
{
	double l0 = edge_length[intrinsic_mesh.edge_handle(h)],
		   l1 = edge_length[intrinsic_mesh.edge_handle(intrinsic_mesh.next_halfedge_handle(h))],
		   l2 = edge_length[intrinsic_mesh.edge_handle(intrinsic_mesh.prev_halfedge_handle(h))];
	return acos((l0 * l0 + l2 * l2 - l1 * l1) / (2 * l0 * l2));
}

void GeodesicPath::dijkstra()
{
	OpenMesh::VProp<double> dist(intrinsic_mesh);
	OpenMesh::VProp<Mesh::HalfedgeHandle> prev(intrinsic_mesh);

	for (auto v : intrinsic_mesh.vertices())
		dist[v] = std::numeric_limits<double>::infinity();
	dist[src] = 0;

	std::priority_queue<Mesh::VertexHandle, std::vector<Mesh::VertexHandle>,
						std::function<bool(Mesh::VertexHandle, Mesh::VertexHandle)>>
		pq([&dist](Mesh::VertexHandle v0, Mesh::VertexHandle v1)
		   { return dist[v0] > dist[v1]; });

	pq.push(src);

	while (!pq.empty())
	{
		auto u = pq.top();
		pq.pop();

		for (auto h : intrinsic_mesh.voh_range(u))
		{
			auto v = h.to();
			double alt_dist = dist[u] + edge_length[h.edge()];
			if (alt_dist < dist[v])
			{
				dist[v] = alt_dist;
				prev[v] = h;
				pq.push(v);
			}
		}
	}

	// Reconstruct the path
	assert(dist[tgt] < std::numeric_limits<double>::infinity());
	path.clear();
	for (auto v = tgt; v != src; v = intrinsic_mesh.from_vertex_handle(prev[v]))
		path.push_back(v);
	path.push_back(src);
	std::reverse(path.begin(), path.end());

	dijkstra_distance = dist[tgt];
}

void GeodesicPath::intrinsic_flip(Mesh::EdgeHandle e)
{
	auto h0 = intrinsic_mesh.smart_halfedge_handle(e, 0);
	auto h1 = intrinsic_mesh.smart_halfedge_handle(e, 1);

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
	auto p1 = Eigen::Vector2d(l01, 0);
	auto p2 = l02 * Eigen::Vector2d(cos(ang0_12), -sin(ang0_12));
	auto p3 = l03 * Eigen::Vector2d(cos(ang0_13), sin(ang0_13));

	// Check if the flip is valid, i.e. the new edge crosses v0--v1
	double x0 = (p2[1] * p3[0] - p2[0] * p3[1]) / (p2[1] - p3[1]);
	assert(x0 > 0 && x0 < l01);

	// Compute values for the new edge
	double l23 = (p2 - p3).norm();
	double ang2_03 = acos((l02 * l02 + l23 * l23 - l03 * l03) / (2 * l02 * l23));
	double ang3_12 = acos((l13 * l13 + l23 * l23 - l12 * l12) / (2 * l13 * l23));

	// Update the mesh (edge length and signpost)
	auto h02 = h0.next();
	auto h13 = h1.next();
	intrinsic_mesh.flip(e);
	edge_length[e] = l23;
	auto h23 = h02.next();
	auto h32 = h13.next();
	assert(h23.edge() == e && h32.edge() == e);

	direction[h23] = fmod(direction[h02.opp()] - ang2_03 / angle_sum[v2] * 2 * M_PI, 2 * M_PI);
	direction[h32] = fmod(direction[h13.opp()] - ang3_12 / angle_sum[v3] * 2 * M_PI, 2 * M_PI);
}

// Flips out the left-hand-side wedge of the two halfedges
void GeodesicPath::flip_out(int i, bool left_hand_side)
{
	// Decide which side of the wedge to flip out
	// ------------------------------------------
	OpenMesh::SmartHalfedgeHandle h0, h1;
	if (left_hand_side)
	{
		h0 = intrinsic_mesh.find_halfedge(path[i - 1], path[i]);
		h1 = intrinsic_mesh.find_halfedge(path[i], path[i + 1]);
	}
	else
	{
		h0 = intrinsic_mesh.find_halfedge(path[i + 1], path[i]);
		h1 = intrinsic_mesh.find_halfedge(path[i], path[i - 1]);
	}

	// Find the first edge with beta < pi to flip
	// ------------------------------------------
	bool should_flip = true;
	double beta = 0;

	while (should_flip)
	{
		should_flip = false;
		auto edge_to_flip = h0.next().edge();

		for (auto h = h0.next().opp(); h.opp() != h1; h = h.next().opp())
		{
			beta = angle(h) + angle(h.opp().next());
			if (beta < M_PI - eps)
			{
				should_flip = true;
				edge_to_flip = intrinsic_mesh.edge_handle(h);
				break;
			}
		}

		if (should_flip)
			intrinsic_flip(edge_to_flip);
	}

	// Update the path
	// ---------------
	// Remove the vertex i in the old path
	path.erase(path.begin() + i, path.begin() + i + 1);
	// Construct the new path
	std::vector<Mesh::VertexHandle> new_path;
	for (auto h = h0.next().opp(); h.opp() != h1; h = h.next().opp())
		new_path.push_back(h.from());
	// Insert the new path into the old path
	if (left_hand_side)
		path.insert(path.begin() + i, new_path.begin(), new_path.end());
	else
		path.insert(path.begin() + i, new_path.rbegin(), new_path.rend());
}

void GeodesicPath::find_geodesic_path()
{
	// Find the wedge with the smallest angle < pi, and flip it out
	double min_wedge_angle = 0;
	while (min_wedge_angle < std::numeric_limits<double>::infinity())
	{
		min_wedge_angle = std::numeric_limits<double>::infinity();

		// Find the smallest wedge angle
		int min_idx = 0;
		bool left_hand_side = true; // which wedge to flip out
		for (int i = 1; i < path.size() - 1; i++)
		{
			auto h0 = intrinsic_mesh.find_halfedge(path[i - 1], path[i]);
			auto h1 = intrinsic_mesh.find_halfedge(path[i], path[i + 1]);

			double wedge_angle_rght = 0,
				   wedge_angle_left = 0;

			for (auto h = h0; h.opp() != h1; h = h.next().opp())
			{
				if (h.is_boundary()) // Flip-out cannot be performed if the wedge contains a boundary
				{
					wedge_angle_left = std::numeric_limits<double>::infinity();
					break;
				}
				wedge_angle_left += angle(h.next());
			}

			for (auto h = h1.opp(); h != h0; h = h.next().opp())
			{
				if (h.is_boundary())
				{
					wedge_angle_rght = std::numeric_limits<double>::infinity();
					break;
				}
				wedge_angle_rght += angle(h.next());
			}

			// If all wedge angles > pi,
			// then min_wedge_angle == infinity, and the loop is terminated
			if (wedge_angle_rght > M_PI - eps && wedge_angle_left > M_PI - eps)
				continue;

			if (std::min(wedge_angle_rght, wedge_angle_left) < min_wedge_angle)
			{
				min_wedge_angle = std::min(wedge_angle_rght, wedge_angle_left);
				min_idx = i;
				left_hand_side = wedge_angle_left < wedge_angle_rght;
			}
		}

		// Flip the smallest wedge
		if (min_wedge_angle < M_PI)
			flip_out(min_idx, left_hand_side);
	}
}

double GeodesicPath::geodesic_distance()
{
	find_geodesic_path();

	double dist = 0;
	for (int i = 0; i < path.size() - 1; i++)
		dist += edge_length[intrinsic_mesh.find_halfedge(path[i], path[i + 1]).edge()];
	return dist;
}

std::vector<Mesh::Point> GeodesicPath::geodesic_path()
{
	std::vector<Mesh::Point> path_points;

	for (int i = 0; i < path.size() - 1; i++)
	{
		// Trace the path v0->v1 on the original mesh
		auto v0 = path[i], v1 = path[i + 1];
		auto h = intrinsic_mesh.find_halfedge(v0, v1);

		auto p0 = intrinsic_mesh.point(v0), p1 = intrinsic_mesh.point(v1);
		path_points.push_back(p0);

		// Check if v0 and v1 are connected in the original mesh
		if (original_mesh.find_halfedge(v0, v1).is_valid())
		{
			path_points.push_back(p1);
			continue;
		}

		// Start tracing
		enum TraceType
		{
			TRACE_FROM_VERTEX,
			TRACE_FROM_EDGE
		} trace_type = TRACE_FROM_VERTEX;

		Mesh::VertexHandle start_vertex = v0;
		Mesh::HalfedgeHandle start_edge;
		double pos_on_edge;
		double dir = direction[h];
		bool has_reached_v1 = false;

		while (!has_reached_v1)
		{
			switch (trace_type)
			{
			case TRACE_FROM_VERTEX:
			{
				// Find the intersected halfedge
				OpenMesh::SmartHalfedgeHandle intersected_halfedge;
				for (auto voh : original_mesh.voh_ccw_range(start_vertex))
				{
					double dir_voh = orig_direction[voh],
						   dir_voh_prev = orig_direction[voh.prev().opp()];
					// dir_voh_prev should be the next of dir_voh in the ccw order
					if (dir_voh > dir_voh_prev)
					{
						if (dir >= dir_voh || dir < dir_voh_prev)
						{
							intersected_halfedge = voh;
							break;
						}
					}
					else if (dir_voh <= dir && dir < dir_voh_prev)
					{
						intersected_halfedge = voh.next();
						break;
					}
				}

				// Find the intersection point
				// Let p0, p1, p2 (ccw) be the three vertices of the triangle,
				// so that p0 is the start_vertex, and p1--p2 is the intersected edge,
				// and p3 be the intersection point.
				double ang0_13 = fmod(dir - orig_direction[intersected_halfedge.prev()], 2 * M_PI);
				double ang1_30 = angle(intersected_halfedge);
				double ang3_01 = M_PI - ang0_13 - ang1_30;
				double l13 = original_mesh.calc_edge_length(intersected_halfedge.prev().edge()) * sin(ang0_13) / sin(ang3_01);

				pos_on_edge = l13 / original_mesh.calc_edge_length(intersected_halfedge.prev().edge());
				assert(pos_on_edge <= 1);

				// Prepare for the next step

				if (pos_on_edge < eps)
				{
					start_vertex = intersected_halfedge.from();
					dir = fmod(orig_direction[intersected_halfedge.prev().opp()] + M_PI, 2 * M_PI);
				}
				else if (pos_on_edge > 1 - eps)
				{
					start_vertex = intersected_halfedge.to();
					dir = fmod(orig_direction[intersected_halfedge.next()] + M_PI, 2 * M_PI);
				}
				else
				{
					trace_type = TRACE_FROM_EDGE;

					dir = ang0_13 + ang1_30;
					start_edge = intersected_halfedge.opp();
				}
			}
			break;
			case TRACE_FROM_EDGE:
			{
				// TODO: finish this part
			}
			break;
			}
		}
	}

	return path_points;
}