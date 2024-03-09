#include <vector>
#include <map>
#include <algorithm>

#include "Mesh.h"

void Mesh::collapseEdge(Mesh::HalfedgeHandle edge, Point contractedPosition)
{
	VertexHandle v0 = to_vertex_handle(edge);
	set_point(v0, contractedPosition);

	collapse(edge);

	garbage_collection();
}

Eigen::Matrix4d Mesh::quadricErrorMatrix(VertexHandle v)
{
	Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();
	Eigen::Vector3d p = eigenPoint(v);

	for (auto vf_it = vf_begin(v); vf_it != vf_end(v); ++vf_it)
	{
		Eigen::Vector3d normal = eigenNormal(vf_it);
		double d = (normal.transpose() * p).value();
		Eigen::Vector4d plane(normal.x(), normal.y(), normal.z(), -d);
		Q += plane * plane.transpose();
	}
	return Q;
}

typedef std::pair<double, Mesh::Point> OptimalPlacement;
OptimalPlacement Mesh::optimalPlacement(HalfedgeHandle edge, const Eigen::Matrix4d &Q)
{
	Eigen::Matrix4d A = Q;
	A.block<1, 3>(3, 0) = Eigen::Vector3d(0, 0, 0).transpose();
	A(3, 3) = 1;

	Eigen::Vector3d midPoint = (eigenPoint(to_vertex_handle(edge)) + eigenPoint(from_vertex_handle(edge))) / 2;
	Eigen::Vector4d midPointHomo(midPoint.x(), midPoint.y(), midPoint.z(), 1);
	double errorMid = (midPointHomo.transpose() * Q * midPointHomo).value();

	if (abs(A.determinant()) <= 1e-10)
		return OptimalPlacement(errorMid, Mesh::Point(midPoint.x(), midPoint.y(), midPoint.z()));

	Eigen::Vector4d x = A.llt().solve(Eigen::Vector4d(0, 0, 0, 1));
	x = x / x[3];
	double error2 = (x.transpose() * Q * x).value();
	if (error2 < errorMid)
		return OptimalPlacement(error2, Mesh::Point(x.x(), x.y(), x.z()));
	else
		return OptimalPlacement(errorMid, Mesh::Point(midPoint.x(), midPoint.y(), midPoint.z()));
}

struct EdgeWithQEM
{
	typedef Mesh::HalfedgeHandle HalfedgeHandle;
	typedef Mesh::Point Point;

	HalfedgeHandle edge;
	double qem;
	Point contractedPosition;

	EdgeWithQEM(HalfedgeHandle edge, OptimalPlacement optimalPlacement) : edge(edge)
	{
		qem = optimalPlacement.first;
		contractedPosition = optimalPlacement.second;
	}

	bool operator<(const EdgeWithQEM &other) const
	{
		return qem > other.qem;
	}
};

// bool Mesh::simplifyQEM(int targetNumVertices)
// {
// 	request_face_normals();
// 	update_face_normals();

// 	request_vertex_status();
// 	request_edge_status();
// 	request_face_status();

// 	for (int i = numVertices(); i > targetNumVertices; i--)
// 	{
// 		std::map<VertexHandle, Eigen::Matrix4d> Qs;
// 		for (auto vertex : vertices())
// 			Qs[vertex] = quadricErrorMatrix(vertex);

// 		double minQEM = std::numeric_limits<double>::max();
// 		HalfedgeHandle minHalfedge = *halfedges().begin();
// 		Point minContractedPosition;

// 		for (auto edge : halfedges())
// 		{
// 			if (!is_collapse_ok(edge))
// 				continue;

// 			auto op = optimalPlacement(edge, Qs[to_vertex_handle(edge)] + Qs[from_vertex_handle(edge)]);
// 			double qem = op.first;
// 			Point contractedPosition = op.second;

// 			if (qem < minQEM)
// 			{
// 				minQEM = qem;
// 				minHalfedge = edge;
// 				minContractedPosition = contractedPosition;
// 			}
// 		}
// 		if (minQEM > 1e-3)
// 		{
// 			std::cout << "min qem too large: " << minQEM << std::endl;
// 			return false;
// 		}
// 		collapseEdge(minHalfedge, minContractedPosition);
// 	}
// 	return true;
// }

struct EdgePriorityQueue
{
	std::vector<EdgeWithQEM> edges;

	void push(EdgeWithQEM edgeWithQEM)
	{
		edges.push_back(edgeWithQEM);
		std::push_heap(edges.begin(), edges.end());
	}

	EdgeWithQEM pop()
	{
		std::pop_heap(edges.begin(), edges.end());
		EdgeWithQEM edgeWithQEM = edges.back();
		edges.pop_back();
		return edgeWithQEM;
	}

	// bool updateEdge(Mesh::HalfedgeHandle edge, Mesh::HalfedgeHandle oppositeHe, OptimalPlacement optimalPlacement)
	// {
	// 	for (auto &e : edges)
	// 	{
	// 		if (e.edge == edge || e.edge == oppositeHe)
	// 		{
	// 			e.qem = optimalPlacement.first;
	// 			e.contractedPosition = optimalPlacement.second;
	// 			std::make_heap(edges.begin(), edges.end());
	// 			return true;
	// 		}
	// 	}
	// 	return false;
	// }

	bool removeEdge(Mesh::HalfedgeHandle edge, Mesh::HalfedgeHandle oppositeHe)
	{
		for (auto e_it = edges.begin(); e_it != edges.end(); e_it++)
		{
			if (e_it->edge == edge || e_it->edge == oppositeHe)
			{
				edges.erase(e_it);
				std::make_heap(edges.begin(), edges.end());
				return true;
			}
		}
		return false;
	}

	bool empty()
	{
		return edges.empty();
	}
};

bool Mesh::simplifyQEM(int targetNumVertices)
{
	request_face_normals();
	update_face_normals();

	request_vertex_status();
	request_edge_status();
	request_face_status();

	std::map<VertexHandle, Eigen::Matrix4d> Qs;
	for (auto vertex : vertices())
		Qs[vertex] = quadricErrorMatrix(vertex);

	EdgePriorityQueue edgePQueue;
	for (auto edge : halfedges())
	{
		if (!is_collapse_ok(edge))
			continue;

		OptimalPlacement op = optimalPlacement(edge, Qs[to_vertex_handle(edge)] + Qs[from_vertex_handle(edge)]);
		edgePQueue.push(EdgeWithQEM(edge, op));
	}

	while (numVertices() > targetNumVertices)
	{
		// std::cout << "Number of vertices: " << numVertices() << std::endl;

		double minQEM = std::numeric_limits<double>::max();
		HalfedgeHandle minHalfedge = *halfedges().begin();
		Point minContractedPosition;

		for (auto edge : halfedges())
		{
			if (!is_collapse_ok(edge))
				continue;

			auto op = optimalPlacement(edge, Qs[to_vertex_handle(edge)] + Qs[from_vertex_handle(edge)]);
			double qem = op.first;
			Point contractedPosition = op.second;

			if (qem < minQEM)
			{
				minQEM = qem;
				minHalfedge = edge;
				minContractedPosition = contractedPosition;
			}
		}
		if (minQEM > 1e-3)
		{
			std::cout << "min qem too large: " << minQEM << std::endl;
			return false;
		}

		VertexHandle vTo = to_vertex_handle(minHalfedge);

		collapseEdge(minHalfedge, minContractedPosition);

		for (auto voh_it : voh_range(vTo))
		{
			VertexHandle v = to_vertex_handle(voh_it);
			Qs[v] = quadricErrorMatrix(v);
		}
	}

	while (numVertices() < targetNumVertices)
	{
		std::cout << "Number of vertices: " << numVertices() << std::endl;

		if (edgePQueue.empty())
		{
			std::cout << "Edge priority queue is empty" << std::endl;
			return true;
		}

		EdgeWithQEM minEdge = edgePQueue.pop();

		VertexHandle vFrom = from_vertex_handle(minEdge.edge);

		VertexHandle vTo = to_vertex_handle(minEdge.edge);
		std::cout << "Number of vertices: " << numVertices() << std::endl;
		// VertexHandle vFrom = from_vertex_handle(minEdge.edge);

		if (numVertices() == 2871)
		{
			std::cout << "Edge: " << minEdge.edge.idx() << std::endl;
			std::cout << "is_collapse_ok: " << is_collapse_ok(minEdge.edge) << std::endl;
			std::cout << "QEM: " << minEdge.qem << std::endl;
			std::cout << "Contracted position: " << minEdge.contractedPosition << std::endl;
		}

		for (auto voh_it : voh_range(vTo))
			edgePQueue.removeEdge(voh_it, opposite_halfedge_handle(voh_it));
		for (auto voh_it : voh_range(vFrom))
			edgePQueue.removeEdge(voh_it, opposite_halfedge_handle(voh_it));

		collapseEdge(minEdge.edge, minEdge.contractedPosition);

		std::cout << minEdge.qem << std::endl;

		Qs.erase(vFrom);
		for (auto voh_it : voh_range(vTo))
		{
			VertexHandle v = to_vertex_handle(voh_it);
			Qs[v] = quadricErrorMatrix(v);
		}
		for (auto voh_it : voh_range(vTo))
			edgePQueue.push(EdgeWithQEM(voh_it, optimalPlacement(voh_it, Qs[to_vertex_handle(voh_it)] + Qs[from_vertex_handle(voh_it)])));
	}

	return true;
}