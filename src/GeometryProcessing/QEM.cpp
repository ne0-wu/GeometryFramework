// This file implements the Quadric Error Metrics (QEM) simplification algorithm.

#include <vector>
#include <map>
#include <algorithm>

#include "Mesh.h"
#include "GeometryProcessing.h"

#ifdef IMPLEMENT_QEM

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
	Eigen::Vector3d p = point(v);

	for (auto f : vf_range(v))
	{
		auto normal = this->normal(f);
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

	Eigen::Vector3d midPoint = (point(to_vertex_handle(edge)) + point(from_vertex_handle(edge))) / 2;
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
			std::cout << "Number of vertices: " << numVertices() << std::endl;
			std::cout << "qem ended" << std::endl;
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

	std::cout << "qem done" << std::endl;
	std::cout << "Number of vertices: " << numVertices() << std::endl;

	return true;
}

#endif