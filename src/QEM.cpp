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

struct edgeWithQEM
{
	typedef Mesh::HalfedgeHandle HalfedgeHandle;
	typedef Mesh::Point Point;

	HalfedgeHandle edge;
	double qem;
	Point contractedPosition;

	edgeWithQEM(HalfedgeHandle edge, OptimalPlacement optimalPlacement) : edge(edge)
	{
		qem = optimalPlacement.first;
		contractedPosition = optimalPlacement.second;
	}

	bool operator<(const edgeWithQEM &other) const
	{
		return qem < other.qem;
	}
};

bool Mesh::simplifyQEM(int targetNumVertices)
{
	request_face_normals();
	update_face_normals();

	request_vertex_status();
	request_edge_status();
	request_face_status();

	for (int i = numVertices(); i > targetNumVertices; i--)
	{
		std::map<VertexHandle, Eigen::Matrix4d> Qs;
		for (auto vertex : vertices())
			Qs[vertex] = quadricErrorMatrix(vertex);

		double minQEM = std::numeric_limits<double>::max();
		HalfedgeHandle minHalfedge = *halfedges().begin();
		Point minContractedPosition;

		bool existsValidEdge = false;
		for (auto edge : halfedges())
		{
			if (!is_collapse_ok(edge))
				continue;
			if (existsValidEdge == false)
				existsValidEdge = true;

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
		if (!existsValidEdge)
		{
			std::cout << "no valid edge" << std::endl;
			return false;
		}
		if (minQEM > 1e-3)
		{
			std::cout << "min qem too large: " << minQEM << std::endl;
			return false;
		}
		collapseEdge(minHalfedge, minContractedPosition);
	}
	return true;
}