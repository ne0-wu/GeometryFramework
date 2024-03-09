#include <map>
#include <algorithm>

#include "Mesh.h"

void Mesh::collapseEdge(Mesh::HalfedgeHandle edge, Point contractedPosition)
{
	VertexHandle v0 = to_vertex_handle(edge);
	VertexHandle v1 = from_vertex_handle(edge);

	set_point(v0, contractedPosition);
	set_point(v1, contractedPosition);

	collapse(edge);

	garbage_collection();
	update_normals();
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

Mesh::Point Mesh::optimalPlacement(const Eigen::Matrix4d &Q, HalfedgeHandle edge)
{
	Eigen::Matrix4d A = Q;
	A.block<1, 3>(3, 0) = Eigen::Vector3d(0, 0, 0).transpose();
	A(3, 3) = 1;
	Eigen::Vector3d midPoint = (eigenPoint(to_vertex_handle(edge)) + eigenPoint(from_vertex_handle(edge))) / 2;
	Eigen::Vector4d midPointHomogeneous(midPoint.x(), midPoint.y(), midPoint.z(), 1);
	double errorMid = (midPointHomogeneous.transpose() * Q * midPointHomogeneous).value();

	if (abs(A.determinant()) <= 1e-10)
	{
		return (point(to_vertex_handle(edge)) + point(from_vertex_handle(edge))) / 2;
	}
	else
	{
		Eigen::Vector4d x = A.lu().solve(Eigen::Vector4d(0, 0, 0, 1));
		x = x / x[3];
		double error2 = (x.transpose() * Q * x).value();
		if (error2 < errorMid)
		{
			return Mesh::Point(x.x(), x.y(), x.z());
		}
		else
		{
			return (point(to_vertex_handle(edge)) + point(from_vertex_handle(edge))) / 2;
		}
	}
}

double Mesh::quadricErrorEdge(const Eigen::Matrix4d &Q, HalfedgeHandle edge)
{
	Eigen::Matrix4d A = Q;
	A.block<1, 3>(3, 0) = Eigen::Vector3d(0, 0, 0).transpose();
	A(3, 3) = 1;
	Eigen::Vector3d midPoint = (eigenPoint(to_vertex_handle(edge)) + eigenPoint(from_vertex_handle(edge))) / 2;
	Eigen::Vector4d midPointHomogeneous(midPoint.x(), midPoint.y(), midPoint.z(), 1);
	double errorMid = (midPointHomogeneous.transpose() * Q * midPointHomogeneous).value();

	if (abs(A.determinant()) <= 1e-10)
	{
		return errorMid;
	}
	else
	{
		Eigen::Vector4d x = A.lu().solve(Eigen::Vector4d(0, 0, 0, 1));
		// std::cout << x << std::endl;
		x = x / x[3];
		double error2 = (x.transpose() * Q * x).value();
		if (error2 < errorMid)
		{
			return error2;
		}
		else
		{
			return errorMid;
		}
	}
}

bool Mesh::simplifyQEM(int targetNumVertices)
{
	request_face_normals();
	request_edge_status();
	request_vertex_status();
	request_face_status();

	for (int i = numVertices(); i > targetNumVertices; i--)
	{
		std::map<VertexHandle, Eigen::Matrix4d> Qs;
		for (auto vertex : vertices())
			Qs[vertex] = quadricErrorMatrix(vertex);
		double minQEM = std::numeric_limits<double>::max();
		HalfedgeHandle minHalfedge = *halfedges().begin();
		bool existsValidEdge = false;
		for (auto edge : halfedges())
		{
			if (!is_collapse_ok(edge))
				continue;
			if (existsValidEdge == false)
				existsValidEdge = true;
			VertexHandle v0 = to_vertex_handle(edge);
			VertexHandle v1 = from_vertex_handle(edge);
			Eigen::Matrix4d Q = Qs[v0] + Qs[v1];
			double qem = quadricErrorEdge(Q, edge);
			if (qem < minQEM)
			{
				minQEM = qem;
				minHalfedge = edge;
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
		collapseEdge(minHalfedge, optimalPlacement(Qs[to_vertex_handle(minHalfedge)] + Qs[from_vertex_handle(minHalfedge)], minHalfedge));
	}
	return true;
}