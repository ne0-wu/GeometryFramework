// #include <queue>
#include <map>
#include <algorithm>
#include <random>

#include "Mesh.h"

void Mesh::collapseEdge(Mesh::HalfedgeHandle edge, Point contractedPosition)
{
	VertexHandle v0 = to_vertex_handle(edge);
	VertexHandle v1 = from_vertex_handle(edge);

	set_point(v0, contractedPosition);
	set_point(v1, contractedPosition);

	collapse(edge);

	this->garbage_collection();
}

Eigen::Matrix4d Mesh::quadricErrorMatrix(VertexHandle v)
{
	Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();
	Eigen::Vector3d p = eigenPoint(v);

	for (auto vf_it = vf_begin(v); vf_it != vf_end(v); ++vf_it)
	{
		if (!vf_it.is_valid())
			continue;
		Eigen::Vector3d normal = eigenNormal(vf_it);
		double d = normal.dot(p);
		Eigen::Vector4d plane(normal.x(), normal.y(), normal.z(), -d);
		auto tmp = (plane * plane.transpose());
		Q += plane * plane.transpose();
	}

	return Q;
}

Mesh::Point Mesh::optimalPlacement(const Eigen::Matrix4d &Q, HalfedgeHandle edge)
{
	Eigen::Matrix4d A = Q;
	A.block<1, 3>(3, 0) = Eigen::Vector3d(0, 0, 0).transpose();
	A(3, 3) = 0;
	Eigen::Vector3d midPoint = (eigenPoint(to_vertex_handle(edge)) + eigenPoint(from_vertex_handle(edge))) / 2;
	Eigen::Vector4d midPointHomogeneous(midPoint.x(), midPoint.y(), midPoint.z(), 1);
	double errorMid = (midPointHomogeneous.transpose() * Q * midPointHomogeneous).value();

	if (abs(A.determinant()) <= 1e-8)
	{
		return (point(to_vertex_handle(edge)) + point(from_vertex_handle(edge))) / 2;
	}
	else
	{
		// Eigen::Vector3d b = -Q.block<3, 1>(0, 3);
		// Eigen::Vector3d x = A.inverse() * b;
		// Eigen::Vector4d xHomogeneous(x.x(), x.y(), x.z(), 1);

		Eigen::Vector4d x = A.inverse() * Eigen::Vector4d(0, 0, 0, 1);
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
		return Mesh::Point(x.x(), x.y(), x.z());
	}
}

double Mesh::quadricErrorEdge(const Eigen::Matrix4d &Q, HalfedgeHandle edge)
{
	Eigen::Matrix4d A = Q;
	A.block<1, 3>(3, 0) = Eigen::Vector3d(0, 0, 0).transpose();
	A(3, 3) = 0;
	Eigen::Vector3d midPoint = (eigenPoint(to_vertex_handle(edge)) + eigenPoint(from_vertex_handle(edge))) / 2;
	Eigen::Vector4d midPointHomogeneous(midPoint.x(), midPoint.y(), midPoint.z(), 1);
	double errorMid = (midPointHomogeneous.transpose() * Q * midPointHomogeneous).value();

	if (abs(A.determinant()) <= 1e-8)
	{
		return errorMid;
	}
	else
	{
		// Eigen::Vector3d b = -Q.block<3, 1>(0, 3);
		// Eigen::Vector3d x = A.inverse() * b;
		// Eigen::Vector4d xHomogeneous(x.x(), x.y(), x.z(), 1);

		Eigen::Vector4d x = A.inverse() * Eigen::Vector4d(0, 0, 0, 1);
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
		return error2;
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
		// std::cout << "remaining vertices: " << numVertices() << std::endl;
		std::map<VertexHandle, Eigen::Matrix4d> Qs;
		for (auto vertex = vertices_begin(); vertex != vertices_end(); ++vertex)
			Qs[vertex] = quadricErrorMatrix(vertex);
		// std::cout << "calculated Qs" << std::endl;
		double minQEM = std::numeric_limits<double>::max();
		HalfedgeHandle minHalfedge = *halfedges().begin();
		bool existsValidEdge = false;
		for (auto edge = halfedges_begin(); edge != halfedges_end(); ++edge)
		{
			if (!is_collapse_ok(edge) || !to_vertex_handle(edge).is_valid() || !from_vertex_handle(edge).is_valid() || is_boundary(edge) || is_boundary(opposite_halfedge_handle(edge)))
				continue;
			if (existsValidEdge == false)
			{
				existsValidEdge = true;
				// std::cout << "there exists a valid edge" << std::endl;
			}
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
		// std::cout << "min qem: " << minQEM << std::endl;
		if (minQEM > 1e-3)
		{
			std::cout << "min qem too large" << std::endl;
			return false;
		}
		collapseEdge(minHalfedge, optimalPlacement(Qs[to_vertex_handle(minHalfedge)] + Qs[from_vertex_handle(minHalfedge)], minHalfedge));
	}
	return true;
}