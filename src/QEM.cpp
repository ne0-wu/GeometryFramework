#include <queue>

#include <Eigen/Core>

#include "Mesh.h"

bool Mesh::removeVertex()
{
	// TODO:
	return false;
}

bool Mesh::collapseEdge()
{
	// TODO:
	return false;
}

Eigen::Matrix4d Mesh::quadricErrorMatrix(VertexHandle v)
{
	Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();
	Eigen::Vector3d p = point(v).data();

	for (auto vf_it = vf_begin(v); vf_it.is_valid(); ++vf_it)
	{
		Eigen::Vector3d normal = this->normal(vf_it).data();
		double d = -normal.dot(p);
		Eigen::Vector4d plane(normal.x(), normal.y(), normal.z(), d);
		Q += plane * plane.transpose();
	}

	return Q;
}

struct QEM
{
	typedef OpenMesh::VertexHandle VertexHandle;

	double qem;
	VertexHandle v1, v2;
	bool isEdge;

	QEM(VertexHandle v, double qem) : v1(v), qem(qem), isEdge(false) {}

	QEM(VertexHandle v1, VertexHandle v2, double qem) : v1(v1), v2(v2), qem(qem), isEdge(true) {}

	bool operator<(const QEM &other) const
	{
		return qem < other.qem;
	}
};

bool Mesh::simplifyQEM(int targetNumVertices)
{
	std::priority_queue<QEM> pq;

	for (auto v : vertices())
	{
		Eigen::Vector3d y(point(v).data());
		for (auto vf_it = vf_begin(v); vf_it.is_valid(); ++vf_it)
		{
			Eigen::Vector3d normal = normal(vf_it);
		}

		// pq.push(QEM(v, 1));
	}

	return false;
}