#include <queue>

#include "Mesh.h"

void Mesh::collapseEdge(Mesh::HalfedgeHandle edge, Eigen::Vector3d contractedPosition)
{
	// 1. Get the two vertices of the edge
	VertexHandle v0 = to_vertex_handle(edge);
	VertexHandle v1 = to_vertex_handle(opposite_halfedge_handle(edge));

	// 2. Move the first vertex to the contracted position
	set_point(v0, Mesh::Point(contractedPosition.x(), contractedPosition.y(), contractedPosition.z()));

	// 3. Collapse the edge
	collapse(edge);

	garbage_collection();
}

Eigen::Matrix4d Mesh::quadricErrorMatrix(VertexHandle v)
{
	Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();
	Eigen::Vector3d p = eigenPoint(v);

	for (auto vf_it = vf_begin(v); vf_it.is_valid(); ++vf_it)
	{
		Eigen::Vector3d normal = eigenNormal(vf_it);
		double d = normal.dot(p);
		Eigen::Vector4d plane(normal.x(), normal.y(), normal.z(), -d);
		Q += plane * plane.transpose();
	}

	return Q;
}

struct QEM
{
	typedef OpenMesh::HalfedgeHandle HalfedgeHandle;

	double qem;
	Eigen::Vector3d contractedPosition;
	HalfedgeHandle edge;

	QEM(HalfedgeHandle edge, double qem, Eigen::Vector3d contractedPosition)
		: edge(edge), qem(qem), contractedPosition(contractedPosition) {}

	bool operator<(const QEM &other) const
	{
		return qem < other.qem;
		return true;
	}
};

bool Mesh::simplifyQEM(int targetNumVertices)
{
	request_face_normals();
	request_edge_status();
	request_halfedge_status();
	request_vertex_status();
	request_face_status();

	for (int i = numVertices(); i > targetNumVertices; i--)
	{
		std::priority_queue<QEM> pq;

		for (auto e_it = edges_begin(); e_it != edges_end(); ++e_it)
		{
			if (status(*e_it).deleted())
				continue;

			VertexHandle v0 = to_vertex_handle(halfedge_handle(*e_it, 0));
			VertexHandle v1 = to_vertex_handle(halfedge_handle(*e_it, 1));

			Eigen::Matrix4d Q = quadricErrorMatrix(v0) + quadricErrorMatrix(v1);
			Eigen::Vector4d p0 = eigenPoint(v0).homogeneous();
			Eigen::Vector4d p1 = eigenPoint(v1).homogeneous();

			// double qem = (p0.transpose() * Q * p0 + p1.transpose() * Q * p1).value();
			double qem = (((p0 + p1) / 2).transpose() * Q * (p0 + p1) / 2).value();
			Eigen::Vector3d contractedPosition = ((p0 + p1) / 2).block<3, 1>(0, 0);

			pq.push(QEM(halfedge_handle(*e_it, 0), qem, contractedPosition));
		}

		while (!is_collapse_ok(pq.top().edge))
			pq.pop();

		QEM qem = pq.top();
		pq.pop();

		collapseEdge(qem.edge, qem.contractedPosition);
	}

	garbage_collection();

	return false;
}