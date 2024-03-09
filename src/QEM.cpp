// #include <queue>
#include <map>
#include <algorithm>
#include <random>

#include "Mesh.h"

void Mesh::collapseEdge(Mesh::HalfedgeHandle edge, Eigen::Vector3d contractedPosition)
{
	VertexHandle v0 = to_vertex_handle(edge);
	VertexHandle v1 = from_vertex_handle(edge);

	set_point(v0, Mesh::Point(contractedPosition.x(), contractedPosition.y(), contractedPosition.z()));
	// set_point(v1, Mesh::Point(contractedPosition.x(), contractedPosition.y(), contractedPosition.z()));

	collapse(edge);

	// std::cout << "edge collapsed" << std::endl;

	garbage_collection();

	if (v0.is_valid())
	{
		std::cout << "v0 is valid" << std::endl;
	}
	else
	{
		std::cout << "v0 is not valid" << std::endl;
	}
	if (v1.is_valid())
	{
		std::cout << "v1 is valid" << std::endl;
	}
	else
	{
		std::cout << "v1 is not valid" << std::endl;
	}
}

void Mesh::collapseEdge(Mesh::HalfedgeHandle edge, Point contractedPosition)
{
	VertexHandle v0 = to_vertex_handle(edge);
	VertexHandle v1 = from_vertex_handle(edge);

	set_point(v0, contractedPosition);
	set_point(v1, contractedPosition);

	collapse(edge);

	this->garbage_collection();

	if (v0.is_valid())
	{
		std::cout << "v0 is valid" << std::endl;
		std::cout << "v0 idx: " << v0.idx() << std::endl;
	}
	else
	{
		std::cout << "v0 is not valid" << std::endl;
	}
	if (v1.is_valid())
	{
		std::cout << "v1 is valid" << std::endl;
		std::cout << "v1 idx: " << v1.idx() << std::endl;
	}
	else
	{
		std::cout << "v1 is not valid" << std::endl;
	}
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

	// std::cout << "print Q:" << std::endl;
	// std::cout << Q << std::endl;

	return Q;
}

// Mesh::Point Mesh::optimalPlacement(const Eigen::Matrix4d &Q, HalfedgeHandle edge)
// {
// 	Eigen::Vector3d p0 = eigenPoint(to_vertex_handle(edge));
// 	Eigen::Vector3d p1 = eigenPoint(from_vertex_handle(edge));
// 	Eigen::Vector4d p0Homogeneous(p0.x(), p0.y(), p0.z(), 1);
// 	Eigen::Vector4d p1Homogeneous(p1.x(), p1.y(), p1.z(), 1);
// 	double a = (p0Homogeneous.transpose() * Q * p0Homogeneous).value();
// 	double b = (p1Homogeneous.transpose() * Q * p1Homogeneous).value();
// 	return a / (a + b) * point(to_vertex_handle(edge)) + b / (a + b) * point(from_vertex_handle(edge));
// }

// double Mesh::quadricErrorEdge(const Eigen::Matrix4d &Q, HalfedgeHandle edge)
// {
// 	Eigen::Vector3d p0 = eigenPoint(to_vertex_handle(edge));
// 	Eigen::Vector3d p1 = eigenPoint(from_vertex_handle(edge));
// 	Eigen::Vector4d p0Homogeneous(p0.x(), p0.y(), p0.z(), 1);
// 	Eigen::Vector4d p1Homogeneous(p1.x(), p1.y(), p1.z(), 1);
// 	return (p0Homogeneous.transpose() * Q * p0Homogeneous + p1Homogeneous.transpose() * Q * p1Homogeneous).value();
// }

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

	// Eigen::Matrix3d A = Q.block<3, 3>(0, 0);
	// Eigen::Vector3d midPoint = (eigenPoint(to_vertex_handle(edge)) + eigenPoint(from_vertex_handle(edge))) / 2;
	// Eigen::Vector4d midPointHomogeneous(midPoint.x(), midPoint.y(), midPoint.z(), 1);
	// double errorMid = (midPointHomogeneous.transpose() * Q * midPointHomogeneous).value();

	// if (abs(A.determinant()) <= 1e-8)
	// {
	// 	return errorMid;
	// }
	// else
	// {
	// 	Eigen::Vector3d b = -Q.block<3, 1>(0, 3);
	// 	Eigen::Vector3d x = A.inverse() * b;
	// 	Eigen::Vector4d xHomogeneous(x.x(), x.y(), x.z(), 1);
	// 	return std::min((xHomogeneous.transpose() * Q * xHomogeneous).value(), errorMid);
	// }

	// Eigen::Vector3d p0 = eigenPoint(to_vertex_handle(edge));
	// Eigen::Vector3d p1 = eigenPoint(from_vertex_handle(edge));
	// auto p = (p0 + p1) / 2;
	// Eigen::Vector4d pHomogeneous(p.x(), p.y(), p.z(), 1);
	// return (pHomogeneous.transpose() * Q * pHomogeneous).value();
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
		return qem > other.qem;
	}
};

template <typename T>
class PriorityQueue
{
public:
	std::vector<T> data;

	PriorityQueue() {}

	void push(const T &value)
	{
		data.push_back(value);
		std::push_heap(data.begin(), data.end());
	}

	void pop()
	{
		std::pop_heap(data.begin(), data.end());
		data.pop_back();
	}

	void update(int i, const T &value)
	{
		data[i] = value;
		std::make_heap(data.begin(), data.end());
	}

	const T &top() const { return data.front(); }

	const T &operator[](int i) const { return data[i]; }

	bool empty() const { return data.empty(); }

	size_t size() const { return data.size(); }

	void heapify() { std::make_heap(data.begin(), data.end()); }
};

bool Mesh::simplifyQEM2(int targetNumVertices)
{
	request_face_normals();
	request_edge_status();
	request_halfedge_status();
	request_vertex_status();
	request_face_status();

	for (auto e_it = edges_begin(); targetNumVertices > 0; ++e_it)
	{
		if (e_it == edges_end())
			e_it = edges_begin();
		if (std::rand() < (targetNumVertices / (float)numVertices()) / 10)
		{
			VertexHandle v0 = to_vertex_handle(halfedge_handle(*e_it, 0));
			VertexHandle v1 = to_vertex_handle(halfedge_handle(*e_it, 1));
			Eigen::Vector3d contractedPosition = (eigenPoint(v0) + eigenPoint(v1)) / 2;
			collapseEdge(halfedge_handle(*e_it, 0), contractedPosition);
			targetNumVertices--;
		}
	}
	return true;
}

bool Mesh::simplifyQEM(int targetNumVertices)
{
	request_face_normals();
	request_edge_status();
	request_vertex_status();
	request_face_status();

	for (int i = numVertices(); i > targetNumVertices; i--)
	{
		std::cout << "remaining vertices: " << numVertices() << std::endl;
		std::map<VertexHandle, Eigen::Matrix4d> Qs;
		for (auto vertex = vertices_begin(); vertex != vertices_end(); ++vertex)
			Qs[vertex] = quadricErrorMatrix(vertex);
		std::cout << "calculated Qs" << std::endl;
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
				std::cout << "there exists a valid edge" << std::endl;
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
		std::cout << "min qem: " << minQEM << std::endl;
		if (minQEM > 1e-3)
		{
			std::cout << "min qem too large" << std::endl;
			return false;
		}
		collapseEdge(minHalfedge, optimalPlacement(Qs[to_vertex_handle(minHalfedge)] + Qs[from_vertex_handle(minHalfedge)], minHalfedge));
	}
	return true;
}

// bool Mesh::simplifyQEM(int targetNumVertices)
// {
// 	request_face_normals();
// 	request_edge_status();
// 	request_halfedge_status();
// 	request_vertex_status();
// 	request_face_status();

// 	// for (int i = numVertices(); i > targetNumVertices; i--)
// 	// {
// 	// 	std::priority_queue<QEM> pq;
// 	// 	for (auto e_it = edges_begin(); e_it != edges_end(); ++e_it)
// 	// 	{
// 	// 		VertexHandle v0 = to_vertex_handle(halfedge_handle(*e_it, 0));
// 	// 		VertexHandle v1 = to_vertex_handle(halfedge_handle(*e_it, 1));
// 	// 		Eigen::Matrix4d Q = quadricErrorMatrix(v0) + quadricErrorMatrix(v1);
// 	// 		Eigen::Vector4d p0 = eigenPoint(v0).homogeneous();
// 	// 		Eigen::Vector4d p1 = eigenPoint(v1).homogeneous();
// 	// 		// double qem = (p0.transpose() * Q * p0 + p1.transpose() * Q * p1).value();
// 	// 		double qem = (((p0 + p1) / 2).transpose() * Q * (p0 + p1) / 2).value();
// 	// 		Eigen::Vector3d contractedPosition = (eigenPoint(v0) + eigenPoint(v1)) / 2;
// 	// 		pq.push(QEM(halfedge_handle(*e_it, 0), qem, contractedPosition));
// 	// 	}
// 	// 	while (!is_collapse_ok(pq.top().edge))
// 	// 		pq.pop();
// 	// 	QEM qem = pq.top();
// 	// 	pq.pop();
// 	// 	collapseEdge(qem.edge, qem.contractedPosition);
// 	// }

// 	PriorityQueue<QEM> pq;

// 	std::map<VertexHandle, Eigen::Matrix4d> Qs;

// 	for (auto v_it = vertices_begin(); v_it != vertices_end(); ++v_it)
// 	{
// 		Qs[*v_it] = quadricErrorMatrix(*v_it);
// 	}

// 	for (auto e_it = edges_begin(); e_it != edges_end(); ++e_it)
// 	{
// 		VertexHandle v0 = to_vertex_handle(halfedge_handle(*e_it, 0));
// 		VertexHandle v1 = from_vertex_handle(halfedge_handle(*e_it, 0));

// 		Eigen::Matrix4d Q = Qs[v0] + Qs[v1];
// 		Eigen::Vector4d p0 = eigenPoint(v0).homogeneous();
// 		Eigen::Vector4d p1 = eigenPoint(v1).homogeneous();

// 		double qem = (((p0 + p1) / 2).transpose() * Q * (p0 + p1) / 2).value();
// 		Eigen::Vector3d contractedPosition = (eigenPoint(v0) + eigenPoint(v1)) / 2;

// 		pq.push(QEM(halfedge_handle(*e_it, 0), qem, contractedPosition));
// 	}

// 	for (int i = numVertices(); i > targetNumVertices; i--)
// 	{
// 		while (!is_collapse_ok(pq.top().edge))
// 			pq.pop();
// 		VertexHandle remainingVertex = to_vertex_handle(pq.top().edge);

// 		collapseEdge(pq.top().edge, pq.top().contractedPosition);
// 		pq.pop();

// 		// update Qs
// 		Qs[remainingVertex] = quadricErrorMatrix(remainingVertex);
// 		for (auto voh_it = voh_begin(remainingVertex); voh_it != voh_end(remainingVertex); ++voh_it)
// 		{
// 			VertexHandle v0 = to_vertex_handle(*voh_it);
// 			Qs[v0] = quadricErrorMatrix(v0);
// 		}

// 		// update pq
// 		for (auto voh_it = voh_begin(remainingVertex); voh_it != voh_end(remainingVertex); ++voh_it)
// 		{
// 			VertexHandle v0 = to_vertex_handle(*voh_it);
// 			for (auto voh_it2 = voh_begin(v0); voh_it2 != voh_end(v0); ++voh_it2)
// 			{
// 				HalfedgeHandle edge = *voh_it2;
// 				for (int i = 0; i < pq.size(); i++)
// 				{
// 					if (pq[i].edge == opposite_halfedge_handle(edge))
// 					{
// 						edge = opposite_halfedge_handle(edge);
// 					}
// 					if (pq[i].edge == edge)
// 					{
// 						VertexHandle v1 = to_vertex_handle(edge);
// 						Eigen::Matrix4d Q = Qs[v0] + Qs[v1];
// 						Eigen::Vector4d p0 = eigenPoint(v0).homogeneous();
// 						Eigen::Vector4d p1 = eigenPoint(v1).homogeneous();
// 						double qem = (((p0 + p1) / 2).transpose() * Q * (p0 + p1) / 2).value();
// 						Eigen::Vector3d contractedPosition = (eigenPoint(v0) + eigenPoint(v1)) / 2;
// 						pq.update(i, QEM(edge, qem, contractedPosition));
// 						break;
// 					}
// 				}
// 			}
// 		}
// 	}

// 	return true;
// }