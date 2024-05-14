// This file contains the implementation of the Quadric Error Metrics (QEM) simplification algorithm.

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include "Mesh.h"
#include "GeometryProcessing.h"

#ifdef IMPLEMENT_QEM_SIMPLIFICATION

Eigen::Matrix4d QEMSimplification::quadricErrorMatrix(Mesh::VertexHandle v)
{
	Eigen::Matrix4d Q = Eigen::Matrix4d::Zero();
	Eigen::Vector3d p = mesh.point(v);

	for (auto f : mesh.vf_range(v))
	{
		Eigen::Vector3d normal = mesh.normal(f);
		double n = normal.norm();
		double d = (normal.transpose() * p).value();
		Eigen::Vector4d plane(normal.x(), normal.y(), normal.z(), -d);
		Q += plane * plane.transpose();
	}

	return Q;
}

QEMSimplification::QEM QEMSimplification::calc_cost(Mesh::HalfedgeHandle edge, Eigen::Matrix4d Q)
{
	Eigen::Matrix4d A = Q;
	A.block<1, 3>(3, 0) = Eigen::Vector3d(0, 0, 0).transpose();
	A(3, 3) = 1;

	Eigen::Vector3d mid = (mesh.point(mesh.to_vertex_handle(edge)) + mesh.point(mesh.from_vertex_handle(edge))) / 2;
	Eigen::Vector4d midHomo(mid.x(), mid.y(), mid.z(), 1);
	double errorMid = (midHomo.transpose() * Q * midHomo).value();

	if (abs(A.determinant()) <= 1e-10)
		return {edge, mid, errorMid};

	Eigen::Vector4d x = A.llt().solve(Eigen::Vector4d(0, 0, 0, 1));
	x = x / x[3];
	double error2 = (x.transpose() * Q * x).value();
	if (error2 < errorMid)
		return {edge, Mesh::Point(x.x(), x.y(), x.z()), error2};
	else
		return {edge, mid, errorMid};
}

QEMSimplification::QEMSimplification(Mesh &input_mesh) : mesh(input_mesh), Qs(mesh), QEMs(mesh)
{
	// Needs vertex/edge/face status attribute in order to delete the items that degenerate.
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_halfedge_status();
	mesh.request_face_status();

	// Needs face normals to compute the quadric error matrices
	mesh.request_face_normals();
	mesh.update_face_normals();

	// Compute quadric error matrices
	for (auto v : mesh.vertices())
		Qs[v] = quadricErrorMatrix(v);

	// Compute quadric error on each edge
	for (auto eh : mesh.edges())
	{
		if (!mesh.is_collapse_ok(eh.halfedge()))
			continue;

		QEMs[eh] = calc_cost(eh.halfedge(), Qs[eh.halfedge().to()] + Qs[eh.halfedge().from()]);
	}
}

void QEMSimplification::collapse1Edge()
{
	QEM minQem;
	for (auto &eh : mesh.edges())
	{
		if (!mesh.is_collapse_ok(eh.halfedge()))
			continue;

		if (minQem.error > QEMs[eh].error)
			minQem = QEMs[eh];
	}

	auto vh2 = mesh.to_vertex_handle(minQem.heh);
	mesh.set_point(vh2, minQem.optimalPlacement);

	mesh.collapse(minQem.heh); // removes the 'from' vertex of the halfedge

	// update affected faces
	// mesh.update_face_normals();
	for (auto voh : mesh.voh_range(vh2))
		mesh.update_normal(voh.face());

	// update affected vertices
	for (auto voh : mesh.voh_range(vh2))
		Qs[voh.to()] = quadricErrorMatrix(voh.to());

	// update affected halfedges
	for (auto voh : mesh.voh_range(vh2))
		for (auto heh : mesh.voh_range(voh.to()))
			if (mesh.is_collapse_ok(heh))
				QEMs[heh.edge()] = calc_cost(heh, Qs[heh.to()] + Qs[heh.from()]);
}

void QEMSimplification::simplify(int targetNumVertices)
{
	for (int i = 0; i < mesh.n_vertices() - targetNumVertices; i++)
		collapse1Edge();

	mesh.garbage_collection();
}

#endif