#include <set>

#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include "GeometryProcessing.h"

Eigen::SparseMatrix<double> SpectralSimplification::laplacian_matrix(const Mesh &mesh)
{
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(mesh.n_halfedges() * 4);

	for (auto h : mesh.halfedges())
	{
		if (h.is_boundary())
			continue;

		auto fr = mesh.point(h.from()),
			 to = mesh.point(h.to()),
			 op = mesh.point(h.next().to());
		double cot_op = 1.0 / tan(acos((fr - op).normalized().dot((to - op).normalized())));

		triplets.push_back(Eigen::Triplet<double>(h.to().idx(), h.to().idx(), cot_op / 2));
		triplets.push_back(Eigen::Triplet<double>(h.from().idx(), h.from().idx(), cot_op / 2));
		triplets.push_back(Eigen::Triplet<double>(h.to().idx(), h.from().idx(), -cot_op / 2));
		triplets.push_back(Eigen::Triplet<double>(h.from().idx(), h.to().idx(), -cot_op / 2));
	}

	Eigen::SparseMatrix<double> L(mesh.n_vertices(), mesh.n_vertices());
	L.setFromTriplets(triplets.begin(), triplets.end());

	return L;
}

Eigen::SparseMatrix<double> SpectralSimplification::diagonal_mass_matrix(const Mesh &mesh)
{
	OpenMesh::FProp<double> area(mesh);
	for (auto f : mesh.faces())
		area[f] = mesh.calc_sector_area(mesh.halfedge_handle(f));

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(mesh.n_vertices());

	for (auto v : mesh.vertices())
	{
		double area_v = 0;
		for (auto vf : mesh.vf_range(v))
			area_v += area[vf];
		triplets.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), area_v / 3));
	}

	Eigen::SparseMatrix<double> M(mesh.n_vertices(), mesh.n_vertices());
	M.setFromTriplets(triplets.begin(), triplets.end());

	return M;
}

Eigen::SparseMatrix<double> SpectralSimplification::diagonal_mass_matrix_inv(const Mesh &mesh)
{
	OpenMesh::FProp<double> area(mesh);
	for (auto f : mesh.faces())
		area[f] = mesh.calc_sector_area(mesh.halfedge_handle(f));

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(mesh.n_vertices());

	for (auto v : mesh.vertices())
	{
		double area_v = 0.0;
		for (auto vf : mesh.vf_range(v))
			area_v += area[vf];
		triplets.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), 3 / area_v));
	}

	Eigen::SparseMatrix<double> M(mesh.n_vertices(), mesh.n_vertices());
	M.setFromTriplets(triplets.begin(), triplets.end());

	return M;
}

SpectralSimplification::SpectralSimplification(Mesh &input_mesh, int k)
	: mesh(input_mesh), k(k),
	  F(mesh.n_vertices(), k),
	  E(mesh), SPMs(mesh)
{
	L = laplacian_matrix(mesh);
	M = diagonal_mass_matrix(mesh);

	// Compute F
	Spectra::SparseSymMatProd<double> op(L);
	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, k, 2 * k + 1);
	eigs.init();
	eigs.compute();
	F = eigs.eigenvectors(); // F is a matrix of size |V| x k

	// Compute Z = M^-1 * L * F
	Eigen::SparseMatrix<double> M_inv = diagonal_mass_matrix_inv(mesh);
	Z = M_inv * L * F;

	// Initialize the metric
	for (auto v : mesh.vertices())
		E[v] = 0.0;

	for (auto e : mesh.edges())
	{
		if (!mesh.is_collapse_ok(e.halfedge()))
			continue;

		SPMs[e] = calc_cost(e.halfedge());
	}
}

SpectralSimplification::SPM SpectralSimplification::calc_cost(Mesh::HalfedgeHandle h)
{
	Mesh::VertexHandle fr = mesh.from_vertex_handle(h),
					   to = mesh.to_vertex_handle(h);

	// Build the restriction matrix Q
	Eigen::SparseMatrix<double> Q(mesh.n_vertices(), mesh.n_vertices());
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(mesh.n_vertices());
	for (auto v : mesh.vertices())
		if (v != fr && v != to)
			triplets.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0));
	triplets.push_back(Eigen::Triplet<double>(to.idx(), to.idx(), 0.5));
	triplets.push_back(Eigen::Triplet<double>(to.idx(), fr.idx(), 0.5));
	Q.setFromTriplets(triplets.begin(), triplets.end());

	// A copy of the mesh after the edge collapse
	Mesh &mesh_before = mesh;
	Mesh mesh_after = mesh_before;
	mesh_after.set_point(to, (mesh_before.point(fr) + mesh_before.point(to)) / 2);
	mesh_after.collapse(h);

	// Extract the One-Ring of the edge
	std::set<Mesh::VertexHandle> one_ring;
	for (auto v : mesh.vv_range(fr))
		one_ring.insert(v);
	for (auto v : mesh.vv_range(to))
		one_ring.insert(v);
	one_ring.erase(fr);
	std::vector<Mesh::VertexHandle> H;
	H.reserve(one_ring.size());
	for (auto v : one_ring)
		H.push_back(v);

	// Compute E of vertices in H
	Eigen::SparseMatrix<double> L_tilde = laplacian_matrix(mesh_after);
	Eigen::SparseMatrix<double> M_tilde_inv = diagonal_mass_matrix_inv(mesh_after);
	Eigen::MatrixXd X = Q * Z - M_tilde_inv * L_tilde * Q * F;

	std::vector<double> E_H_after(H.size());
	for (int i = 0; i < H.size(); ++i)
		E_H_after[i] = pow(X.row(H[i].idx()).norm(), 2) / M_tilde_inv.coeff(H[i].idx(), H[i].idx());

	std::cout << X.rows() << " " << X.cols() << std::endl;
	std::cout << X.row(to.idx()).norm() << std::endl;

	// Compute the cost
	double cost = 0;
	for (int i = 0; i < H.size(); ++i)
		cost += E_H_after[i] - E[H[i]];

	return {Q, cost, H, E_H_after};
}

void SpectralSimplification::collapse_edge()
{
	Mesh::EdgeHandle min_edge;
	double min_cost = std::numeric_limits<double>::infinity();
	for (auto e : mesh.edges())
	{
		if (mesh.is_collapse_ok(e.halfedge()))
			continue;

		if (SPMs[e].cost < min_cost)
		{
			min_cost = SPMs[e].cost;
			min_edge = e;
		}
	}

	Mesh::VertexHandle remaining_vertex = mesh.to_vertex_handle(mesh.halfedge_handle(min_edge));
	mesh.set_point(remaining_vertex, (mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(min_edge))) + mesh.point(remaining_vertex)) / 2);

	mesh.collapse(mesh.halfedge_handle(min_edge));

	// update affected vertices
	for (auto v : SPMs[min_edge].H)
		E[v] = SPMs[min_edge].E_H[v.idx()];

	// update affected edges
	for (auto voh : mesh.voh_range(remaining_vertex))
		for (auto he : mesh.voh_range(voh.to()))
			if (mesh.is_collapse_ok(he))
				SPMs[he.edge()] = calc_cost(he);

	// update matrices
	P = SPMs[min_edge].Q * P;
	F = SPMs[min_edge].Q * F;
	Z = SPMs[min_edge].Q * Z;
}

void SpectralSimplification::simplify(int target_num_vertices)
{
	for (int i = 0; i < mesh.n_vertices() - target_num_vertices; i++)
		collapse_edge();
}

// SpectralSimplification::SPM::SPM(Mesh mesh, Mesh::HalfedgeHandle h)
// 	: Q(mesh.n_vertices(), mesh.n_vertices())
// {
// 	Mesh::VertexHandle fr = mesh.from_vertex_handle(h),
// 					   to = mesh.to_vertex_handle(h);

// 	// Build the restriction matrix Q
// 	std::vector<Eigen::Triplet<double>> triplets;
// 	triplets.reserve(mesh.n_vertices() + 1);
// 	for (auto v : mesh.vertices())
// 		if (v != fr && v != to)
// 			triplets.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1.0));
// 	triplets.push_back(Eigen::Triplet<double>(to.idx(), to.idx(), 0.5));
// 	triplets.push_back(Eigen::Triplet<double>(to.idx(), fr.idx(), 0.5));
// 	Q.setFromTriplets(triplets.begin(), triplets.end());

// 	// A copy of the mesh after the edge collapse
// 	Mesh mesh_copy = mesh;
// 	mesh_copy.set_point(to, (mesh.point(fr) + mesh.point(to)) / 2);
// 	mesh_copy.collapse(h);

// 	// Extract the One-Ring of the edge
// 	std::set<Mesh::VertexHandle> one_ring;
// }