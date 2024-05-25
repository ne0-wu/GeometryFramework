#include <set>

#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include "GeometryProcessing.h"

#include "Utils/TickTock.h"

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

	Eigen::SparseMatrix<double> L(num_vertices_original, num_vertices_original);
	L.setFromTriplets(triplets.begin(), triplets.end());

	return L;
}

Eigen::DiagonalMatrix<double, Eigen::Dynamic> SpectralSimplification::diagonal_mass_matrix(const Mesh &mesh)
{
	OpenMesh::FProp<double> area_f(mesh);
	for (auto f : mesh.faces())
		area_f[f] = mesh.calc_sector_area(mesh.halfedge_handle(f));

	Eigen::VectorXd mass = Eigen::VectorXd::Zero(num_vertices_original);

	for (auto v : mesh.vertices())
	{
		double area_v = 0;
		for (auto f : mesh.vf_range(v))
			area_v += area_f[f];
		mass[v.idx()] = area_v / 3;
	}

	Eigen::DiagonalMatrix<double, Eigen::Dynamic> M(mass);

	return M;
}

Eigen::SparseMatrix<double> SpectralSimplification::local_laplacian_matrix(const Mesh &mesh, std::vector<Mesh::VertexHandle> query_vertices)
{
	std::set<OpenMesh::SmartHalfedgeHandle> halfedges;
	for (auto v : query_vertices)
		for (auto h : mesh.voh_range(v))
		{
			halfedges.insert(h);
			halfedges.insert(h.opp());
		}

	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(halfedges.size() * 4);

	for (auto h : halfedges)
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

	Eigen::SparseMatrix<double> L(num_vertices_original, num_vertices_original);
	L.setFromTriplets(triplets.begin(), triplets.end());

	return L;
}

SpectralSimplification::SpectralSimplification(Mesh &input_mesh, int k)
	: mesh(input_mesh), num_vertices_original(input_mesh.n_vertices()), k(k),
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
	Z = M.inverse() * L * F;

	// Initialize the metric
	for (auto v : mesh.vertices())
		E[v] = 0.0;

	// Compute the metric for each edge
	TickTock tt("calculate cost on all edges");
	for (auto e : mesh.edges())
		if (mesh.is_collapse_ok(e.halfedge()))
			SPMs[e] = calc_cost(e.halfedge());

	tt.tock();

	P = Eigen::MatrixXd::Identity(num_vertices_original, num_vertices_original);
}

SpectralSimplification::SPM SpectralSimplification::calc_cost(Mesh::HalfedgeHandle h)
{
	Mesh::VertexHandle fr = mesh.from_vertex_handle(h),
					   to = mesh.to_vertex_handle(h);

	// Build the restriction matrix Q
	Eigen::SparseMatrix<double> Q(num_vertices_original, num_vertices_original);
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(num_vertices_original);
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
	mesh_after.garbage_collection();

	// Extract the One-Ring of the edge
	std::set<Mesh::VertexHandle> one_ring;
	for (auto v : mesh_before.vv_range(fr))
		one_ring.insert(v);
	for (auto v : mesh_before.vv_range(to))
		one_ring.insert(v);
	one_ring.erase(fr);
	std::vector<Mesh::VertexHandle> H;
	H.reserve(one_ring.size());
	for (auto v : one_ring)
		H.push_back(v);

	// Compute E of vertices in H
	auto L_tilde = local_laplacian_matrix(mesh_after, H);
	auto M_tilde = diagonal_mass_matrix(mesh_after);
	Eigen::MatrixXd X = Q * Z - M_tilde.inverse() * L_tilde * Q * F;

	std::vector<double> E_H_after(H.size());
	for (int i = 0; i < H.size(); ++i)
		E_H_after[i] = pow(X.row(H[i].idx()).norm(), 2) * M_tilde.diagonal()[H[i].idx()];

	// Compute the cost
	double cost = 0;
	for (int i = 0; i < H.size(); ++i)
		cost += E_H_after[i] - E[H[i]];
	cost -= E[fr];

	return {Q, cost, H, E_H_after};
}

void SpectralSimplification::collapse_edge()
{
	Mesh::EdgeHandle min_edge;
	double min_cost = std::numeric_limits<double>::infinity();
	for (auto e : mesh.edges())
	{
		if (!mesh.is_collapse_ok(e.halfedge()))
			continue;

		if (SPMs[e].cost < min_cost)
		{
			min_cost = SPMs[e].cost;
			min_edge = e;
		}
	}

	Mesh::VertexHandle to = mesh.to_vertex_handle(mesh.halfedge_handle(min_edge)),
					   fr = mesh.from_vertex_handle(mesh.halfedge_handle(min_edge));
	mesh.set_point(to, (mesh.point(fr) + mesh.point(to)) / 2);

	mesh.collapse(mesh.halfedge_handle(min_edge));
	mesh.garbage_collection();

	// update affected vertices
	for (int i = 0; i < SPMs[min_edge].H.size(); ++i)
		E[SPMs[min_edge].H[i]] = SPMs[min_edge].E_H[i];

	// update affected edges
	for (auto voh : mesh.voh_range(to))
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
	TickTock tt("collapse 1 edge");
	for (int i = 0; i < num_vertices_original - target_num_vertices; i++)
	{
		tt.tick();
		collapse_edge();
		tt.tock();
	}
}