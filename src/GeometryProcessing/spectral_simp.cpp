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

SpectralSimplification::SpectralSimplification(Mesh &input_mesh, int k)
	: mesh(input_mesh), k(k)
{
	L = laplacian_matrix(mesh);
	M = diagonal_mass_matrix(mesh);

	// Compute F
	// ---------
	Spectra::SparseSymMatProd<double> op(L);
	Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigs(op, k, 2 * k + 1);
	eigs.init();
	eigs.compute();
	F = eigs.eigenvectors();

	//
}