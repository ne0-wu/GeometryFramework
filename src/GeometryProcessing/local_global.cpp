// This file contains the implementation of the Local-Global Parameterization algorithm.

#include <Eigen/Sparse>

#include "GeometryProcessing.h"
#include "Utils/TickTock.h"

#ifdef IMPLEMENT_LGP

typedef Eigen::Triplet<double> Triplet;

struct SVD22
{
	// For an input matrix A, the SVD decomposition is given by A = U * S * V^T
	Eigen::Matrix2d U;
	Eigen::Matrix2d V;
	Eigen::Matrix2d S;

	SVD22(Eigen::Matrix2d &A)
	{
		// This algorithm actually provides the signed SVD decomposition
		// where sigma 2 could be negative
		double E = (A(0, 0) + A(1, 1)) / 2,
			   F = (A(0, 0) - A(1, 1)) / 2,
			   G = (A(1, 0) + A(0, 1)) / 2,
			   H = (A(1, 0) - A(0, 1)) / 2;

		double Q = sqrt(E * E + H * H), R = sqrt(F * F + G * G);

		double S1 = Q + R, S2 = Q - R;

		double T1 = atan2(G, F), T2 = atan2(H, E);

		double theta = (T2 - T1) / 2, phi = (T2 + T1) / 2;

		U = Eigen::Rotation2Dd(phi).toRotationMatrix();
		V = Eigen::Rotation2Dd(theta).toRotationMatrix().transpose();
		S << S1, 0.0, 0.0, S2;
	}

	Eigen::Matrix2d bestFitARAP()
	{
		return U * V.transpose();
	}

	Eigen::Matrix2d bestFitASAP()
	{
		double meanS = (S(0, 0) + S(1, 1)) / 2;
		Eigen::Matrix2d S;
		S << meanS, 0.0, 0.0, meanS;
		return U * S * V.transpose();
	}
};

// Eigen::Matrix2d bestFitARAP(Eigen::Matrix2d &J)
// {
// 	SVD22 svd(J);
// 	return svd.U * svd.V.transpose();
// }

// Eigen::Matrix2d bestFitASAP(Eigen::Matrix2d &J)
// {
// 	SVD22 svd(J);
// 	double meanS = (svd.S(0, 0) + svd.S(1, 1)) / 2;
// 	Eigen::Matrix2d S;
// 	S << meanS, 0.0, 0.0, meanS;
// 	return svd.U * S * svd.V.transpose();
// }

double Tutte::oppositeAngle(Mesh::HalfedgeHandle heh) const
{
	auto p0 = mesh.point(mesh.from_vertex_handle(heh));
	auto p1 = mesh.point(mesh.to_vertex_handle(heh));
	auto p2 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));

	return acos(((p1 - p2).normalized()).dot((p0 - p2).normalized()));
}

void Tutte::computeCotangents()
{
	cotangents.resize(mesh.numHalfEdges(), 0.0);

	for (auto heh : mesh.halfedges())
	{
		if (heh.is_boundary() || !heh.is_valid())
			continue;
		cotangents[heh.idx()] = 1.0 / tan(oppositeAngle(heh));
	}
}

void Tutte::findBoundaryVertices()
{
	boundaryVertices.clear();

	Mesh::HalfedgeHandle firstBoundaryHalfedge;

	for (auto heh : mesh.halfedges())
		if (mesh.is_boundary(heh))
		{
			firstBoundaryHalfedge = heh;
			break;
		}

	// OpenMesh Doc: If you are on a boundary, the next halfedge is guaranteed to be also a boundary halfedge.
	for (auto heh = firstBoundaryHalfedge; heh.is_valid(); heh = mesh.next_halfedge_handle(heh))
	{
		boundaryVertices.push_back(mesh.to_vertex_handle(heh));
		if (mesh.to_vertex_handle(heh) == mesh.from_vertex_handle(firstBoundaryHalfedge))
			break;
	}
}

Eigen::SparseMatrix<double> Tutte::laplacianUniform()
{
	std::vector<Triplet> triplets;
	triplets.reserve(mesh.numHalfEdges() + mesh.numVertices());

	for (auto v : mesh.vertices())
	{
		double valence = mesh.valence(v);
		for (auto vv : mesh.vv_range(v))
			triplets.push_back(Triplet(v.idx(), vv.idx(), -1.0 / valence));
		triplets.push_back(Triplet(v.idx(), v.idx(), 1.0));
	}

	Eigen::SparseMatrix<double> L(mesh.numVertices(), mesh.numVertices());
	L.setFromTriplets(triplets.begin(), triplets.end());

	return L;
}

Eigen::SparseMatrix<double> Tutte::laplacianCotangent()
{
	std::vector<Triplet> triplets;
	triplets.reserve(mesh.numHalfEdges() * 4);

	for (auto heh : mesh.halfedges())
	{
		if (heh.is_boundary() || !heh.is_valid())
			continue;

		double cot = cotangents[heh.idx()];
		triplets.push_back(Triplet(heh.from().idx(), heh.to().idx(), cot));
		triplets.push_back(Triplet(heh.to().idx(), heh.from().idx(), cot));
		triplets.push_back(Triplet(heh.from().idx(), heh.from().idx(), -cot));
		triplets.push_back(Triplet(heh.to().idx(), heh.to().idx(), -cot));
	}

	Eigen::SparseMatrix<double> L(mesh.numVertices(), mesh.numVertices());
	L.setFromTriplets(triplets.begin(), triplets.end());

	return L;
}

void Tutte::tutte()
{
	findBoundaryVertices();

	switch (laplacianType)
	{
	case LaplacianType::UNIFORM:
		laplacian = laplacianUniform();
		break;
	case LaplacianType::COTANGENT:
		computeCotangents();
		laplacian = laplacianCotangent();
		break;
	}

	int numVertices = mesh.numVertices();

	// Fix the boundary vertices to the unit circle
	Eigen::MatrixX2d b(numVertices, 2);
	b.setZero();
	for (int i = 0; i < boundaryVertices.size(); i++)
	{
		int idx_vi = boundaryVertices[i].idx();
		b(idx_vi, 0) = cos(2 * M_PI * i / boundaryVertices.size());
		b(idx_vi, 1) = sin(2 * M_PI * i / boundaryVertices.size());
	}

	// Construct the matrix A for the linear system Ax = b
	std::vector<Triplet> triplets;
	triplets.reserve(laplacian.nonZeros());

	// Add the boundary constraints
	for (int i = 0; i < boundaryVertices.size(); i++)
		triplets.push_back(Triplet(boundaryVertices[i].idx(), boundaryVertices[i].idx(), 1.0));

	// Extract triplets from the laplacian matrix
	std::vector<int> boundaryVertexIndices(boundaryVertices.size());
	for (int i = 0; i < boundaryVertices.size(); i++)
		boundaryVertexIndices[i] = boundaryVertices[i].idx();
	std::sort(boundaryVertexIndices.begin(), boundaryVertexIndices.end());
	for (int k = 0; k < laplacian.outerSize(); ++k)
	{
		// Skip the rows of boundary vertices
		// Assume that the matrix is symmetric
		if (std::binary_search(boundaryVertexIndices.begin(), boundaryVertexIndices.end(), k))
			continue;

		for (Eigen::SparseMatrix<double>::InnerIterator it(laplacian, k); it; ++it)
		{
			// Add the triplet to the list
			if (laplacian.Flags & Eigen::RowMajorBit)
				triplets.push_back(Triplet(it.row(), it.col(), it.value())); // Row major
			else
				triplets.push_back(Triplet(it.col(), it.row(), it.value())); // Column major
		}
	}

	// Construct the sparse matrix A
	Eigen::SparseMatrix<double> A(numVertices, numVertices);
	A.setFromTriplets(triplets.begin(), triplets.end());

	// Solve the linear system
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);

	uv = solver.solve(b);
}

void Tutte::flatten()
{
	tutte();
}

void LocalGlobal::localGlobal()
{
	// Initial guess
	setLaplacianType(LaplacianType::COTANGENT);
	tutte(); // also computes the cotangents and the laplacian

	// Precompute the solver
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(laplacian);

	// Directly flatten the faces to the 2D plane
	std::vector<Eigen::Matrix2d> triangleXs(mesh.numFaces());
	for (auto f : mesh.faces())
	{
		auto heh = mesh.halfedge_handle(f);
		auto v0 = mesh.from_vertex_handle(heh),
			 v1 = mesh.to_vertex_handle(heh),
			 v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));

		auto x0 = mesh.point(v0),
			 x1 = mesh.point(v1),
			 x2 = mesh.point(v2);

		Eigen::Vector2d e01, e02;
		double theta = acos((x1 - x0).normalized().dot((x2 - x0).normalized()));
		e01 << (x1 - x0).norm(), 0;
		e02 = (x2 - x0).norm() * Eigen::Vector2d(cos(theta), sin(theta));

		// Use row vectors
		triangleXs[f.idx()] << e01(0), e01(1),
			e02(0), e02(1);
	}

	for (int iter = 1; iter <= numIter; iter++)
	{
		// Local step
		Eigen::MatrixX2d rhs(mesh.numVertices(), 2);
		rhs.setZero();

		for (auto f : mesh.faces())
		{
			// Get the 2D triangle in the UV plane, from the last iteration
			auto heh01 = mesh.halfedge_handle(f);
			auto v0 = mesh.from_vertex_handle(heh01),
				 v1 = mesh.to_vertex_handle(heh01),
				 v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh01));

			Eigen::Vector2d u0 = uv.row(v0.idx()), u1 = uv.row(v1.idx()), u2 = uv.row(v2.idx());

			auto e01 = u1 - u0, e02 = u2 - u0;

			Eigen::Matrix2d triangleU;
			triangleU << e01(0), e01(1),
				e02(0), e02(1);

			// Compute the best fit rotation from the Jacobi matrix
			// x and u are row vectors, so x * J == u, J = x^-1 * u
			Eigen::Matrix2d J = triangleXs[f.idx()].inverse() * triangleU;

			Eigen::Matrix2d R;
			switch (target)
			{
			case LocalGlobalTarget::ARAP:
				R = SVD22(J).bestFitARAP();
				break;
			case LocalGlobalTarget::ASAP:
				R = SVD22(J).bestFitASAP();
				break;
			}

			// Compute the right hand side of the linear system
			Eigen::Matrix<double, 1, 2> x0, x1, x2;
			x0 << 0, 0;
			x1 << triangleXs[f.idx()](0, 0), triangleXs[f.idx()](0, 1);
			x2 << triangleXs[f.idx()](1, 0), triangleXs[f.idx()](1, 1);

			auto heh12 = mesh.next_halfedge_handle(heh01);
			auto heh20 = mesh.next_halfedge_handle(heh12);

			rhs.block<1, 2>(v0.idx(), 0) += cotangents[heh01.idx()] * ((x0 - x1) * R) + cotangents[heh20.idx()] * ((x0 - x2) * R);
			rhs.block<1, 2>(v1.idx(), 0) += cotangents[heh12.idx()] * ((x1 - x2) * R) + cotangents[heh01.idx()] * ((x1 - x0) * R);
			rhs.block<1, 2>(v2.idx(), 0) += cotangents[heh20.idx()] * ((x2 - x0) * R) + cotangents[heh12.idx()] * ((x2 - x1) * R);
		}

		// Global step
		uv = solver.solve(rhs);
	}
}

void LocalGlobal::flatten()
{
	localGlobal();
}

#endif