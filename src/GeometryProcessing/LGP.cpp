// This file implements the Local-Global Parameterization algorithm.

#include <algorithm>

#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/SparseLU>

#include "GeometryProcessing.h"

struct Triplet : public Eigen::Triplet<double>
{
	Triplet(int i, int j, double v) : Eigen::Triplet<double>(i, j, v) {}
};

double oppositeAngle(const Mesh &mesh, Mesh::HalfedgeHandle heh)
{
	auto v0 = mesh.from_vertex_handle(heh);
	auto v1 = mesh.to_vertex_handle(heh);
	auto v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));

	auto p0 = mesh.point(v0);
	auto p1 = mesh.point(v1);
	auto p2 = mesh.point(v2);

	auto e0 = p1 - p2;
	auto e1 = p0 - p2;

	return acos((e0.normalized()).dot(e1.normalized()));
}

Eigen::SparseMatrix<double> cotLaplacian(const Mesh &mesh)
{
	int numVertices = mesh.numVertices();

	std::vector<Triplet> triplets;
	triplets.reserve(numVertices * 7);

	// auto oppositeAngle = [&mesh](Mesh::HalfedgeHandle heh)
	// {
	// 	auto v0 = mesh.from_vertex_handle(heh);
	// 	auto v1 = mesh.to_vertex_handle(heh);
	// 	auto v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));
	// 	auto p0 = mesh.point(v0);
	// 	auto p1 = mesh.point(v1);
	// 	auto p2 = mesh.point(v2);
	// 	auto e0 = p1 - p2;
	// 	auto e1 = p0 - p2;
	// 	return acos((e0.normalized()).dot(e1.normalized()));
	// };

	for (auto vi : mesh.vertices())
	{
		int i = vi.idx();
		double sumCotAlphaBeta = 0.0;

		for (auto voh : mesh.voh_range(vi))
		{
			auto vj = mesh.to_vertex_handle(voh);
			int j = vj.idx();

			double cotAlpha = 0.0, cotBeta = 0.0;
			if (voh.is_valid() && !mesh.is_boundary(voh))
				cotAlpha = 1.0 / tan(oppositeAngle(mesh, voh));

			if (mesh.opposite_halfedge_handle(voh).is_valid() && !mesh.opposite_halfedge_handle(voh).is_boundary())
				cotBeta = 1.0 / tan(oppositeAngle(mesh, mesh.opposite_halfedge_handle(voh)));

			triplets.push_back(Triplet(i, j, cotAlpha + cotBeta));

			sumCotAlphaBeta += cotAlpha + cotBeta;
		}

		triplets.push_back(Triplet(i, i, -sumCotAlphaBeta));
	}

	Eigen::SparseMatrix<double> L(numVertices, numVertices);
	L.setFromTriplets(triplets.begin(), triplets.end());

	return L;
}

std::vector<int> boundaryVertices(const Mesh &mesh)
{
	std::vector<int> boundaryVertices;

	OpenMesh::SmartHalfedgeHandle firstBoundaryHalfedge;

	for (auto heh : mesh.halfedges())
		if (mesh.is_boundary(heh))
		{
			firstBoundaryHalfedge = heh;
			break;
		}

	// If you are on a boundary, the next halfedge is guaranteed to be also a boundary halfedge.
	for (auto heh = firstBoundaryHalfedge; heh.is_valid(); heh = mesh.next_halfedge_handle(heh))
	{
		boundaryVertices.push_back(mesh.to_vertex_handle(heh).idx());
		if (mesh.to_vertex_handle(heh) == mesh.from_vertex_handle(firstBoundaryHalfedge))
			break;
	}

	return boundaryVertices;
}

struct SVD22
{
	// For an input matrix A, the SVD decomposition is given by A = U * S * V^T

	Eigen::Matrix2d U;
	Eigen::Matrix2d V;
	Eigen::Matrix2d S;

	SVD22(Eigen::Matrix2d &A)
	{
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
};

Eigen::Matrix2d bestFitARAP(Eigen::Matrix2d &J)
{
	SVD22 svd(J);
	return svd.U * svd.V.transpose();
}

Eigen::Matrix2d bestFitASAP(Eigen::Matrix2d &J)
{
	SVD22 svd(J);
	double meanS = (svd.S(0, 0) + svd.S(1, 1)) / 2;
	Eigen::Matrix2d S;
	S << meanS, 0.0, 0.0, meanS;
	return svd.U * S * svd.V.transpose();
}

// Tutte Parameterization with given laplacian matrix
Eigen::MatrixX2d tutteParameterization(const Mesh &mesh, Eigen::SparseMatrix<double> laplacian)
{
	int numVertices = mesh.numVertices();
	std::vector<int> boundary = boundaryVertices(mesh);

	// Fix the boundary vertices to the unit circle
	Eigen::MatrixX2d b(numVertices + boundary.size(), 2);
	b.setZero();
	for (int i = 0; i < boundary.size(); i++)
	{
		b(numVertices + i, 0) = cos(2 * M_PI * i / boundary.size());
		b(numVertices + i, 1) = sin(2 * M_PI * i / boundary.size());
	}

	// Construct the matrix A for the linear system Ax = b
	// Start by extracting triplets from the laplacian matrix
	std::vector<Triplet> triplets;
	triplets.reserve(laplacian.nonZeros());

	for (int k = 0; k < laplacian.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(laplacian, k); it; ++it)
		{
			// Skip the rows of boundary vertices
			if (std::find(boundary.begin(), boundary.end(), it.row()) != boundary.end())
				continue;

			// Add the triplet to the list
			triplets.push_back(Triplet(it.row(), it.col(), it.value()));
		}

	// Add the boundary constraints
	for (int i = 0; i < boundary.size(); i++)
		triplets.push_back(Triplet(numVertices + i, boundary[i], 1.0));

	// Construct the sparse matrix A
	Eigen::SparseMatrix<double> A(numVertices + boundary.size(), numVertices);
	A.setFromTriplets(triplets.begin(), triplets.end());

	// Solve the linear system
	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr(A);
	return qr.solve(b);
}

Eigen::MatrixX2d localGlobalParameterization(Mesh &mesh, int numIter)
{
	Eigen::SparseMatrix<double> laplacian = cotLaplacian(mesh);
	Eigen::SimplicialLDLT solver(laplacian);

	// Initial guess
	Eigen::MatrixX2d U = tutteParameterization(mesh, laplacian);

	// Directly flattern the faces to the 2D plane
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

	// Pre compute the cotangent weights
	std::vector<double> cotWeights(mesh.numHalfEdges(), 0.0);
	for (auto heh : mesh.halfedges())
	{
		if (heh.is_boundary() || !heh.is_valid())
			continue;
		cotWeights[heh.idx()] = 1.0 / tan(oppositeAngle(mesh, heh));
	}

	for (int iter = 0; iter < numIter; iter++)
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

			Eigen::Vector2d u0 = U.row(v0.idx()),
							u1 = U.row(v1.idx()),
							u2 = U.row(v2.idx());

			auto e01 = u1 - u0, e02 = u2 - u0;

			Eigen::Matrix2d triangleU;
			triangleU << e01(0), e01(1),
				e02(0), e02(1);

			// Compute the best fit rotation from the Jacobi matrix
			// x and u are row vectors, so x * J == u, J = x^-1 * u
			Eigen::Matrix2d J = triangleXs[f.idx()].inverse() * triangleU;
			Eigen::Matrix2d R = bestFitARAP(J);

			// Compute the right hand side of the linear system
			Eigen::Matrix<double, 1, 2> x0, x1, x2;
			x0 << 0, 0;
			x1 << triangleXs[f.idx()](0, 0), triangleXs[f.idx()](0, 1);
			x2 << triangleXs[f.idx()](1, 0), triangleXs[f.idx()](1, 1);

			auto heh12 = mesh.next_halfedge_handle(heh01);
			auto heh20 = mesh.next_halfedge_handle(heh12);

			rhs.block<1, 2>(v0.idx(), 0) += cotWeights[heh01.idx()] * ((x0 - x1) * R) + cotWeights[heh20.idx()] * ((x0 - x2) * R);
			rhs.block<1, 2>(v1.idx(), 0) += cotWeights[heh12.idx()] * ((x1 - x2) * R) + cotWeights[heh01.idx()] * ((x1 - x0) * R);
			rhs.block<1, 2>(v2.idx(), 0) += cotWeights[heh20.idx()] * ((x2 - x0) * R) + cotWeights[heh12.idx()] * ((x2 - x1) * R);
		}

		// Global step
		U = solver.solve(rhs);
	}

	return U;
}