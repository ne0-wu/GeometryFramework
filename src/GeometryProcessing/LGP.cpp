// This file implements the Local-Global Parameterization algorithm.

#include <Eigen/Sparse>

#include "GeometryProcessing.h"

Eigen::SparseMatrix<double> cotLaplacian(const Mesh &mesh)
{
	int numVertices = mesh.numVertices();

	struct Triplet : public Eigen::Triplet<double>
	{
		Triplet(int i, int j, double v) : Eigen::Triplet<double>(i, j, v) {}
	};

	std::vector<Triplet> triplets;
	triplets.reserve(numVertices * 7);

	auto oppositeAngle = [&mesh](Mesh::HalfedgeHandle heh)
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
	};

	for (auto vi : mesh.vertices())
	{
		int i = vi.idx();
		double sumCotAlphaBeta = 0.0;

		for (auto voh : mesh.voh_range(vi))
		{
			auto vj = mesh.to_vertex_handle(voh);
			int j = vj.idx();

			double cotAlpha = 1.0 / tan(oppositeAngle(voh));
			double cotBeta = 0.0;
			if (mesh.opposite_halfedge_handle(voh).is_valid())
				cotBeta = 1.0 / tan(oppositeAngle(mesh.opposite_halfedge_handle(voh)));

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

Eigen::MatrixX2d tutteParameterization(const Mesh &mesh, Eigen::SparseMatrix<double> laplacian)
{
	int numVertices = mesh.numVertices();

	Eigen::MatrixX2d X(numVertices, 2);
	X.setZero();

	for (int bdrIdx : boundaryVertices(mesh))
	{
		X(bdrIdx, 0) = mesh.point(Mesh::VertexHandle(bdrIdx)).x();
		X(bdrIdx, 1) = mesh.point(Mesh::VertexHandle(bdrIdx)).y();
	}

	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(laplacian);

	X = solver.solve(X);

	return X;
}

Eigen::MatrixX2d localGlobalParameterization(Mesh &mesh, int numIter = 5)
{
	Eigen::SparseMatrix<double> laplacian = cotLaplacian(mesh);

	// Initial guess
	Eigen::MatrixX2d X = tutteParameterization(mesh, laplacian);

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

		triangleXs[f.idx()] << e01(0), e02(0),
			e01(1), e02(1);
	}

	for (int iter = 0; iter < numIter; iter++)
	{
		// Local step
		Eigen::MatrixX2d rhs;
		rhs.setZero();

		for (auto f : mesh.faces())
		{
			auto heh = mesh.halfedge_handle(f);
			auto v0 = mesh.from_vertex_handle(heh),
				 v1 = mesh.to_vertex_handle(heh),
				 v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(heh));

			Eigen::Matrix2d triangleU;
			triangleU << X(v1.idx(), 0) - X(v0.idx(), 0), X(v2.idx(), 0) - X(v0.idx(), 0),
				X(v1.idx(), 1) - X(v0.idx(), 1), X(v2.idx(), 1) - X(v0.idx(), 1);

			Eigen::Matrix2d J = triangleU * triangleXs[f.idx()].inverse();

			Eigen::Matrix2d R = bestFitARAP(J);

			Eigen::Matrix2d triangleUU = R * triangleXs[f.idx()];

			rhs.block<1, 2>(v0.idx(), 0) +=
		}

		// Global step
	}
}
