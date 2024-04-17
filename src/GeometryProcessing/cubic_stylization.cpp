// This file contains the implementation of the cubic stylization algorithm.

#include "GeometryProcessing.h"

CubicStylization::CubicStylization(Mesh &input_mesh, double lambda, int numIter)
	: mesh(input_mesh), lambda(lambda), numIter(numIter),
	  barycentric_area(mesh), Rs(mesh), cotangents(mesh),
	  laplacian(mesh.n_vertices(), mesh.n_vertices()), rhs(mesh.n_vertices(), 3), V(mesh.n_vertices(), 3)
{
	// Request and comput normals
	mesh.request_face_normals();
	mesh.request_vertex_normals();
	mesh.update_normals();

	// Compute barycentric area
	OpenMesh::FProp<double> face_area(mesh);
	for (auto f : mesh.faces())
		face_area(f) = mesh.calc_sector_area(mesh.halfedge_handle(f)) / 3;
	for (auto v : mesh.vertices())
	{
		double area = 0;
		for (auto f : mesh.vf_range(v))
			area += face_area(f);
		barycentric_area[v] = area;
		V.row(v.idx()) = mesh.point(v); // Initialize V as identical to input mesh
	}

	// Compute cotangents and laplacian
	std::vector<Eigen::Triplet<double>> triplets;
	triplets.reserve(mesh.n_halfedges() * 4);
	for (auto he : mesh.halfedges())
	{
		if (he.is_boundary())
			continue;

		auto fr = mesh.point(he.from()),
			 to = mesh.point(he.to()),
			 op = mesh.point(he.next().to());
		double cot_op = 1.0 / tan(acos((fr - op).normalized().dot((to - op).normalized())));

		cotangents[he] = cot_op;

		triplets.push_back(Eigen::Triplet<double>(he.to().idx(), he.to().idx(), cot_op));
		triplets.push_back(Eigen::Triplet<double>(he.from().idx(), he.from().idx(), cot_op));
		triplets.push_back(Eigen::Triplet<double>(he.to().idx(), he.from().idx(), -cot_op));
		triplets.push_back(Eigen::Triplet<double>(he.from().idx(), he.to().idx(), -cot_op));
	}
	laplacian.setFromTriplets(triplets.begin(), triplets.end());
	laplacian_solver.compute(laplacian);
}

void CubicStylization::local()
{
	const double eps_abs = 1e-5,
				 eps_rel = 1e-3,
				 mu = 10.0,
				 tau_incr = 2.0,
				 tau_decr = 2.0;

	for (auto vi : mesh.vertices())
	{
		// Reusable variables
		Eigen::Vector3d normal = mesh.normal(vi);
		int deg = vi.valence();
		double ai = barycentric_area(vi);

		// The matrix A in the ADMM constraint Ax + Bz = c,
		// where x is R flattened, and B = I, c = 0.
		Eigen::MatrixXd A(3, 9);
		A.setZero();
		A.block<1, 3>(0, 0) = normal.transpose();
		A.block<1, 3>(1, 3) = normal.transpose();
		A.block<1, 3>(2, 6) = normal.transpose();

		Eigen::Matrix3d R;

		double rho = 1e-4;
		Eigen::Vector3d z, z_prev, u;
		u.setZero();
		z = mesh.normal(vi);

		// Pre calculating of Mi
		// Mi = D * Mi_diag * D_tilde^T
		// In every iteration, the last column of D_tilde and the last element of Mi_diag are updated
		Eigen::Matrix3Xd D(3, deg + 1),
			D_tilde(3, deg + 1);
		Eigen::MatrixXd Mi_diag(deg + 1, deg + 1);
		Mi_diag.setZero();
		int j = 0;
		for (auto heh : mesh.voh_range(vi))
		{
			auto vj = mesh.to_vertex_handle(heh);
			D.col(j) = mesh.point(vj) - mesh.point(vi);
			D_tilde.col(j) = V.row(vj.idx()) - V.row(vi.idx());
			Mi_diag(j, j) = (cotangents[heh] + cotangents[heh.opp()]);
			j++;
		}
		D.col(deg) = mesh.normal(vi);

		// Find the optimal R with ADMM
		for (int numIter = 0; numIter < 100; numIter++)
		{
			// Update R
			Mi_diag(deg, deg) = rho;
			D_tilde.col(deg) = z - u;
			Eigen::Matrix3d Mi = D * Mi_diag * D_tilde.transpose();
			Eigen::JacobiSVD svd(Mi, Eigen::ComputeFullU | Eigen::ComputeFullV);
			Eigen::Matrix3d U = svd.matrixU(), V = svd.matrixV();
			// Flip a row of U so that R is positively oriented
			if ((V * U.transpose()).determinant() < 0)
				U.col(2) *= -1;
			R = V * U.transpose();

			// Update z
			z_prev = z;
			Eigen::Vector3d x = R * normal + u;
			double k = lambda * ai / rho;
			for (int j : {0, 1, 2})
			{
				double coe = 1 - k / abs(x(j));
				coe = coe > 0 ? coe : 0;
				z(j) = coe * x(j);
			}

			// Update u
			u += R * normal - z;

			// Update rho
			double r = (R * normal - z).norm(),			   // Primal residual
				s = (A.transpose() * (z - z_prev)).norm(); // Dual residual

			if (r > mu * s)
				rho *= tau_incr, u /= tau_incr;
			else if (s > mu * r)
				rho /= tau_decr, u *= tau_decr;

			// Stop condition
			if (r < eps_abs && s < eps_rel)
				break;
		}

		Rs[vi] = R;
	}
}

void CubicStylization::global()
{
	rhs.setZero();
	for (auto vi : mesh.vertices())
		for (auto he : mesh.voh_range(vi))
		{
			auto vj = he.to();
			double wij = cotangents[he] + cotangents[he.opp()];
			Eigen::Matrix3d R = (Rs[vi] + Rs[vj]) / 2;
			rhs.row(vi.idx()) += wij * R * (mesh.point(vi) - mesh.point(vj));
		}

	V = laplacian_solver.solve(rhs);
}

void CubicStylization::local_global()
{
	for (int i = 1; i < numIter; i++)
	{
		local();
		global();
	}
}

void CubicStylization::stylize()
{
	local_global();
}

Mesh CubicStylization::get_stylized_mesh()
{
	Mesh stylized_mesh = mesh;
	for (auto vi : mesh.vertices())
		stylized_mesh.set_point(vi, V.row(vi.idx()).transpose());
	return stylized_mesh;
}