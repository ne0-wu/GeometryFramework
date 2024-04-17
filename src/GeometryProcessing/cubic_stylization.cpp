// This file contains the implementation of the cubic stylization algorithm.

#include "GeometryProcessing.h"

CubicStylization::CubicStylization(Mesh &input_mesh, double lambda = 0.1, int numIter)
	: mesh(input_mesh), lambda(lambda), numIter(numIter),
	  barycentric_area(mesh), Rs(mesh), cotangents(mesh)
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
		V.row(v.idx()) = mesh.point(v); // Initialize V
	}

	// Compute cotangents
	for (auto he : mesh.halfedges())
	{
		if (he.is_boundary())
			continue;

		auto p0 = mesh.point(he.from()),
			 p1 = mesh.point(he.to()),
			 p2 = mesh.point(he.next().to());
		cotangents[he] = 1.0 / tan(acos((p0 - p2).normalized().dot((p1 - p2).normalized())));
	}
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

		Eigen::Matrix3d R;
		// R.setIdentity();
		// {
		// 	Eigen::Matrix3Xd D(3, vi.valence()),
		// 		D_tilde(3, vi.valence());
		// 	int j = 0;
		// 	for (auto heh : mesh.voh_range(vi))
		// 	{
		// 		auto vj = mesh.to_vertex_handle(heh);
		// 		D.col(j) = mesh.point(vj) - mesh.point(vi);
		// 		D_tilde.col(j) = V.row(vj.idx()) - V.row(vi.idx());
		// 		j++;
		// 	}
		// 	// Solve for J * D * D^T == D_tilde * D^T
		// 	Eigen::Matrix3d J = ((D.transpose()).colPivHouseholderQr().solve(D_tilde.transpose())).transpose();
		// 	Eigen::JacobiSVD<Eigen::Matrix3d> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
		// 	R = svd.matrixU() * svd.matrixV().transpose();
		// }

		double rho = 1e-4;
		Eigen::Vector3d z, u;
		u.setZero();
		z = mesh.normal(vi);

		// Pre calculating of Mi
		// Mi = D * Mi_diag * D_tilde^T
		// In every iteration, the last column of D_tilde and the last element of Mi_diag are updated
		Eigen::Matrix3Xd D(3, deg + 1),
			D_tilde(3, deg + 1);
		Eigen::MatrixXd Mi_diag(deg + 1);
		int j = 0;
		for (auto heh : mesh.voh_range(vi))
		{
			auto vj = mesh.to_vertex_handle(heh);
			D.col(j) = mesh.point(vj) - mesh.point(vi);
			D_tilde.col(j) = V.row(vj.idx()) - V.row(vi.idx());
			Mi_diag(j, j) = (cotangents[heh] + cotangents[heh.opp()]) / 2;
		}
		D.col(deg) = mesh.normal(vi);

		for (int numIter = 0; numIter < 100; numIter++)
		{
			// Update R
			Mi_diag(deg, deg) = rho;
			D_tilde.col(deg) = z - u;
			Eigen::Matrix3d Mi = D * Mi_diag * D_tilde.transpose();
			Eigen::JacobiSVD svd(Mi, Eigen::ComputeFullU | Eigen::ComputeFullV);
			R = svd.matrixV() * svd.matrixU().transpose();

			// Update z
			Eigen::Vector3d x = R * normal + u;
			double k = lambda * ai / rho;
			for (int j : {0, 1, 2})
			{
				double coe = 1 - k / abs(x(j));
				coe = coe > 0 ? coe : 0;
				z(j) = coe * x(j);
			}

			// Update u
			u = u + R * normal - z;

			// Update rho
			Eigen::Vector3d r, s; // primal and dual residuals
			r = R * normal - z;
			// TODO: compute the dual residual
		}
	}
}