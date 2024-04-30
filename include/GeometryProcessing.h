#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <OpenMesh/Core/Utils/PropertyManager.hh>

#include "Mesh.h"
#include "PointCloud.h"

#define IMPLEMENT_QEM_SIMPLIFICATION
// #define IMPLEMENT_LOCAL_GLOBAL_PARA
// #define IMPLEMENT_POISSON_RECON
#define IMPLEMENT_CUBIC_STYLIZATION

class QEMSimplification
{
public:
	QEMSimplification(Mesh &mesh);

	void simplify(int targetNumVertices);
	Mesh &getMesh() { return mesh; };

private:
	struct QEM
	{
		Mesh::HalfedgeHandle heh;
		Mesh::Point optimalPlacement;
		double error = std::numeric_limits<double>::infinity();

		bool operator<(const QEM &rhs) const { return error < rhs.error; }
	};

	Mesh mesh; // The mesh to simplify

	// Store reusable variables as mesh properties
	OpenMesh::VProp<Eigen::Matrix4d> Qs; // Quadric error matrices
	OpenMesh::EProp<QEM> QEMs;			 // Quadric error on each edge

	Eigen::Matrix4d quadricErrorMatrix(Mesh::VertexHandle v);
	QEM optimalPlacement(Mesh::HalfedgeHandle edge, Eigen::Matrix4d Q);

	void collapse1Edge();
};

void poissonSurfaceReconstruction(PointCloud &pointCloud);

class Parameterization
{
public:
	Parameterization(Mesh const &mesh) : mesh(mesh) {}
	virtual void flatten() = 0;
	Eigen::MatrixX2d get_uv() { return uv; }

protected:
	Mesh const &mesh;	 // input mesh
	Eigen::MatrixX2d uv; // output uv coordinates
};

class Tutte : public Parameterization
{
public:
	enum class LaplacianType
	{
		UNIFORM,
		COTANGENT
	};

	Tutte(Mesh const &mesh, LaplacianType laplacianType = LaplacianType::UNIFORM)
		: Parameterization(mesh), laplacianType(laplacianType) {}

	void set_laplacian_type(LaplacianType type)
	{
		laplacianType = type;
	}

	void flatten();

protected:
	// Reusable variables
	std::vector<double> cotangents;
	std::vector<Mesh::VertexHandle> boundary_vertices;

	// Tutte parameterization
	void tutte(Eigen::SparseMatrix<double> const &laplacian);

	LaplacianType laplacianType;
	Eigen::SparseMatrix<double> laplacian_uniform();
	Eigen::SparseMatrix<double> laplacian_cotangent();

	Eigen::SparseMatrix<double> laplacian;

	// Opposite angle of a halfedge
	double opposite_angle(Mesh::HalfedgeHandle heh) const;

	// Pre-compute cotangent
	void compute_cot();

	// Find boundary vertices
	void find_boundary_vertices();

	void tutte();
};

class LocalGlobal : public Tutte
{
public:
	enum class LocalGlobalTarget
	{
		ARAP,
		ASAP
	};

	LocalGlobal(Mesh const &mesh, LocalGlobalTarget target = LocalGlobalTarget::ARAP, int numIter = 100)
		: Tutte(mesh), target(target), numIter(numIter) {}

	void set_num_iter(int numIter) { this->numIter = numIter; }
	int get_num_iter() { return numIter; }

	void flatten();

protected:
	int numIter;
	LocalGlobalTarget target;

	void localGlobal();
};

class CubicStylization
{
public:
	CubicStylization(Mesh &mesh, double lambda = 0.1, int numIter = 100);

	void stylize();
	Mesh get_stylized_mesh();

private:
	Mesh &mesh;
	double lambda;
	int numIter;

	OpenMesh::VProp<double> barycentric_area;
	OpenMesh::VProp<Eigen::Matrix3d> Rs;
	OpenMesh::HProp<double> cotangents;

	Eigen::MatrixX3d V;
	Eigen::SparseMatrix<double> laplacian;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> laplacian_solver;
	Eigen::MatrixX3d rhs;

	void local();
	void global();
	void local_global();
};

class GeodesicPath
{
public:
	GeodesicPath(Mesh const &input_mesh, Mesh::VertexHandle source, Mesh::VertexHandle target);

private:
	const Mesh &original_mesh;
	Mesh mesh;					 // The mesh to compute geodesic path on
	Mesh::VertexHandle src, tgt; // Source and target vertices

	std::vector<Mesh::VertexHandle> path; // The geodesic path
	// std::vector<Mesh::HalfedgeHandle> path; // The geodesic path

	// Dijkstra's algorithm
	void dijkstra();

	// Signpost data structure
	OpenMesh::EProp<double> edge_length;
	OpenMesh::HProp<double> direction;
	OpenMesh::VProp<double> angle_sum; // For updating direction

	double angle(Mesh::HalfedgeHandle h) const;

	void intrinsic_flip(Mesh::EdgeHandle e);

	void flip_out(int i, bool left_hand_side);

	void find_geodesic_path();
};