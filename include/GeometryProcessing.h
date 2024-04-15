#include <Eigen/Core>
#include <Eigen/Sparse>

#include "Mesh.h"
#include "PointCloud.h"

void poissonSurfaceReconstruction(PointCloud &pointCloud);

class TutteParameterization
{
public:
	enum class LaplacianType
	{
		UNIFORM,
		COTANGENT
	};

	TutteParameterization(Mesh const &mesh) : mesh(mesh) {}

	void setLaplacianType(LaplacianType type)
	{
		laplacianType = type;
	}

	void flatten();

	Eigen::MatrixX2d getUV() const { return uv; }

protected:
	Mesh const &mesh;	 // input mesh
	Eigen::MatrixX2d uv; // output uv coordinates

	// Reusable variables
	std::vector<double> cotangents;
	std::vector<Mesh::VertexHandle> boundaryVertices;

	// Tutte parameterization
	void tutte(Eigen::SparseMatrix<double> const &laplacian);

	LaplacianType laplacianType = LaplacianType::COTANGENT;
	Eigen::SparseMatrix<double> laplacianUniform();
	Eigen::SparseMatrix<double> laplacianCotangent();

	Eigen::SparseMatrix<double> laplacian;

	// Opposite angle of a halfedge
	double oppositeAngle(Mesh::HalfedgeHandle heh) const;

	// Pre-compute cotangent
	void computeCotangents();

	// Find boundary vertices
	void findBoundaryVertices();

	void tutte();
};

class LocalGlobalParameterization : public TutteParameterization
{
protected:
	void localGlobal();

public:
	enum class LocalGlobalTarget
	{
		ARAP,
		ASAP
	} target;
	int numIter;

	LocalGlobalParameterization(Mesh const &mesh, LocalGlobalTarget target = LocalGlobalTarget::ARAP, int numIter = 10)
		: TutteParameterization(mesh), target(target), numIter(numIter) {}

	void flatten();
};