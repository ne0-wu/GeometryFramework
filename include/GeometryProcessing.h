#include <Eigen/Core>
#include <Eigen/Sparse>

#include "Mesh.h"
#include "PointCloud.h"

void poissonSurfaceReconstruction(PointCloud &pointCloud);

enum class LocalGlobalTarget
{
    ARAP,
    ASAP
};

Eigen::MatrixX2d localGlobalParameterization(Mesh &mesh, int numIter = 100, LocalGlobalTarget target = LocalGlobalTarget::ARAP);

class Parameterization
{
public:
    enum class LaplacianType
    {
        UNIFORM,
        COTANGENT
    };

private:
    Mesh const &mesh;    // input mesh
    Eigen::MatrixX2d uv; // output uv coordinates

    // Tutte parameterization
    void tutte(Eigen::SparseMatrix<double> const &laplacian);

    LaplacianType laplacianType = LaplacianType::UNIFORM;
    Eigen::SparseMatrix<double> laplacianUniform();
    Eigen::SparseMatrix<double> laplacianCotangent();

    Eigen::SparseMatrix<double> laplacian;

    // Pre-compute cotangent
    std::vector<double> cotangents;
    void computeCotangents();

public:
    Parameterization(Mesh const &mesh) : mesh(mesh) {}

    void setLaplacianType(LaplacianType type)
    {
        laplacianType = type;
    }

    void flatten();

    Eigen::MatrixX2d getUV() const { return uv; }
};