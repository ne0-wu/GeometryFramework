#include "Mesh.h"
#include "PointCloud.h"

void poissonSurfaceReconstruction(PointCloud &pointCloud);

enum class LocalGlobalTarget
{
    ARAP,
    ASAP
};

Eigen::MatrixX2d localGlobalParameterization(Mesh &mesh, int numIter = 100, LocalGlobalTarget target = LocalGlobalTarget::ARAP);