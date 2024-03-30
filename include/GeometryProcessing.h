#include "Mesh.h"
#include "PointCloud.h"

void poissonSurfaceReconstruction(PointCloud &pointCloud);

Eigen::MatrixX2d localGlobalParameterization(Mesh &mesh, int numIter = 5);