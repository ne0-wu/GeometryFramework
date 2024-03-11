// This file implements the Poisson Surface Reconstruction (PSR) algorithm.

#include <random>

#include "Mesh.h"

PointCloud Mesh::pointCloud(int numPoints, bool useNormals)
{
    PointCloud pointCloud;
    pointCloud.points.reserve(numPoints); // Reserve space for the points

    int eachFace = numPoints / n_faces() + 1; // Number of points per face
    int remaining = numPoints % n_faces();    // Remaining points

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0, 1);

    // Generate random points
    for (auto f : faces())
    {
        auto he_h = halfedge_handle(f);
        Eigen::Vector3d p0 = eigenPoint(to_vertex_handle(he_h));
        Eigen::Vector3d p1 = eigenPoint(to_vertex_handle(next_halfedge_handle(he_h)));
        Eigen::Vector3d p2 = eigenPoint(to_vertex_handle(next_halfedge_handle(next_halfedge_handle(he_h))));
        Eigen::Vector3d n0 = eigenNormal(to_vertex_handle(he_h));
        Eigen::Vector3d n1 = eigenNormal(to_vertex_handle(next_halfedge_handle(he_h)));
        Eigen::Vector3d n2 = eigenNormal(to_vertex_handle(next_halfedge_handle(next_halfedge_handle(he_h))));

        for (int i = 0; i < eachFace; i++)
        {
            // Random barycentric coordinates
            double r0 = dis(gen);
            double r1 = dis(gen);
            if (r0 + r1 > 1)
            {
                r0 = 1 - r0;
                r1 = 1 - r1;
            }
            double r2 = 1 - r0 - r1;

            // Interpolate point and normal
            Eigen::Vector3d point = r0 * p0 + r1 * p1 + r2 * p2;
            Eigen::Vector3d normal = r0 * n0 + r1 * n1 + r2 * n2;

            pointCloud.points.push_back(point);
            if (useNormals)
                pointCloud.normals.push_back(normal.normalized());
        }

        remaining--;
        if (remaining == 0)
            eachFace--;
    }

    return pointCloud;
}