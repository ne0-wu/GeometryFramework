// This file contains the implementation of the cubic stylization algorithm.

#include "GeometryProcessing.h"

CubicStylization::CubicStylization(Mesh &input_mesh, double lambda = 0.1, int numIter)
    : mesh(input_mesh), barycentric_area(mesh), Rs(mesh),
      lambda(lambda), numIter(numIter)
{
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
}

void CubicStylization::local()
{
    for (auto v : mesh.vertices())
    {
        // Solve R with ADMM
    }
}