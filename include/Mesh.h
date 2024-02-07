#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::VectorT<double, 3> Point;
	// typedef OpenMesh::Vec3d Normal;
	// VertexAttributes(OpenMesh::Attributes::Normal);
	// FaceAttributes(OpenMesh::Attributes::Normal);
};

class Mesh : public OpenMesh::TriMesh_ArrayKernelT<MyTraits>
{
private:
public:
	Mesh(std::string &filename);
	Mesh(std::vector<double> vertexList, std::vector<unsigned int> faceList);
};