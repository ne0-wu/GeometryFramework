#pragma once

#include <iostream>
#include <string>
#include <vector>

#define OM_STATIC_BUILD
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/reader/OBJReader.hh>
// #include <OpenMesh/Core/IO/writer/OBJWriter.hh>

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
	Mesh(const std::string &filename);
	Mesh(const std::vector<double> &vertexList, const std::vector<unsigned int> &faceList);

	std::vector<float> vertexList();
	const std::vector<unsigned int> faceList();
};