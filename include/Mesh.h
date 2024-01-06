#pragma once

#include <iostream>
#include <string>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;
	// VertexAttributes(OpenMesh::Attributes::Normal);
	// FaceAttributes(OpenMesh::Attributes::Normal);
};

class Mesh : public OpenMesh::TriMesh_ArrayKernelT<MyTraits>
{
private:
	std::vector<std::vector<int>> faceList;
	std::vector<std::vector<float>> vertexList;

public:
	/**
	 * @brief Constructs a Mesh object from a file.
	 *
	 * @param filename The name of the file to load the mesh from.
	 */
	Mesh(std::string &filename);

	/**
	 * @brief Constructs a Mesh object from vertex and face data.
	 *
	 * @param vertices A vector of vectors representing the vertex positions.
	 * @param faces A vector of vectors representing the face indices.
	 */
	Mesh(const std::vector<std::vector<int>> &vertices, const std::vector<std::vector<float>> &faces);
};