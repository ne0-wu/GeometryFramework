#include "Mesh.h"

Mesh::Mesh(std::string &filename)
{
	if (!OpenMesh::IO::read_mesh(*this, filename))
	{
		std::cerr << "Error loading mesh from file " << filename << std::endl;
		exit(EXIT_FAILURE);
	}
}

Mesh::Mesh(const std::vector<std::vector<int>> &faceList, const std::vector<std::vector<float>> &vertexList)
{
	for (const auto &face : faceList)
	{
		std::vector<VertexHandle> faceVertices;
		for (const auto &vertexIndex : face)
		{
			VertexHandle vh = this->add_vertex(Point(vertexList[vertexIndex][0], vertexList[vertexIndex][1], vertexList[vertexIndex][2]));
			faceVertices.push_back(vh);
		}
		this->add_face(faceVertices);
	}
}