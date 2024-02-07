#include "Mesh.h"

Mesh::Mesh(std::string &filename)
{
	if (!OpenMesh::IO::read_mesh(*this, filename))
	{
		std::cerr << "Error loading mesh from file " << filename << std::endl;
		exit(EXIT_FAILURE);
	}
}

Mesh::Mesh(std::vector<double> vertexList, std::vector<unsigned int> faceList)
{
	for (int i = 0; i < vertexList.size(); i += 3)
		this->add_vertex(Mesh::Point(vertexList[i], vertexList[i + 1], vertexList[i + 2]));
	for (int i = 0; i < faceList.size(); i += 3)
		this->add_face(Mesh::VertexHandle(faceList[i]), Mesh::VertexHandle(faceList[i + 1]), Mesh::VertexHandle(faceList[i + 2]));
}