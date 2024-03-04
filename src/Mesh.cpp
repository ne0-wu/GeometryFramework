#include "Mesh.h"

Mesh::Mesh(const std::string &filename)
{
	// OpenMesh::IO::_OBJReader_();
	if (!OpenMesh::IO::read_mesh(*this, filename))
	{
		std::cerr << "Error loading mesh from file " << filename << std::endl;
		exit(EXIT_FAILURE);
	}
}

Mesh::Mesh(const std::vector<double> &vertexList, const std::vector<unsigned int> &faceList)
{
	for (int i = 0; i < vertexList.size(); i += 3)
		this->add_vertex(Mesh::Point(vertexList[i], vertexList[i + 1], vertexList[i + 2]));
	for (int i = 0; i < faceList.size(); i += 3)
		this->add_face(Mesh::VertexHandle(faceList[i]), Mesh::VertexHandle(faceList[i + 1]), Mesh::VertexHandle(faceList[i + 2]));
}

std::vector<float> Mesh::vertexList()
{
	std::vector<float> vertexList;
	for (Mesh::VertexIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it)
	{
		Mesh::Point p = this->point(*v_it);
		vertexList.push_back(p[0]);
		vertexList.push_back(p[1]);
		vertexList.push_back(p[2]);
	}
	return vertexList;
}

const std::vector<unsigned int> Mesh::faceList()
{
	std::vector<unsigned int> faceList;
	for (Mesh::FaceIter f_it = this->faces_begin(); f_it != this->faces_end(); ++f_it)
	{
		Mesh::FaceVertexIter fv_it;
		for (fv_it = this->fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			faceList.push_back(fv_it->idx());
		}
	}
	return faceList;
}