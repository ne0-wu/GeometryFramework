#ifndef OM_STATIC_BUILD
#define OM_STATIC_BUILD
#endif

#include <OpenMesh/Core/IO/MeshIO.hh>

#include "Mesh.h"

Mesh::Mesh(const std::string &filename)
{
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

void Mesh::save(const std::string &filename)
{
	garbage_collection();
	if (!OpenMesh::IO::write_mesh(*this, filename))
	{
		std::cerr << "Error saving mesh to file " << filename << std::endl;
		exit(EXIT_FAILURE);
	}
}

std::vector<float> Mesh::vertexListFloat()
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

std::vector<double> Mesh::vertexList()
{
	std::vector<double> vertexList(numVertices() * 3);
	for (Mesh::VertexIter v_it = this->vertices_begin(); v_it != this->vertices_end(); ++v_it)
	{
		Mesh::Point p = this->point(*v_it);
		vertexList[3 * v_it->idx()] = p[0];
		vertexList[3 * v_it->idx() + 1] = p[1];
		vertexList[3 * v_it->idx() + 2] = p[2];
	}
	return vertexList;
}

const std::vector<unsigned int> Mesh::faceList()
{
	std::vector<unsigned int> faceList(numFaces() * 3);
	for (Mesh::FaceIter f_it = this->faces_begin(); f_it != this->faces_end(); ++f_it)
	{
		Mesh::ConstFaceVertexIter fv_it = this->cfv_iter(*f_it);
		faceList[3 * f_it->idx()] = fv_it->idx();
		++fv_it;
		faceList[3 * f_it->idx() + 1] = fv_it->idx();
		++fv_it;
		faceList[3 * f_it->idx() + 2] = fv_it->idx();
	}
	return faceList;
}

Mesh::Point Mesh::boundingBoxMin()
{
	double minX, minY, minZ;
	minX = minY = minZ = INFINITY;
	std::for_each(vertices_begin(), vertices_end(), [&](Mesh::VertexHandle v_it)
				  { Mesh::Point p = point(v_it); 
					minX = std::min(minX, p[0]);
					minY = std::min(minY, p[1]);
					minZ = std::min(minZ, p[2]); });
	return Mesh::Point(minX, minY, minZ);
}

Mesh::Point Mesh::boundingBoxMax()
{
	double maxX, maxY, maxZ;
	maxX = maxY = maxZ = -INFINITY;
	std::for_each(vertices_begin(), vertices_end(), [&](Mesh::VertexHandle v_it)
				  { Mesh::Point p = point(v_it); 
					maxX = std::max(maxX, p[0]);
					maxY = std::max(maxY, p[1]);
					maxZ = std::max(maxZ, p[2]); });
	return Mesh::Point(maxX, maxY, maxZ);
}

Mesh::Point Mesh::center()
{
	return (boundingBoxMin() + boundingBoxMax()) / 2;
}

Mesh::Point Mesh::barycenter()
{
	Mesh::Point center(0, 0, 0);
	std::for_each(vertices_begin(), vertices_end(), [&](Mesh::VertexHandle v_it)
				  { center += point(v_it); });
	center /= numVertices();
	return center;
}

void Mesh::move(Mesh::Point delta)
{
	std::for_each(vertices_begin(), vertices_end(), [&](Mesh::VertexHandle v_it)
				  { set_point(v_it, point(v_it) + delta); });
}

void Mesh::moveCenterTo(Mesh::Point destiny)
{
	move(destiny - center());
}

void Mesh::moveBarycenterTo(Mesh::Point destiny)
{
	move(destiny - barycenter());
}

void Mesh::moveCenterToOrigin()
{
	moveCenterTo(Mesh::Point(0, 0, 0));
}

void Mesh::scale(double scaler)
{
	Mesh::Point center = this->center();
	std::for_each(vertices_begin(), vertices_end(), [&](Mesh::VertexHandle v_it)
				  { set_point(v_it, (point(v_it) - center) * scaler + center); });
}

void Mesh::resize(Mesh::Point min, Mesh::Point max)
{
	Mesh::Point destCenter = (min + max) / 2, currCenter = this->center();
	auto destScale = max - destCenter;
	auto currScale = this->center() - currCenter;
	std::for_each(vertices_begin(), vertices_end(), [&](Mesh::VertexHandle v_it)
				  { set_point(v_it, (point(v_it) - currCenter) * destScale / currScale + destCenter); });
}

void Mesh::resize(double minX, double minY, double minZ, double maxX, double maxY, double maxZ)
{
	resize(Mesh::Point(minX, minY, minZ), Mesh::Point(maxX, maxY, maxZ));
}

void Mesh::fitInto(Mesh::Point min, Mesh::Point max)
{
	Mesh::Point destCenter = (min + max) / 2, currCenter = center(),
				destScale = max - destCenter, currScale = boundingBoxMax() - center();
	double scaler = (destScale / currScale).min();
	std::for_each(vertices_begin(), vertices_end(), [&](Mesh::VertexHandle v_it)
				  { set_point(v_it, (point(v_it) - currCenter) * scaler + destCenter); });
}

void Mesh::fitInto(double minX, double minY, double minZ, double maxX, double maxY, double maxZ)
{
	fitInto(Mesh::Point(minX, minY, minZ), Mesh::Point(maxX, maxY, maxZ));
}

void Mesh::fitIntoUnitCube()
{
	fitInto(Mesh::Point(0, 0, 0), Mesh::Point(1, 1, 1));
}

void Mesh::fitIntoUnitBall()
{
	moveCenterToOrigin();
	double maxRadius = 0;
	std::for_each(vertices_begin(), vertices_end(), [&](Mesh::VertexHandle v_it)
				  { maxRadius = std::max(maxRadius, point(v_it).norm()); });
	scale(1 / maxRadius);
}