#pragma once

#include <iostream>
#include <string>
#include <vector>

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
public:
	Mesh(const std::string &filename);
	Mesh(const std::vector<double> &vertexList, const std::vector<unsigned int> &faceList);

	std::vector<double> vertexList();
	std::vector<float> vertexListFloat();
	const std::vector<unsigned int> faceList();

	int numVertices();
	int numFaces();

	Point boundingBoxMin();
	Point boundingBoxMax();

	Point center();		// center of bounding box
	Point barycenter(); // center of mass

	void move(Point delta);

	void moveCenterTo(Point destiny);
	void moveBarycenterTo(Point destiny);

	void moveCenterToOrigin();

	void scale(double scaler);

	void resize(Point min, Point max);
	void resize(double minX, double minY, double minZ, double maxX, double maxY, double maxZ);

	void fitInto(Mesh::Point min, Mesh::Point max);
	void fitInto(double minX, double minY, double minZ, double maxX, double maxY, double maxZ);
	void fitIntoUnitCube();
	void fitIntoUnitBall();

	bool simplifyQEM(int nVertices);

private:
	bool removeVertex();
	bool collapseEdge();
};