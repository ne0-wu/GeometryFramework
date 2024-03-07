#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include "GLMesh.h"

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
	GLMesh glmesh;

public:
	Mesh(const std::string &filename);
	Mesh(const std::vector<double> &vertexList, const std::vector<unsigned int> &faceList);

	std::vector<double> vertexList();
	std::vector<float> vertexListFloat();
	const std::vector<unsigned int> faceList();

	int numVertices();
	int numFaces();

	// Geometry information
	// --------------------

	Point boundingBoxMin();
	Point boundingBoxMax();

	Point center();		// center of bounding box
	Point barycenter(); // center of mass

	void move(Point delta);

	void moveCenterTo(Point destiny);
	void moveBarycenterTo(Point destiny);

	void moveCenterToOrigin();

	// Basic geometry transformations
	// ------------------------------

	void scale(double scaler);

	void resize(Point min, Point max);
	void resize(double minX, double minY, double minZ, double maxX, double maxY, double maxZ);

	void fitInto(Mesh::Point min, Mesh::Point max);
	void fitInto(double minX, double minY, double minZ, double maxX, double maxY, double maxZ);
	void fitIntoUnitCube();
	void fitIntoUnitBall();

	// Mesh simplification
	// -------------------
	bool simplifyQEM(int nVertices);
	bool removeVertex();
	bool collapseEdge();
};