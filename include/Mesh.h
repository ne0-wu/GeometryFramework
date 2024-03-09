#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "GLMesh.h"

struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;

	VertexAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);
};

class Mesh : public OpenMesh::TriMesh_ArrayKernelT<MyTraits>
{
private:
	GLMesh glmesh;

public:
	Mesh(const std::string &filename);
	Mesh(const std::vector<double> &vertexList, const std::vector<unsigned int> &faceList);

	void save(const std::string &filename);

	std::vector<double> vertexList();
	std::vector<float> vertexListFloat();
	const std::vector<unsigned int> faceList();

	int numVertices()
	{
		return this->n_vertices();
	}

	int numFaces()
	{
		return this->n_faces();
	}

	// Geometry information
	// --------------------

	Eigen::Vector3d eigenPoint(VertexHandle v)
	{
		return Eigen::Vector3d(point(v)[0], point(v)[1], point(v)[2]);
	}

	Eigen::Vector3d eigenNormal(VertexHandle v)
	{
		return Eigen::Vector3d(normal(v)[0], normal(v)[1], normal(v)[2]);
	}

	Eigen::Vector3d eigenNormal(FaceHandle f)
	{
		return Eigen::Vector3d(normal(f)[0], normal(f)[1], normal(f)[2]).normalized();
	}

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
	bool simplifyQEM(int targetNumVertices);
	bool simplifyQEM2(int targetNumVertices);
	void collapseEdge(HalfedgeHandle edge, Eigen::Vector3d contractedPosition);
	void collapseEdge(HalfedgeHandle edge, Point contractedPosition);
	Eigen::Matrix4d quadricErrorMatrix(VertexHandle v);
	double quadricErrorEdge(const Eigen::Matrix4d &Q, HalfedgeHandle edge);
	Point optimalPlacement(const Eigen::Matrix4d &Q, HalfedgeHandle edge);
};