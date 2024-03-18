#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Core>
#include <Eigen/Dense>

#include "GLMesh.h"

struct PointCloud
{
	std::vector<Eigen::Vector3d> points;
	std::vector<Eigen::Vector3d> normals;

	void reserve(int numPoints)
	{
		points.reserve(numPoints);
		normals.reserve(numPoints);
	}
};

struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;

	VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
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

	int numVertices() { return this->n_vertices(); }

	int numFaces() { return this->n_faces(); }

	// Render
	// ------
	void initializeGLMesh() { glmesh = GLMesh(vertexListFloat(), faceList()); }

	void initializeShader(const std::string &vertexShaderSource, const std::string &fragmentShaderSource, bool fromFile = true)
	{
		glmesh.initializeShader(vertexShaderSource, fragmentShaderSource, fromFile);
	}

	int getShaderID() { return glmesh.getShaderID(); }

	void bindShader() { glmesh.bindShader(); }

	void draw()
	{
		glmesh.setMesh(vertexListFloat(), faceList());

		int vertexColorLocation = glGetUniformLocation(getShaderID(), "color");
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glUniform4f(vertexColorLocation, 1.0f, 1.0f, 1.0f, 1.0f);
		glmesh.draw();

		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(2.0f);
		glUniform4f(vertexColorLocation, 0.0f, 0.0f, 0.0f, 1.0f);
		glmesh.draw();
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
	Eigen::Matrix4d quadricErrorMatrix(VertexHandle v);
	std::pair<double, Point> optimalPlacement(HalfedgeHandle edge, const Eigen::Matrix4d &Q);
	void collapseEdge(HalfedgeHandle edge, Point contractedPosition);

	// Generate point cloud
	// --------------------
	PointCloud generatePointCloud(int numPoints, bool useNormals = false);
};

// debug
void testOctree(PointCloud pointCloud);