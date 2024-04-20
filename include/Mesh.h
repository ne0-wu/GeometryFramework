#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>

#include <Eigen/Dense>

struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef Eigen::Vector3d Point;
	typedef Eigen::Vector3d Normal;
	typedef Eigen::Vector2d TexCoord2D;

	VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	EdgeAttributes(OpenMesh::Attributes::Status);
};

class Mesh : public OpenMesh::TriMesh_ArrayKernelT<MyTraits>
{
public:
	Mesh(const std::string &filename);
	Mesh(const std::vector<double> &vertexList, const std::vector<unsigned int> &faceList);

	void save(const std::string &filename);

	std::vector<double> vertexList();
	std::vector<float> vertexListFloat();
	const std::vector<unsigned int> faceList();

	int numVertices() const { return this->n_vertices(); }

	int numFaces() const { return this->n_faces(); }

	int numEdges() const { return this->n_edges(); }

	int numHalfEdges() const { return this->n_halfedges(); }

	// Smart access to mesh elements
	// -----------------------------
	template <typename... Args>
	OpenMesh::SmartVertexHandle smart_vertex_handle(Args &&...args)
	{
		return OpenMesh::make_smart(vertex_handle(std::forward<Args>(args)...), *this);
	}

	template <typename... Args>
	OpenMesh::SmartHalfedgeHandle smart_halfedge_handle(Args &&...args)
	{
		return OpenMesh::make_smart(halfedge_handle(std::forward<Args>(args)...), *this);
	}

	template <typename... Args>
	OpenMesh::SmartEdgeHandle smart_edge_handle(Args &&...args)
	{
		return OpenMesh::make_smart(edge_handle(std::forward<Args>(args)...), *this);
	}

	template <typename... Args>
	OpenMesh::SmartFaceHandle smart_face_handle(Args &&...args)
	{
		return OpenMesh::make_smart(face_handle(std::forward<Args>(args)...), *this);
	}

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

	// Boundary
	// --------
	// TODO: Need to consider the case of multiple boundary components
	std::vector<Mesh::VertexHandle> boundary_vertices();

	// Basic geometry transformations
	// ------------------------------
	void scale(double scaler);

	void resize(Point min, Point max);
	void resize(double minX, double minY, double minZ, double maxX, double maxY, double maxZ);

	void fitInto(Mesh::Point min, Mesh::Point max);
	void fitInto(double minX, double minY, double minZ, double maxX, double maxY, double maxZ);
	void fitIntoUnitCube();
	void fitIntoUnitBall();
};