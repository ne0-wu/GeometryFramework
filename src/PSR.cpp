// This file implements the Poisson Surface Reconstruction (PSR) algorithm.

#include <random>
#include <queue>

#include <Eigen/Sparse>

#include "Mesh.h"

PointCloud Mesh::generatePointCloud(int numPoints, bool useNormals)
{
	update_normals();

	PointCloud pointCloud;
	pointCloud.reserve(numPoints); // Reserve space

	int eachFace = numPoints / n_faces() + 1; // Number of points per face
	int remaining = numPoints % n_faces();	  // Remaining points

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0, 1);

	// Generate random points
	for (auto f : faces())
	{
		auto v_it = fv_iter(f);
		Eigen::Vector3d p0 = eigenPoint(v_it);
		Eigen::Vector3d n0 = eigenNormal(v_it++);
		Eigen::Vector3d p1 = eigenPoint(v_it);
		Eigen::Vector3d n1 = eigenNormal(v_it++);
		Eigen::Vector3d p2 = eigenPoint(v_it);
		Eigen::Vector3d n2 = eigenNormal(v_it);

		for (int i = 0; i < eachFace; i++)
		{
			// Random barycentric coordinates
			double r0 = dis(gen);
			double r1 = dis(gen);
			if (r0 + r1 > 1)
			{
				r0 = 1 - r0;
				r1 = 1 - r1;
			}
			double r2 = 1 - r0 - r1;

			// Interpolate point and normal
			Eigen::Vector3d point = r0 * p0 + r1 * p1 + r2 * p2;
			Eigen::Vector3d normal = r0 * n0 + r1 * n1 + r2 * n2;

			pointCloud.points.push_back(point);
			if (useNormals)
				pointCloud.normals.push_back(normal.normalized());
		}

		remaining--;
		if (remaining == 0)
			eachFace--;
	}

	return pointCloud;
}

int childIndex(Eigen::Vector3d center, Eigen::Vector3d point)
{
	int index = 0;
	if (point.x() >= center.x())
		index |= 1;
	if (point.y() >= center.y())
		index |= 2;
	if (point.z() >= center.z())
		index |= 4;
	return index;
}

Eigen::Vector3d childCenter(Eigen::Vector3d center, double halfSize, int index)
{
	Eigen::Vector3d childCenter = center;
	if (index & 1)
		childCenter.x() += halfSize;
	else
		childCenter.x() -= halfSize;
	if (index & 2)
		childCenter.y() += halfSize;
	else
		childCenter.y() -= halfSize;
	if (index & 4)
		childCenter.z() += halfSize;
	else
		childCenter.z() -= halfSize;
	return childCenter;
}

struct Node
{
	// Node box
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	double size = 1.0;

	// Tree structure
	std::vector<Node> children;
	int depth = 0;

	// Contained point
	int index = -1;
	Eigen::Vector3d data;

	// Neighbor nodes for trilinear interpolation, size (if not empty) should be 8
	std::vector<Node *> neighbors;

	// Constructors
	Node() = default;

	Node(const Eigen::Vector3d &center, double size, int depth = 0)
		: center(center), size(size), depth(depth) {}

	// Check if node is a leaf
	bool isLeaf() const { return children.empty(); }

	// Check if node is empty
	bool isEmpty() const { return index == -1; }

	// Split node into 8 children
	void splitNode()
	{
		// Create children
		double halfSize = size / 2;
		children.reserve(8);
		for (int i = 0; i < 8; i++)
			children.emplace_back(childCenter(center, halfSize, i), halfSize, depth + 1);

		// Move point to children
		children[childIndex(center, data)].insertPoint(data, index);

		// Clear point
		index = -1;
	}

	// Insert point into the octree
	bool insertPoint(const Eigen::Vector3d &point, int index)
	{
		if (isLeaf()) // Node is a leaf
		{
			if (isEmpty()) // Node is empty
			{
				this->index = index;
				this->data = point;
				return true;
			}
			else // Node already contains a point
			{
				if (this->index == index)
					return true;

				// Split node and insert current point to children
				splitNode();

				// Insert new point
				return insertPoint(point, index);
			}
		}
		else // Node is not a leaf, insert point into children
			return children[childIndex(center, point)].insertPoint(point, index);
	}

	// Find leaf node containing the point
	Node *findLeaf(const Eigen::Vector3d &point)
	{
		if (isLeaf())
			return this;
		else
			return children[childIndex(center, point)].findLeaf(point);
	}

	// Base function for FEM
	double baseFunc(Eigen::Vector3d x)
	{
		auto spline = [](double t)
		{
			if (abs(t) <= 0.5)
				return 0.75 - t * t;
			else if (abs(t) <= 1.5)
				return 0.5 * pow(1.5 - abs(t), 2);
			else
				return 0.0;
		};

		x = (x - center) / size;
		return spline(x.x()) * spline(x.y()) * spline(x.z()) / pow(size, 3);
	}

	// Debug
	void print()
	{
		// Indent
		auto indent = [](int depth)
		{
			for (int i = 0; i < depth; i++)
				std::cout << "    ";
		};

		if (isLeaf()) // Leaf node
		{
			if (isEmpty()) // Empty leaf node
			{
				indent(depth);
				std::cout << "Node: " << center.transpose() << " Size: " << size << " Depth: " << depth << " EMPTY" << std::endl;
			}
			else // Non-empty leaf node
			{
				indent(depth);
				std::cout << "Node: " << center.transpose() << " Size: " << size << " Depth: " << depth << std::endl;
				indent(depth);
				bool isInside = (data.x() < center.x() + size) && (data.x() > center.x() - size) &&
								(data.y() < center.y() + size) && (data.y() > center.y() - size) &&
								(data.z() < center.z() + size) && (data.z() > center.z() - size);
				std::cout << "└-Point: " << data.transpose() << " Index: " << index << " is inside: " << isInside << std::endl;
			}
		}
		else // Non-leaf node
		{
			indent(depth);
			std::cout << "Node: " << center.transpose() << " Size: " << size << " Depth: " << depth << std::endl;
		}

		// Print children
		for (auto &child : children)
			child.print();
	}
};

class Octree
{
private:
	Node root;

	// Trade time with space
	int maxDepth = 0;
	std::vector<Node *> leafNodes;
	std::vector<Node *> nonEmptyLeafNodes;

	// Max depth of the octree
	void updateMaxDepth()
	{
		maxDepth = 0;
		int depth = 0;
		traverse([&depth](Node &node)
				 {
				if (node.depth > depth)
					depth = node.depth; });
		maxDepth = depth;
	}

	void updateLeafNodes()
	{
		leafNodes.clear();
		traverse([&](Node &node)
				 {
				if (node.isLeaf())
					leafNodes.push_back(&node); });
	}

	void updateNonEmptyLeafNodes()
	{
		nonEmptyLeafNodes.clear();
		traverse([&](Node &node)
				 {
				if (node.isLeaf() && !node.isEmpty())
					nonEmptyLeafNodes.push_back(&node); });
	}

public:
	Octree(const PointCloud &pointCloud, Eigen::Vector3d center, double size)
		: root(center, size)
	{
		for (int i = 0; i < pointCloud.points.size(); i++)
			insertPoint(pointCloud.points[i], i);
	}

	// Traverse the octree with a callback function
	void traverse(std::function<void(Node &)> callback)
	{
		std::queue<Node *> queue;
		queue.push(&root);
		while (!queue.empty())
		{
			Node *node = queue.front();
			queue.pop();

			callback(*node);

			if (!node->isLeaf())
				for (auto &child : node->children)
					queue.push(&child);
		}
	}

	int getMaxDepth() { return maxDepth; }
	auto getLeafNodeList() { return leafNodes; }
	auto getNonEmptyLeafNodeList() { return nonEmptyLeafNodes; }

	// Refine the octree
	void refine()
	{
		updateMaxDepth();

		// Split nodes until all non-empty leaf nodes are at the maximum depth
		std::queue<Node *> queue;
		queue.push(&root);
		while (!queue.empty())
		{
			Node *node = queue.front();
			queue.pop();

			if (node->isLeaf() && !node->isEmpty() && node->depth < maxDepth)
				node->splitNode();

			if (!node->isLeaf())
				for (auto &child : node->children)
					queue.push(&child);
		}

		// TODO: Split nodes until all neighbors of non-empty leaf nodes are at the maximum depth
		updateNonEmptyLeafNodes();
		for (auto &node : nonEmptyLeafNodes)
		{
			;
		}

		updateLeafNodes();
	}

	// TODO: Laplacian matrix
	Eigen::SparseMatrix<double> laplacianMatrix()
	{
		int n = leafNodes.size();
		Eigen::SparseMatrix<double> L(n, n);
		L.reserve(50 * n);

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j)
				{
				}
				else
				{
				}
			}
		}
	}

	// TODO: Indicator function

	// Debug
	void print() { root.print(); }

private:
	void insertPoint(Eigen::Vector3d point, int index)
	{
		root.insertPoint(point, index);
	}
};

void testOctree(PointCloud pointCloud)
{
	Octree tree(pointCloud, Eigen::Vector3d::Zero(), 1.0);
	tree.print();

	tree.refine();

	tree.print();
}