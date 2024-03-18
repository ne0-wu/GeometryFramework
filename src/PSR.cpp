// This file implements the Poisson Surface Reconstruction (PSR) algorithm.

#include <random>
#include <queue>

#include "Mesh.h"

PointCloud Mesh::pointCloud(int numPoints, bool useNormals)
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

struct Node
{
	// Node box
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	double size = 1.0;

	// Tree structure
	std::vector<Node> children;

	// Contained point
	int index = -1;
	Eigen::Vector3d data;

	// Constructors
	Node() = default;

	Node(const Eigen::Vector3d &center, double size)
		: center(center), size(size) {}

	// Check if node is a leaf
	bool isLeaf() const { return children.empty(); }

	// Check if node is empty
	bool isEmpty() const { return index == -1; }

	// Insert points
	void splitNode()
	{
		double halfSize = size / 2;
		children.reserve(8);
		for (int i = 0; i < 8; i++)
		{
			Eigen::Vector3d childCenter = center;
			if (i & 1)
				childCenter.x() += halfSize;
			else
				childCenter.x() -= halfSize;
			if (i & 2)
				childCenter.y() += halfSize;
			else
				childCenter.y() -= halfSize;
			if (i & 4)
				childCenter.z() += halfSize;
			else
				childCenter.z() -= halfSize;
			children.emplace_back(childCenter, halfSize);
		}
	}

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

				// Split node and insert both points
				splitNode();
				bool result = true;
				result &= insertPoint(point, index);
				result &= insertPoint(this->data, this->index);

				// Clear point
				this->index = -1;

				return result;
			}
		}
		else // Node is not a leaf, insert point into children
		{
			int childIndex = 0;
			if (point.x() >= center.x())
				childIndex |= 1;
			if (point.y() >= center.y())
				childIndex |= 2;
			if (point.z() >= center.z())
				childIndex |= 4;

			return children[childIndex].insertPoint(point, index);
		}
	}

	// Function
	// TODO:

	// Debug
	void print(int depth = 0)
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
				std::cout << "Node: " << center.transpose() << " " << size << " EMPTY" << std::endl;
			}
			else // Non-empty leaf node
			{
				indent(depth);
				std::cout << "Node: " << center.transpose() << " " << size << " " << index << std::endl;
				indent(depth);
				bool isInside = (data.x() < center.x() + size) && (data.x() > center.x() - size) &&
								(data.y() < center.y() + size) && (data.y() > center.y() - size) &&
								(data.z() < center.z() + size) && (data.z() > center.z() - size);
				std::cout << "└-Point: " << data.transpose() << " is inside: " << isInside << std::endl;
			}
		}
		else // Non-leaf node
		{
			indent(depth);
			std::cout << "Node: " << center.transpose() << " " << size << std::endl;
		}

		// Print children
		for (auto &child : children)
			child.print(depth + 1);
	}
};

class Octree
{
private:
	Node root;

public:
	Octree(const PointCloud &pointCloud, Eigen::Vector3d center, double size)
		: root(center, size)
	{
		for (int i = 0; i < pointCloud.points.size(); i++)
			insertPoint(pointCloud.points[i], i);
	}

	// Get all leaf nodes
	std::vector<Node> getLeafNodes()
	{
		std::vector<Node> leafNodes;
		std::queue<Node *> queue;
		queue.push(&root);
		while (!queue.empty())
		{
			Node *node = queue.front();
			queue.pop();
			if (node->isLeaf())
				leafNodes.push_back(*node);
			else
				for (auto &child : node->children)
					queue.push(&child);
		}
		return leafNodes;
	}

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
	// tree.print();

	auto leafNodes = tree.getLeafNodes();
	int numNonEmptyLeafNodes = 0;
	std::cout << "Leaf nodes: " << leafNodes.size() << std::endl;
	for (auto &node : leafNodes)
	{
		node.print();
		if (!node.isEmpty())
			numNonEmptyLeafNodes++;
	}
	std::cout << "Non-empty leaf nodes: " << numNonEmptyLeafNodes << std::endl;
}