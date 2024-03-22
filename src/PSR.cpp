// This file implements the Poisson Surface Reconstruction (PSR) algorithm.

#include <random>
#include <queue>

#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "Mesh.h"

#define GAUSS_QUADRAQURE_N 7

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

/*

		Index of children

		  4-------------------6
		/ |                 / |
	  /   |       back    /   |
	/     |             /     |       z
  /       |           /       |       |
5-------------------7         |       O--y
|         |         |         |      /
|       front       |         |     x
|         |         |         |
|         ０---------|--------2
|       /           |       /
|     /             |     /
|   /               |   /
| /                 | /
1-------------------3

*/

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

class GaussianQuadrature
{
public:
	std::vector<double> x;
	std::vector<double> w;

	GaussianQuadrature(int n = 5)
	{
		switch (n)
		{
		case 1:
			x = {0.0};
			w = {2.0};
			break;
		case 2:
			x = {-0.5773502691896257, 0.5773502691896257};
			w = {1.0, 1.0};
			break;
		case 3:
			x = {-0.7745966692414834, 0.0, 0.7745966692414834};
			w = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
			break;
		case 4:
			x = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
			w = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
			break;
		case 5:
			x = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640};
			w = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891};
			break;
		case 6:
			x = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.9324695142031521};
			w = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704};
			break;
		case 7:
			x = {-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585};
			w = {0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.2797053914892766, 0.1294849661688697};
			break;
		case 8:
			x = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
			w = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
			break;
		case 9:
			x = {-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0.0, 0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261};
			w = {0.0812743883615744, 0.1806481606948574, 0.2606106964029354, 0.3123470770400029, 0.3302393550012598, 0.3123470770400029, 0.2606106964029354, 0.1806481606948574, 0.0812743883615744};
			break;
		case 10:
			x = {-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312, 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717};
			w = {0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963, 0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2190863625159820, 0.1494513491505806, 0.0666713443086881};
			break;
		}
	}
};

struct Intersection
{
	bool isIntersecting = false;
	Eigen::Vector3d min;
	Eigen::Vector3d max;

	Intersection(Eigen::Vector3d leftCenter, double leftRadius,
				 Eigen::Vector3d rightCenter, double rightRadius)
	{
		min = (leftCenter.array() - leftRadius).cwiseMax(rightCenter.array() - rightRadius);
		max = (leftCenter.array() + leftRadius).cwiseMin(rightCenter.array() + rightRadius);
		isIntersecting = (min.array() < max.array()).all();
	}
};

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
	Eigen::Vector<double, 8> neighborWeights;

	// Support of base function of neighbor functions
	Eigen::Vector3d supportCenter;
	double supportSize;

	// Constructors
	Node() = default;

	Node(const Eigen::Vector3d &center, double size, int depth = 0)
		: center(center), size(size), depth(depth) {}

	// Check if node is a leaf
	bool isLeaf() const { return children.empty(); }

	// Check if node is empty
	bool isEmpty() const { return index == -1; }

	// Split node into 8 children
	// Returns the pointer to the child which is supposed to contain the point
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
	Node *findLeaf(const Eigen::Vector3d point)
	{
		if (isLeaf())
			return this;
		else
			return children[childIndex(center, point)].findLeaf(point);
	}

	// Intersection of support of base function of two nodes
	Intersection operator&(const Node &other)
	{
		// Intersection intersection;
		// Eigen::Vector3d min = center.array() - size * 1.5;
		// Eigen::Vector3d max = center.array() + size * 1.5;
		// Eigen::Vector3d otherMin = other.center.array() - other.size * 1.5;
		// Eigen::Vector3d otherMax = other.center.array() + other.size * 1.5;

		// intersection.min = min.cwiseMax(otherMin);
		// intersection.max = max.cwiseMin(otherMax);

		// if ((intersection.min.array() < intersection.max.array()).all())
		// 	intersection.isIntersecting = true;

		return Intersection(center, size * 1.5, other.center, other.size * 1.5);
	}

	// Base function for FEM
	double baseFunc(Eigen::Vector3d x)
	{
		assert(isLeaf());

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

	// Gradient of base function for FEM
	Eigen::Vector3d gradBaseFunc(Eigen::Vector3d x)
	{
		assert(isLeaf());

		auto dSpline = [](double t)
		{
			if (abs(t) <= 0.5)
				return -2.0 * t;
			else if (abs(t) <= 1.5)
				return -2.0 * (1.5 - abs(t)) * (t > 0 ? 1 : -1);
			else
				return 0.0;
		};

		x = (x - center) / size;
		return Eigen::Vector3d(dSpline(x.x()), dSpline(x.y()), dSpline(x.z())) / size;
	}

	// Stiffness matrix element
	double stiffnessMatrixElement(Node &other, Intersection &intersection)
	{
		if (!intersection.isIntersecting)
			return 0.0;

		// Integrate base function product over the intersection
		Eigen::Vector3d min = intersection.min;
		Eigen::Vector3d max = intersection.max;

		GaussianQuadrature gauss(GAUSS_QUADRAQURE_N);

		auto f = [&](Eigen::Vector3d &x)
		{
			return gradBaseFunc(x).dot(other.gradBaseFunc(x));
		};

		double element = 0.0;
		for (int i = 0; i < GAUSS_QUADRAQURE_N; i++)
			for (int j = 0; j < GAUSS_QUADRAQURE_N; j++)
				for (int k = 0; k < GAUSS_QUADRAQURE_N; k++)
				{
					Eigen::Vector3d x = (min + max) + (max - min).cwiseProduct(Eigen::Vector3d(gauss.x[i], gauss.x[j], gauss.x[k]));
					x /= 2.0;
					element += gauss.w[i] * gauss.w[j] * gauss.w[k] * f(x);
				}

		element *= (max.x() - min.x()) * (max.y() - min.y()) * (max.z() - min.z());

		return element;
	}

	// Divergence of V
	// V is the smoothed normal field
	double divV(Eigen::Vector3d q)
	{
		double result = 0.0;
		if (abs((supportCenter - q).array()).maxCoeff() < supportSize)
			for (int i = 0; i < 8; i++)
			{
				auto neighbor = neighbors[i];
				if (abs((neighbor->supportCenter - q).array()).maxCoeff() < neighbor->supportSize)
					result += neighborWeights[i] * (neighbor->gradBaseFunc(q).dot(data));
			}
		return result;
	}

	// Traverse the octree by depth first
	void depthFirst(std::function<void(Node &)> callback)
	{
		callback(*this);
		if (!isLeaf())
			for (auto &child : children)
				child.depthFirst(callback);
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
				std::cout << "Node: " << center.transpose()
						  << " Size: " << size
						  << " Depth: " << depth
						  << " EMPTY" << std::endl;
			}
			else // Non-empty leaf node
			{
				indent(depth);
				std::cout << "Node: " << center.transpose()
						  << " Size: " << size
						  << " Depth: " << depth
						  << std::endl;
				indent(depth);
				bool isInside = (data.x() < center.x() + size) && (data.x() > center.x() - size) &&
								(data.y() < center.y() + size) && (data.y() > center.y() - size) &&
								(data.z() < center.z() + size) && (data.z() > center.z() - size);
				std::cout << "└-Point: " << data.transpose()
						  << " Index: " << index
						  << " is inside: " << isInside
						  << std::endl;

				printNeighbors();
			}
		}
		else // Non-leaf node
		{
			indent(depth);
			std::cout << "Node: " << center.transpose()
					  << " Size: " << size
					  << " Depth: " << depth
					  << std::endl;
		}
	}

	void printNeighbors()
	{
		// Indent
		auto indent = [](int depth)
		{
			for (int i = 0; i < depth; i++)
				std::cout << "    ";
		};

		double sumOfWeights = neighborWeights.sum();
		indent(depth);
		std::cout << "  Sum of weights: " << sumOfWeights << std::endl;

		for (auto &neighbor : neighbors)
		{
			indent(depth);
			std::cout << "  Neighbor: " << neighbor->center.transpose()
					  << " Size: " << neighbor->size
					  << " Depth: " << neighbor->depth
					  << " Distance: " << (neighbor->center - center).norm()
					  << std::endl;
		}
	}
};

struct Octree
{
	Node root;

	// Trade time with space
	int maxDepth = 0;
	std::vector<Node *> leafNodes;
	std::vector<Node *> nonEmptyLeafNodes;

	void updateMaxDepth()
	{
		maxDepth = 0;
		int depth = 0;
		breadthFirst([&depth](Node &node)
					 {
				if (node.depth > depth)
					depth = node.depth; });
		maxDepth = depth;
	}

	void updateLeafNodes()
	{
		leafNodes.clear();
		breadthFirst([&](Node &node)
					 {
				if (node.isLeaf())
					leafNodes.push_back(&node); });
	}

	void updateNonEmptyLeafNodes()
	{
		nonEmptyLeafNodes.clear();
		breadthFirst([&](Node &node)
					 {
				if (node.isLeaf() && !node.isEmpty())
					nonEmptyLeafNodes.push_back(&node); });
	}

	// Traverse the octree by depth first
	void depthFirst(std::function<void(Node &)> callback)
	{
		root.depthFirst(callback);
	}

	// Traverse the octree by breadth first
	void breadthFirst(std::function<void(Node &)> callback)
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

	Octree(const PointCloud &pointCloud, Eigen::Vector3d center, double size)
		: root(center, size)
	{
		for (int i = 0; i < pointCloud.points.size(); i++)
			insertPoint(pointCloud.points[i], i);
	}

	// Split nodes until all neighbors of non-empty leaf nodes are at the maximum depth
	// According to the original paper, the neighbors of a node
	// are 8 nodes that are closest to the sample point of the node.
	void refineNeighbors(Node &node)
	{
		assert(node.isLeaf() && !node.isEmpty() && node.neighbors.empty());

		node.neighbors.reserve(8);

		// Imagine a virtual parent 'node', that contains the neighbors of the current node.
		// Note that this virtual node is very likely not a real node in the octree.
		// The virtual node is used to find the neighbors of the current node.
		// We construct it manually and use the center of its children to find the neighbors.
		Eigen::Vector3d c = (node.data / (2 * node.size)).array().round() * (2 * node.size);
		double s = node.size * 2;
		Node virtualNode(c, s);

		for (int i = 0; i < 8; i++)
		{
			// Use the center of the children of the virtual node to find the neighbors
			Eigen::Vector3d targetPoint = childCenter(c, s / 2, i);
			Node *leaf = root.findLeaf(targetPoint);

			// Split nodes until the neighbor node is at the maximum depth
			while (leaf->depth < maxDepth)
			{
				leaf->splitNode();
				leaf = leaf->findLeaf(targetPoint);
			}

			// Add the neighbor node to the neighbors of the current node
			node.neighbors.push_back(leaf);

			// Compute and save trilinear weights
			auto v = (node.data - targetPoint) / (node.size * 2);
			node.neighborWeights[i] = abs(v.x() * v.y() * v.z());

			// Support of base function of neighbor functions
			node.supportCenter = node.center;
			node.supportSize = node.size * 1.5;
		}
	}

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

		// Split nodes until all neighbors of non-empty leaf nodes are at the maximum depth
		updateNonEmptyLeafNodes();
		for (auto &node : nonEmptyLeafNodes)
			refineNeighbors(*node);

		updateLeafNodes();
	}

	// Stiffness matrix
	// Note that the leaf nodes are sorted by Breadth-first order,
	// so larger elements are more likely to be at the start of the matrix
	Eigen::SparseMatrix<double> stiffnessMatrix()
	{
		int n = leafNodes.size();
		Eigen::SparseMatrix<double> L(n, n);
		L.reserve(50 * n);

		for (int i = 0; i < n; i++)
			for (int j = 0; j <= i; j++)
			{
				auto intersection = *leafNodes[i] & *leafNodes[j];
				if (intersection.isIntersecting)
				{
					double element = leafNodes[i]->stiffnessMatrixElement(*leafNodes[j], intersection);
					L.insert(i, j) = element;
					if (i != j)
						L.insert(j, i) = element;
				}
			}

		L.makeCompressed();

		return L;
	}

	// // Divergence of V
	// // V is the smoothed normal field
	// double divV(Eigen::Vector3d q)
	// {
	// 	double result = 0.0;
	// 	for (auto &node : nonEmptyLeafNodes)
	// 		if (abs((node->supportCenter - q).array()).maxCoeff() < node->supportSize)
	// 			for (int i = 0; i < 8; i++)
	// 			{
	// 				auto neighbor = node->neighbors[i];
	// 				if (abs((neighbor->supportCenter - q).array()).maxCoeff() < neighbor->supportSize)
	// 					result += node->neighborWeights[i] * (neighbor->gradBaseFunc(q).dot(node->data));
	// 			}
	// 	return result;
	// }

	// Load vector
	Eigen::SparseVector<double> loadVector()
	{
		Eigen::SparseVector<double> b(leafNodes.size());

		for (int i = 0; i < leafNodes.size(); i++)
		{
			auto &node = *leafNodes[i];

			// Integrate (base function * divV) over the support of the base function
			double element = 0.0;
			for (auto &sample : nonEmptyLeafNodes)
			{
				Intersection Intersection(node.center, node.size * 1.5,
										  sample->supportCenter, sample->supportSize);
				if (Intersection.isIntersecting)
				{
					Eigen::Vector3d min = Intersection.min;
					Eigen::Vector3d max = Intersection.max;

					GaussianQuadrature gauss(GAUSS_QUADRAQURE_N);

					auto f = [&](Eigen::Vector3d &x)
					{
						return node.baseFunc(x) * sample->divV(x);
					};

					for (int i = 0; i < GAUSS_QUADRAQURE_N; i++)
						for (int j = 0; j < GAUSS_QUADRAQURE_N; j++)
							for (int k = 0; k < GAUSS_QUADRAQURE_N; k++)
							{
								Eigen::Vector3d x = (min + max) + (max - min).cwiseProduct(Eigen::Vector3d(gauss.x[i], gauss.x[j], gauss.x[k]));
								x /= 2.0;
								element += gauss.w[i] * gauss.w[j] * gauss.w[k] * f(x);
							}
				}
			}

			b.insert(i) = element;
		}
		return b;
	}

	// TODO: Indicator function
	double indicatorFunction(Eigen::Vector3d q)
	{
		double result = 0.0;
		for (auto &node : nonEmptyLeafNodes)
			if ((node->center - q).cwiseAbs().maxCoeff() < node->size * 1.5)
				result += node->baseFunc(q) * x[i];
		return result;
	}

	// Debug
	void printDepthFirst()
	{
		depthFirst([](Node &node)
				   { node.print(); });
	}

	void printBreadthFirst()
	{
		breadthFirst([](Node &node)
					 { node.print(); });
	}

private:
	void insertPoint(Eigen::Vector3d point, int index)
	{
		root.insertPoint(point, index);
	}
};

void testGaussianQuadrature(int n)
{
	std::cout << "Testing Gauss Quadrature with " << n << " points" << std::endl;

	GaussianQuadrature gauss(n);
	auto x = gauss.x;
	auto w = gauss.w;

	// Accurate results
	double sinAccurate = 0.0;							 // Accurate result of integrating sin(x) from -1 to 1
	double cosAccurate = std::sin(1.0) - std::sin(-1.0); // Accurate result of integrating cos(x) from -1 to 1
	double xAccurate = 0.0;								 // Accurate result of integrating x from -1 to 1
	double x2Accurate = 2.0 / 3.0;						 // Accurate result of integrating x^2 from -1 to 1
	double x3Accurate = 0.0;							 // Accurate result of integrating x^3 from -1 to 1
	double expAccurate = std::exp(1.0) - std::exp(-1.0); // Accurate result of integrating e^x from -1 to 1

	// Test with annoying function
	auto test = [&x, &w](std::function<double(double)> f, double accurate)
	{
		double result = 0.0;
		for (int i = 0; i < x.size(); i++)
			result += w[i] * f(x[i]);
		double logError = std::log(std::abs(result - accurate)) / std::log(10);
		std::cout << "Result: " << result << " (Error: " << logError << ")" << std::endl;
	};

	// Test with sin
	test([](double x)
		 { return std::sin(x); },
		 sinAccurate);

	// Test with cos
	test([](double x)
		 { return std::cos(x); },
		 cosAccurate);

	// Test with x
	test([](double x)
		 { return x; },
		 xAccurate);

	// Test with x^2
	test([](double x)
		 { return x * x; },
		 x2Accurate);

	// Test with x^3
	test([](double x)
		 { return x * x * x; },
		 x3Accurate);

	// Test with e^x
	test([](double x)
		 { return std::exp(x); },
		 expAccurate);
}

// Spy a sparse matrix by writing it to a PNG file

// #define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

enum SpyType
{
	CONTINUOUS,
	DISCRETE,
};

void spy(Eigen::SparseMatrix<double> matrix, std::string filename, SpyType spyType = DISCRETE)
{
	int width = matrix.cols();
	int height = matrix.rows();

	std::vector<unsigned char> image(width * height, 255);

	double max = 0.0;
	for (int k = 0; k < matrix.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it)
		{
			max = std::max(max, it.value());
		}

	for (int k = 0; k < matrix.outerSize(); ++k)
		for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, k); it; ++it)
		{
			int i = it.row();
			int j = it.col();
			double value = it.value();
			switch (spyType)
			{
			case CONTINUOUS:
				image[i * width + j] = (unsigned char)(255 - value / max * 255);
				break;
			case DISCRETE:
				image[i * width + j] = 0;
				break;
			}
		}

	stbi_flip_vertically_on_write(true);
	stbi_write_png(filename.data(), width, height, 1, image.data(), 0);
}

void testOctree(PointCloud pointCloud)
{

	// for (int n : {1, 2, 3, 4, 5, 6, 7, 8, 9, 10})
	// {
	// 	testGaussianQuadrature(n);
	// }

	Octree tree(pointCloud, Eigen::Vector3d::Zero(), 1.0);
	// tree.printDepthFirst();

	tree.refine();

	auto L = tree.stiffnessMatrix();
	auto b = tree.loadVector();

	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
	solver.compute(L);
	Eigen::VectorXd x = solver.solve(b);

	std::cout << x.transpose() << std::endl;

	// spy(L, "output1.png");
	// spy(L, "output2.png", CONTINUOUS);

	// std::cout << L << std::endl;
}