// This file implements the Poisson Surface Reconstruction (PSR) algorithm.

#include <random>
#include <queue>
#include <fstream>

#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#define MC_IMPLEM_ENABLE
#define MC_CPP_USE_DOUBLE_PRECISION
#include <MC.h>

#include "Mesh.h"

// Utils

#define GAUSS_QUADRAQURE_N 7

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

PointCloud Mesh::generatePointCloud()
{
	update_normals();

	PointCloud pointCloud;
	pointCloud.reserve(numVertices());

	for (auto v : vertices())
	{
		pointCloud.points.push_back(eigenPoint(v));
		pointCloud.normals.push_back(eigenNormal(v));
	}

	return pointCloud;
}

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

// Find which child the point belongs to
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

// For a given parent center and half size, find the center of the child of certain index
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

struct PointCloudData
{
	Eigen::Vector3d point;
	Eigen::Vector3d normal;

	PointCloudData(Eigen::Vector3d point, Eigen::Vector3d normal)
		: point(point), normal(normal) {}
};

struct Node
{
	// Node box
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	double halfWidth = 1.0;

	// Tree structure
	std::vector<Node> children;
	int depth = 0;

	// Contained point
	// int index = -1;
	// Eigen::Vector3d point;

	std::vector<PointCloudData> data;

	// Neighbor nodes for trilinear interpolation, size (if not empty) should be 8
	std::vector<Node *> neighbors;
	Eigen::Vector<double, 8> neighborWeights;

	// Support of base function of neighbor functions
	Eigen::Vector3d supportCenter;
	double supportSize;

	// Constructors
	Node() = default;

	Node(const Eigen::Vector3d &center, double halfWidth, int depth = 0)
		: center(center), halfWidth(halfWidth), depth(depth) {}

	// Check if node is a leaf
	bool isLeaf() const { return children.empty(); }

	// Check if node is empty
	// bool isEmpty() const { return index == -1; }
	bool isEmpty() const { return data.empty(); }

	// Split node into 8 children
	// Returns the pointer to the child which is supposed to contain the point
	// void splitNode()
	// {
	// 	// Create children
	// 	double halfSize = size / 2;
	// 	children.reserve(8);
	// 	for (int i = 0; i < 8; i++)
	// 		children.emplace_back(childCenter(center, halfSize, i), halfSize, depth + 1);

	// 	// Move point to children
	// 	children[childIndex(center, point)].insertPoint(point, index);

	// 	// Clear point
	// 	index = -1;
	// }

	void splitNode()
	{
		// Create children
		double childHalfWidth = halfWidth / 2;
		children.reserve(8);
		for (int i = 0; i < 8; i++)
			children.emplace_back(childCenter(center, childHalfWidth, i), childHalfWidth, depth + 1);

		// Move point to children
		while (!data.empty())
		{
			auto pcd = data.back();
			data.pop_back();
			children[childIndex(center, pcd.point)].insertPoint(pcd);
		}
	}

	// Insert point into the octree
	bool insertPoint(const PointCloudData &pcd)
	{
		if (isLeaf()) // Node is a leaf
		{
			if (isEmpty()) // Node is empty
			{
				data.push_back(pcd);
				return true;
			}
			else // Node already contains a point
			{
				// Split node and insert current point to children
				splitNode();

				// Insert new point
				return insertPoint(pcd);
			}
		}
		else // Node is not a leaf, insert point into children
			return children[childIndex(center, pcd.point)].insertPoint(pcd);
	}

	// Find the leaf node containing the point
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
		return Intersection(center, halfWidth * 1.5, other.center, other.halfWidth * 1.5);
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

		x = (x - center) / halfWidth;
		return spline(x.x()) * spline(x.y()) * spline(x.z()) / pow(halfWidth, 3);
	}

	// Gradient of base function for FEM
	Eigen::Vector3d gradBaseFunc(Eigen::Vector3d x)
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

		auto dSpline = [](double t)
		{
			if (abs(t) <= 0.5)
				return -2.0 * t;
			else if (abs(t) <= 1.5)
				return -2.0 * (1.5 - abs(t)) * (t > 0 ? 1 : -1);
			else
				return 0.0;
		};

		x = (x - center) / halfWidth;
		return Eigen::Vector3d(dSpline(x.x()) * spline(x.y()) * spline(x.z()),
							   dSpline(x.y()) * spline(x.x()) * spline(x.z()),
							   dSpline(x.z()) * spline(x.x()) * spline(x.y())) /
			   pow(halfWidth, 4);
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
					Eigen::Vector3d x = (min + max) / 2.0 + ((max - min) / 2.0).cwiseProduct(Eigen::Vector3d(gauss.x[i], gauss.x[j], gauss.x[k]));
					element += gauss.w[i] * gauss.w[j] * gauss.w[k] * f(x);
				}

		element *= (max - min).prod() / 8.0;

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
					// result += neighborWeights[i] * (neighbor->gradBaseFunc(q).dot(point));
					result += neighborWeights[i] * (neighbor->gradBaseFunc(q).dot(data[0].normal));
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
						  << " Width: " << halfWidth
						  << " Depth: " << depth
						  << " EMPTY" << std::endl;
			}
			else // Non-empty leaf node
			{
				indent(depth);
				std::cout << "Node: " << center.transpose()
						  << " Width: " << halfWidth
						  << " Depth: " << depth
						  << std::endl;
				indent(depth);
				bool isInside = (center - data[0].point).cwiseAbs().maxCoeff() < halfWidth;
				std::cout << "└-Point: " << data[0].point.transpose()
						  << " is inside: " << isInside
						  << std::endl;

				printNeighbors();
			}
		}
		else // Non-leaf node
		{
			indent(depth);
			std::cout << "Node: " << center.transpose()
					  << " Width: " << halfWidth
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
					  << " Width: " << neighbor->halfWidth
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

	Eigen::VectorXd x;

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

	Octree(const PointCloud &pointCloud, Eigen::Vector3d center, double halfWidth)
		: root(center, halfWidth)
	{
		for (int i = 0; i < pointCloud.points.size(); i++)
			insertPoint(PointCloudData(pointCloud.points[i], pointCloud.normals[i]));
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
		Eigen::Vector3d c = (node.data[0].point / (2 * node.halfWidth)).array().round() * (2 * node.halfWidth);
		double s = node.halfWidth * 2;

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
			auto v = (node.data[0].point - targetPoint) / (node.halfWidth * 2);
			node.neighborWeights[i] = abs(v.x() * v.y() * v.z());

			// Support of base function of neighbor functions
			node.supportCenter = node.center;
			node.supportSize = node.halfWidth * 3.5;
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
		updateNonEmptyLeafNodes();
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

	// Load vector
	Eigen::SparseVector<double> loadVector()
	{
		Eigen::SparseVector<double> b(leafNodes.size());

		for (int i = 0; i < leafNodes.size(); i++)
		{
			auto &node = *leafNodes[i];

			// Integrate (base function * div V) over the support of the base function
			double element = 0.0;
			for (auto &sample : nonEmptyLeafNodes)
			{
				double element2 = 0.0;
				Intersection Intersection(node.center, node.halfWidth * 1.5,
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
								element2 += gauss.w[i] * gauss.w[j] * gauss.w[k] * f(x);
							}

					element2 *= (max - min).prod() / 8.0;
				}

				element += element2;
			}

			b.insert(i) = element;
		}
		return b;
	}

	// TODO: Indicator function
	double indicatorFunction(Eigen::Vector3d q)
	{
		double result = 0.0;
		for (int i = 0; i < leafNodes.size(); i++)
			// if ((leafNodes[i]->center - q).cwiseAbs().maxCoeff() < leafNodes[i]->size * 1.5)
			// if (abs((leafNodes[i]->center - q).x()) < leafNodes[i]->size * 1.5 &&
			// 	abs((leafNodes[i]->center - q).y()) < leafNodes[i]->size * 1.5 &&
			// 	abs((leafNodes[i]->center - q).z()) < leafNodes[i]->size * 1.5)
			result += leafNodes[i]->baseFunc(q) * x[i];
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
	// void insertPoint(Eigen::Vector3d point, int index)
	// {
	// 	root.insertPoint(point, index);
	// }

	void insertPoint(PointCloudData pcd)
	{
		root.insertPoint(pcd);
	}
};

void testOctree(PointCloud pointCloud)
{
	std::cout
		<< "Build octree" << std::endl;
	Octree tree(pointCloud, Eigen::Vector3d::Zero(), 1.0);
	// tree.printDepthFirst();

	tree.updateLeafNodes();
	tree.updateNonEmptyLeafNodes();
	std::cout << "Leaf status: " << tree.leafNodes.size() << " " << tree.nonEmptyLeafNodes.size() << std::endl;

	std::cout << "Refine octree" << std::endl;
	tree.refine();

	std::cout << "Leaf status: " << tree.leafNodes.size() << " " << tree.nonEmptyLeafNodes.size() << std::endl;

	// tree.printBreadthFirst();

	std::cout << "Calculate stiffness matrix" << std::endl;

	tree.updateLeafNodes();
	tree.updateNonEmptyLeafNodes();

	auto L = tree.stiffnessMatrix();
	std::cout << "Non-zero elements: " << L.nonZeros() << std::endl;

	std::cout << "Calculate load vector" << std::endl;
	auto b = tree.loadVector();

	std::cout << "Solve linear system" << std::endl;
	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
	solver.compute(L);

	tree.x = solver.solve(b);

	std::cout << b.transpose() << std::endl;
	std::cout << tree.x.transpose() << std::endl;

	if (solver.info() != Eigen::Success)
	{
		// Handle the case when the solution process fails.
		std::cout << "Failed to solve the linear system" << std::endl;
		return;
	}
	else
	{
		std::cout << "Solved" << std::endl;
	}

	std::cout << "Generate mesh with Marching Cubes" << std::endl;

	MC::mcMesh mesh;

	int n = pow(2, tree.maxDepth);
	double *field = new double[n * n * n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
			{
				Eigen::Vector3d q = (Eigen::Vector3d(i, j, k) / n) * 2 - Eigen::Vector3d::Ones();
				field[i * n * n + j * n + k] = tree.indicatorFunction(q);
				std::cout << field[i * n * n + j * n + k] << " ";
			}

	std::cout << "Constructed field" << std::endl;

	MC::marching_cube(field, n, n, n, mesh);

	std::ofstream out;
	out.open("test.obj");
	if (out.is_open() == false)
		return;
	out << "g "
		<< "Obj" << std::endl;
	for (size_t i = 0; i < mesh.vertices.size(); i++)
		out << "v " << mesh.vertices.at(i).x << " " << mesh.vertices.at(i).y << " " << mesh.vertices.at(i).z << '\n';
	for (size_t i = 0; i < mesh.vertices.size(); i++)
		out << "vn " << mesh.normals.at(i).x << " " << mesh.normals.at(i).y << " " << mesh.normals.at(i).z << '\n';
	for (size_t i = 0; i < mesh.indices.size(); i += 3)
	{
		out << "f " << mesh.indices.at(i) + 1 << "//" << mesh.indices.at(i) + 1
			<< " " << mesh.indices.at(i + 1) + 1 << "//" << mesh.indices.at(i + 1) + 1
			<< " " << mesh.indices.at(i + 2) + 1 << "//" << mesh.indices.at(i + 2) + 1
			<< '\n';
	}
	out.close();

	std::cout << "Output mesh to test.obj" << std::endl;

	// spy(L, "output1.png");
	// spy(L, "output2.png", CONTINUOUS);

	// std::cout << L << std::endl;
}