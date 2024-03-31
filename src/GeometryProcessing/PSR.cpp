// This file implements the Poisson Surface Reconstruction (PSR) algorithm.

#include <random>
#include <queue>
#include <fstream>

#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>

#define MC_IMPLEM_ENABLE
#define MC_CPP_USE_DOUBLE_PRECISION
#include <MarchingCube/MC.h>

#include "GeometryProcessing.h"

// Utils

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

double integrate3D(Eigen::Vector3d min, Eigen::Vector3d max,
				   std::function<double(const Eigen::Vector3d &)> f,
				   GaussianQuadrature gauss = GaussianQuadrature(5))
{
	double result = 0.0;
	int n = gauss.x.size();

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
			{
				Eigen::Vector3d x = (min + max) / 2.0 + ((max - min) / 2.0).cwiseProduct(Eigen::Vector3d(gauss.x[i], gauss.x[j], gauss.x[k]));
				result += gauss.w[i] * gauss.w[j] * gauss.w[k] * f(x);
			}

	result *= (max - min).prod() / 8.0;

	return result;
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

// Spline function and its derivative
double spline(double t)
{
	if (abs(t) <= 0.5)
		return 0.75 - t * t;
	else if (abs(t) <= 1.5)
		return 0.5 * pow(1.5 - abs(t), 2);
	else
		return 0.0;
}

double dSpline(double t)
{
	if (abs(t) <= 0.5)
		return -2.0 * t;
	else if (abs(t) <= 1.5)
		return -2.0 * (1.5 - abs(t)) * (t > 0 ? 1 : -1);
	else
		return 0.0;
}

struct Node
{
	// Node box
	Eigen::Vector3d center = Eigen::Vector3d::Zero();
	double halfWidth = 1.0;

	// Tree structure
	std::vector<Node> children;
	int depth = 0;

	int leafIndex = -1;

	// Contained point
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
	// Used in calculating the indicator function
	double baseFunc(Eigen::Vector3d x)
	{
		assert(isLeaf());

		x = (x - center) / halfWidth;
		return spline(x.x()) * spline(x.y()) * spline(x.z()) / pow(halfWidth, 3);
	}

	// Gradient of base function for FEM
	// Used in construncting the stiffness matrix and the load vector
	Eigen::Vector3d gradBaseFunc(Eigen::Vector3d x)
	{
		assert(isLeaf());

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

		auto f = [&](const Eigen::Vector3d &x)
		{
			return gradBaseFunc(x).dot(other.gradBaseFunc(x));
		};

		return integrate3D(min, max, f);
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
					result += neighborWeights[i] * (neighbor->gradBaseFunc(q).dot(data[0].normal));
			}
		return result;
	}

	// Traverse the octree by depth first
	void depthFirst(
		std::function<void(Node &)> callback,
		std::function<bool(Node &)> condition = [](Node &node)
		{ return true; })
	{
		if (condition(*this))
		{
			callback(*this);
			for (auto &child : children)
				child.depthFirst(callback, condition);
		}
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
		int leafIndex = 0;
		breadthFirst([&](Node &node)
					 {
				if (node.isLeaf())
				{
					leafNodes.push_back(&node); 
					node.leafIndex = leafIndex++;
				} });
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
	void depthFirst(
		std::function<void(Node &)> callback,
		std::function<bool(Node &)> condition = [](Node &node)
		{ return true; })
	{
		root.depthFirst(callback, condition);
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
		for (int i = 0; i < pointCloud.data.size(); i++)
			insertPoint(PointCloudData(pointCloud.data[i].point, pointCloud.data[i].normal));
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
		{
			auto supportOverlap = [&](Node &node)
			{
				return (*leafNodes[i] & node).isIntersecting;
			};
			auto calculateElement = [&](Node &node)
			{
				if (node.isLeaf())
				{
					// double element = leafNodes[i]->stiffnessMatrixElement(node, *leafNodes[i] & node);
					L.insert(i, node.leafIndex) = leafNodes[i]->stiffnessMatrixElement(node, *leafNodes[i] & node);
				}
			};
			depthFirst(calculateElement, supportOverlap);
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

					auto f = [&](const Eigen::Vector3d &x)
					{
						return node.baseFunc(x) * sample->divV(x);
					};

					element += integrate3D(min, max, f);
				}
			}

			b.insert(i) = element;
		}
		return b;
	}

	// Indicator function
	double indicatorFunction(Eigen::Vector3d q)
	{
		double result = 0.0;

		auto isInSupport = [&](Node &node)
		{
			return abs((node.center - q).array()).maxCoeff() < node.halfWidth * 1.5;
		};
		auto calculateIndicatorFunc = [&](Node &node)
		{
			if (node.isLeaf())
				result += node.baseFunc(q) * x[node.leafIndex];
		};

		depthFirst(calculateIndicatorFunc, isInSupport);

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
	void insertPoint(PointCloudData pcd)
	{
		root.insertPoint(pcd);
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
};

void poissonSurfaceReconstruction(PointCloud &pointCloud)
{
	// Build and refine octree
	std::cout << "Build octree" << std::endl;
	Octree tree(pointCloud, Eigen::Vector3d::Zero(), 1.0);

	tree.updateLeafNodes();
	tree.updateNonEmptyLeafNodes();
	std::cout << "Leaf status: " << tree.leafNodes.size() << " " << tree.nonEmptyLeafNodes.size() << std::endl;

	std::cout << "Refine octree" << std::endl;
	tree.refine();

	std::cout << "Leaf status: " << tree.leafNodes.size() << " " << tree.nonEmptyLeafNodes.size() << std::endl;

	// FEM
	std::cout << "Calculate stiffness matrix" << std::endl;

	auto L = tree.stiffnessMatrix();
	std::cout << "Non-zero elements: " << L.nonZeros() << std::endl;

	std::cout << "Calculate load vector" << std::endl;
	auto b = tree.loadVector();

	std::cout << "Solve linear system" << std::endl;

	// Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double>> solver;
	solver.compute(L);

	tree.x = solver.solve(b);

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

	// Marching Cubes
	std::cout << "Construct indicator function field" << std::endl;
	int n = 100;
	double *field = new double[n * n * n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
			{
				Eigen::Vector3d q = (Eigen::Vector3d(i, j, k) / n) * 2 - Eigen::Vector3d::Ones();
				field[i * n * n + j * n + k] = tree.indicatorFunction(q);
			}

	std::cout << "Generate mesh with Marching Cubes" << std::endl;
	MC::mcMesh mesh;
	MC::marching_cube(field, n, n, n, mesh);

	// Output mesh to file
	std::ofstream out;
	out.open("output.obj");
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

	std::cout << "Output mesh to output.obj" << std::endl;
}