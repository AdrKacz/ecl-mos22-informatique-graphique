#ifndef DEF_TRIANGLEMESH
#define DEF_TRIANGLEMESH

#define USE_BVH false

#include <string.h>
#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <list>

#include "TriangleIndices.h"
#include "../Ray/Ray.h"
#include "../Vector/Vector.h"
#include "../Object/Object.h"

class BoundingBox
{
public:
	Vector m, M;
	BoundingBox();
	BoundingBox(const Vector&, const Vector&);
	~BoundingBox();
};

class Node
{
public:
	BoundingBox box;

	Node* left;
	Node* right;

	unsigned int i1;
	unsigned int i2;
	
	Node();
	~Node();
};


class TriangleMesh : public Object
{
public:
	TriangleMesh();
	~TriangleMesh();

	Node root_box;
	
	bool intersect_with_triangle(const TriangleIndices&, const Ray&);
	bool intersect_with_triangle(const TriangleIndices&, const Ray&, Vector&, Vector&, double&);
	bool intersect_with_bounding_box(const Ray&, const BoundingBox&);

	virtual bool intersect(const Ray&);
    virtual bool intersect(const Ray&, Vector&, Vector&);
    virtual bool intersect(const Ray&, Vector&, Vector&, double&);
	bool intersect_bounding_volume_hierarchy(const Ray&, Vector&, Vector&, double&);

	double center_of_triangle(unsigned int, unsigned int);

	void readOBJ(const char* obj);
	void init_bounding_box();
	void build_bounding_volume_hierarchy(Node*, unsigned int, unsigned int);
	BoundingBox create_bounding_box(int, int);

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	
};

#endif