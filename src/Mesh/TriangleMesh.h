#ifndef DEF_TRIANGLEMESH
#define DEF_TRIANGLEMESH

#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>

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

class TriangleMesh : public Object
{
public:
	TriangleMesh();
	~TriangleMesh();

	BoundingBox box;
	
	bool intersect_with_triangle(const TriangleIndices&, const Ray&);
	bool intersect_with_triangle(const TriangleIndices&, const Ray&, Vector&, Vector&, double&);
	bool intersect_with_bounding_box(const Ray&);

	virtual bool intersect(const Ray&);
    virtual bool intersect(const Ray&, Vector&, Vector&);
    virtual bool intersect(const Ray&, Vector&, Vector&, double&);

	void readOBJ(const char* obj);
	void init_bounding_box();
	BoundingBox create_bounding_box(int, int);

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	
};

#endif