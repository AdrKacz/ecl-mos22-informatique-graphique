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

class TriangleMesh : public Object
{
public:
	TriangleMesh();
	~TriangleMesh();

	// Bounding box
	Vector m, M;

	// void computer_bounding_box();
	
	bool intersect_with_triangle(const TriangleIndices&, const Ray&);
	bool intersect_with_triangle(const TriangleIndices&, const Ray&, Vector&, Vector&, double&);
	bool intersect_with_bounding_box(const Ray&);

	virtual bool intersect(const Ray&);
    virtual bool intersect(const Ray&, Vector&, Vector&);
    virtual bool intersect(const Ray&, Vector&, Vector&, double&);

	void readOBJ(const char* obj);
	void init_bounding_box();

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	
};

#endif