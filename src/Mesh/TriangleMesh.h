#ifndef DEF_TRIANGLEMESH
#define DEF_TRIANGLEMESH

#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>

#include "TriangleIndices.h"
#include "../Vector/Vector.h"

class TriangleMesh {
public:
	~TriangleMesh();
	TriangleMesh();
	
	bool intersect_with_triangle(const TriangleIndices&, const Vector&, const Vector&);
	void readOBJ(const char* obj);

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;
	
};

#endif