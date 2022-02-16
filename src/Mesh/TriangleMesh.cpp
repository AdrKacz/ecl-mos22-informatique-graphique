#include "TriangleMesh.h"

#include <iostream>
#include <string>

BoundingBox::BoundingBox()
{
}

BoundingBox::BoundingBox(const Vector& p_m, const Vector& p_M)
{
	m = p_m;
	M = p_M;
}

BoundingBox::~BoundingBox()
{
}

TriangleMesh::TriangleMesh() 
{
}

TriangleMesh::~TriangleMesh()
{
}

bool TriangleMesh::intersect_with_triangle(const TriangleIndices& triangle, const Ray& ray)
{
	double t;
	Vector p, n;
	return intersect_with_triangle(triangle, ray, p, n, t);

}

bool TriangleMesh::intersect_with_triangle(const TriangleIndices& triangle, const Ray& ray, Vector& P, Vector& N, double& t)
{
	Vector vertice_i = vertices[triangle.vtxi];
	Vector vertice_j = vertices[triangle.vtxj];
	Vector vertice_k = vertices[triangle.vtxk];

	Vector e_1 = vertice_j - vertice_i;
	Vector e_2 = vertice_k - vertice_i;

	N = e_1.cross(e_2);

	Vector i_to_origin = (ray.C - vertice_i);
	Vector i_to_origin_cross_u = i_to_origin.cross(ray.u);
	double u_dot_n = ray.u.dot(N);

	double beta  = - e_2.dot(i_to_origin_cross_u) / u_dot_n;
	double gamma  = + e_1.dot(i_to_origin_cross_u) / u_dot_n;
	double alpha = 1. - beta - gamma;
	t = - i_to_origin.dot(N) / u_dot_n;
	if (!(0 <= alpha && alpha <= 1 && 0 <= beta && beta <= 1 && 0 <= gamma && gamma <= 1 && t >= 0)) {
		return false;
	}
	P = vertice_i + e_1 * beta + e_2 * gamma;
	N.normalize();
	return true;
}

bool TriangleMesh::intersect_with_bounding_box(const Ray& r)
{
	// std::string outputm = std::string("m> " + box.m.to_string() + "\n");
	// std::string outputM = std::string("M> " + box.M.to_string() + "\n--- --- ---\n");
	// std::cout << std::string(outputm + outputM);

	double tx_alpha = (box.m[0] - r.C[0]) / r.u[0];
	double tx_beta = (box.M[0] - r.C[0]) / r.u[0];
	double tx1 = std::min(tx_alpha, tx_beta);
	double tx2 = std::max(tx_alpha, tx_beta);

	double ty_alpha = (box.m[1] - r.C[1]) / r.u[1];
	double ty_beta = (box.M[1] - r.C[1]) / r.u[1];
	double ty1 = std::min(ty_alpha, ty_beta);
	double ty2 = std::max(ty_alpha, ty_beta);

	double tz_alpha = (box.m[2] - r.C[2]) / r.u[2];
	double tz_beta = (box.M[2] - r.C[2]) / r.u[2];
	double tz1 = std::min(tz_alpha, tz_beta);
	double tz2 = std::max(tz_alpha, tz_beta);

	double tmax_min = std::max(tx1, std::max(ty1, tz1));
	if (tmax_min < 0)
	{
		return false;
	}

	double tmin_max = std::min(tx2, std::min(ty2, tz2));
	if (tmax_min <= tmin_max)
	{
		return true;
	}
	return false;
}

bool TriangleMesh::intersect(const Ray& r)
{
	for (int i = 0; i < indices.size(); i++)
	{
		if (intersect_with_triangle(indices[i], r))
		{
			return true;
		}
	}
	return false;
}

bool TriangleMesh::intersect(const Ray& r, Vector& P, Vector& N)
{
    double T;
    return intersect(r, P, N, T);
}

bool TriangleMesh::intersect(const Ray& r, Vector& P, Vector& N, double& T)
{
	if (!intersect_with_bounding_box(r)) {
		return false;
	}
    bool has_intersected = false;
    T = std::numeric_limits<double>::max();
    for (int i = 0; i < indices.size(); i++)
    {
		double t;
		Vector p, n;
		if (intersect_with_triangle(indices[i], r, p, n, t))
		{
			has_intersected = true;
			if (t < T) {
				T = t;
				P = p;
				N = n * +1.;
			}
		}
    }

    return has_intersected; 
}

void TriangleMesh::init_bounding_box() {
	box = create_bounding_box(0, vertices.size());
}

BoundingBox TriangleMesh::create_bounding_box(int from_index, int to_index) {
	double double_min = std::numeric_limits<double>::min();
	double double_max = std::numeric_limits<double>::max();

	Vector m = Vector(double_max, double_max, double_max);
	Vector M = Vector(double_min, double_min, double_min);

	for (unsigned int i = from_index; i < to_index; i++)
	{
		m[0] = std::min(m[0], vertices[i][0]);
		m[1] = std::min(m[1], vertices[i][1]);
		m[2] = std::min(m[2], vertices[i][2]);

		M[0] = std::max(M[0], vertices[i][0]);
		M[1] = std::max(M[1], vertices[i][1]);
		M[2] = std::max(M[2], vertices[i][2]);
	}

	return BoundingBox(m, M);
}

void TriangleMesh::readOBJ(const char* obj) {

	char matfile[255];
	char grp[255];

	FILE* f;
	f = fopen(obj, "r");
	int curGroup = -1;
	while (!feof(f)) {
		char line[255];
		if (!fgets(line, 255, f)) break;

		std::string linetrim(line);
		linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
		strcpy(line, linetrim.c_str());

		if (line[0] == 'u' && line[1] == 's') {
			sscanf(line, "usemtl %[^\n]\n", grp);
			curGroup++;
		}

		if (line[0] == 'v' && line[1] == ' ') {
			Vector vec;

			Vector col;
			if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
				col[0] = std::min(1., std::max(0., col[0]));
				col[1] = std::min(1., std::max(0., col[1]));
				col[2] = std::min(1., std::max(0., col[2]));

				vertices.push_back(vec);
				vertexcolors.push_back(col);

			} else {
				sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				vertices.push_back(vec);
			}
		}
		if (line[0] == 'v' && line[1] == 'n') {
			Vector vec;
			sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
			normals.push_back(vec);
		}
		if (line[0] == 'v' && line[1] == 't') {
			Vector vec;
			sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
			uvs.push_back(vec);
		}
		if (line[0] == 'f') {
			TriangleIndices t;
			int i0, i1, i2, i3;
			int j0, j1, j2, j3;
			int k0, k1, k2, k3;
			int nn;
			t.group = curGroup;

			char* consumedline = line + 1;
			int offset;

			nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
			if (nn == 9) {
				if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
				if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
				if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
				if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
				if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
				if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
				if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
				if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
				if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
				indices.push_back(t);
			} else {
				nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
				if (nn == 6) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
					if (nn == 3) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
						if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
						if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
						indices.push_back(t);
					}
				}
			}

			consumedline = consumedline + offset;

			while (true) {
				if (consumedline[0] == '\n') break;
				if (consumedline[0] == '\0') break;
				nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
				TriangleIndices t2;
				t2.group = curGroup;
				if (nn == 3) {
					if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
					if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
					if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
					if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
					if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
					if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
					if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
					if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
					if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
					indices.push_back(t2);
					consumedline = consumedline + offset;
					i2 = i3;
					j2 = j3;
					k2 = k3;
				} else {
					nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
					if (nn == 2) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						indices.push_back(t2);
					} else {
						nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
							if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
							if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
							consumedline = consumedline + offset;
							i2 = i3;
							k2 = k3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u%n", &i3, &offset);
							if (nn == 1) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								indices.push_back(t2);
							} else {
								consumedline = consumedline + 1;
							}
						}
					}
				}
			}

		}

	}
	fclose(f);

}