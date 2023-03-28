#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"
#include<windows.h>
#include<winuser.h>
#include <strstream>
#include <algorithm>
using namespace std;

struct vec2d {
	float u = 0;
	float v = 0;
};

struct vec3d {
	float x = 0;
	float y = 0;
	float z = 0;
	float w = 1;




	vec3d() {
		x = y = z = 0; w = 1;
	}

	vec3d(float a, float b, float c) {
		x = a; y = b; z = c; w = 1;
	}

	float length() const {
		return sqrtf(x * x + y * y + z * z);
	}
	//vect start//
	vec3d normal() {
		return *this / length();
	}

	vec3d Normilize() {
		return *this /= length();
	}

	float dotprod(const vec3d& b) { // a.dotprod(b) returs the dot product of a and b
		return (x * b.x + y * b.y + z * b.z);
	}

	

	vec3d cross(const vec3d& b) {
		return{ this->y * b.z - this->z * b.y,
			this->z * b.x - this->x * b.z,
		this->x * b.y - this->y * b.x };
	}

	vec3d& operator+=(const vec3d& rhs) {
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		return *this;
	}

	vec3d& operator-=(const vec3d& rhs) {
		this->x -= rhs.x;
		this->y -= rhs.y;
		this->z -= rhs.z;
		return *this;
	}

	vec3d& operator*=(const float rhs) {
		this->x *= rhs;
		this->y *= rhs;
		this->z *= rhs;
		return *this;
	}

	vec3d& operator/=(const float rhs) {
		this->x /= rhs;
		this->y /= rhs;
		this->z /= rhs;
		return *this;
	}

	vec3d operator+(const vec3d& rhs) {
		vec3d r = { 0,0,0 };
		r.x =
			this->x + rhs.x;
		r.y =
			this->y + rhs.y;
		r.z =
			this->z + rhs.z;
		return r;
	}

	vec3d operator-(const vec3d& rhs) {
		vec3d r;
		r.x =
			this->x - rhs.x;
		r.y =
			this->y - rhs.y;
		r.z =
			this->z - rhs.z;
		return r;
	}

	vec3d operator*(const float rhs) {
		vec3d r;
		r.x =
			this->x * rhs;
		r.y =
			this->y * rhs;
		r.z =
			this->z * rhs;
		return r;
	}

	vec3d operator/(const float rhs) {
		vec3d r;
		r.x =
			this->x / rhs;
		r.y =
			this->y / rhs;
		r.z =
			this->z / rhs;
		return r;
	}

	bool operator==(const vec3d& rhs) {
		return (
			this->x == rhs.x &&
			this->y == rhs.y &&
			this->z == rhs.z
			);
	}

};


struct triangle {
	vec3d p[3];
	vec2d t[3];
	int criticalPoints = 0;
	float shade;
	float r = 0;
	float g = 0;
	float b = 0;

};



struct mesh {
	vector<triangle> tris;
	vector<vec3d> verts;

	float r = 0;
	float g = 0;
	float b = 0;



	bool loadFromObjectFile(string sFilename)
	{
		ifstream f(sFilename);
		if (!f.is_open())
			return false;



		while (!f.eof()) {
			char line[128];
			f.getline(line, 128);

			strstream s;
			s << line;
			char junk;

			if (line[0] == 'v') {
				vec3d v;
				s >> junk >> v.x >> v.y >> v.z;
				verts.push_back(v);
			}

			if (line[0] == 'f') {
				int f[3];
				s >> junk >> f[0] >> f[1] >> f[2];
				tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
			}
		}






		return true;
	}


	void setColour(float ri, float gi, float bi) {
		r = ri; g = gi; b = bi;
	}

};

struct meshNode { // mesh tree for meshes that have moving parts
	vector<triangle> tri;
	vector<vec3d>criticalPoints;
	vec3d toOrigin; //for rotation mabey not needed (critical point to origin)
	vec3d position; // if top of tree then place in world(top's root is the world root) else then its the displacement needed for its local root
	float Rotx = 0;
	float Roty = 0;
	float Rotz = 0;

	meshNode() {}

	meshNode(vector<triangle> model, float rX, float rY, float rZ, vec3d translation) {
		for (auto& obj : model) {
			tri.push_back(obj);
		}
		Rotx = rX; Roty = rY; Rotz = rZ; position = translation;

	}

	meshNode(vector<triangle> model, float rX, float rY, float rZ, vec3d translation, vector<vec3d>cP) {
		for (auto& obj : model) {
			tri.push_back(obj);
		}
		Rotx = rX; Roty = rY; Rotz = rZ; position = translation;
		for (auto& obj : cP) {
			criticalPoints.push_back(obj);
		}
	}
};

struct meshTree {
	meshNode myMeshNode;
	vector<meshTree> subMeshNodes;
	boolean hasSubs;

	meshTree() {}

	meshTree(meshNode p) {
		myMeshNode = p;
		hasSubs = false;
	}

	meshTree(meshNode p, vector<meshTree> v) {
		myMeshNode = p;
		subMeshNodes = v;
		hasSubs = true;
	}

	bool loadFromMeshTreeFile(string sFilename)
	{
		ifstream f(sFilename);
		if (!f.is_open())
			return false;



		while (!f.eof()) {
			char line[128];
			f.getline(line, 128);

			strstream s;
			s << line;
			char junk;

			if (line[0] == 'n') {
				vec3d v;
				s >> junk >> v.x >> v.y >> v.z;
				//verts.push_back(v);
			}

			if (line[0] == 's') {
				int f[3];
				s >> junk >> f[0] >> f[1] >> f[2];
				//tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
			}
			if (line[0] == 'b') {
				int f[3];
				s >> junk >> f[0] >> f[1] >> f[2];
				//tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
			}
		}






		return true;
	} // TODO

	void addSub(meshTree subTree) {
		hasSubs = true;
		subMeshNodes.push_back(subTree);

	}


};

struct aabb {
	float minx, miny, minz, maxx, maxy, maxz;
	aabb() { minx = 0; miny = 0; minz = 0; maxx = 0; maxy = 0; maxz = 0; };
	aabb(float a, float b, float c, float d, float e, float f) {
		minx = a; miny = b; minz = c; maxx = d; maxy = e; maxz = f;
	}

}; //TODO make better hitboxs

struct entity {
	vec3d position;
	float Rotx = 0;
	float Roty = 0;
	float Rotz = 0;
	float scale = 1;
	vec3d speed = { 0,0,0 };
	bool rasterFirst = false;
	mesh model;
	int hitboxType = 0;
	aabb aabbHitbox;

	entity() {  }
	entity(float a, float b, float c, string filename) {
		position.x = a; position.y = b; position.z = c; model.loadFromObjectFile(filename);
	}

	entity(float a, float b, float c, string filename, float rx, float ry, float rz) {
		position.x = a; position.y = b; position.z = c; model.loadFromObjectFile(filename);
		Rotx = rx; Roty = ry; Rotz = rz;
	}
	void setColour(float a, float b, float c) {
		model.setColour(a, b, c);
	}

	meshNode createNode() {

	}
};

struct cursor {
	vec3d position = { 1,1.5,-1 };
	vec3d normal = { 0,1,1 };

	mesh cursorPoint() {

		//normal = normal.normal();

		vector<triangle> tris;

		entity point(position.x, position.y, position.z, "cube.obj", 0, 0, 0);
		point.setColour(100, 255, 100);
		point.scale = .01;
		point.setColour(100, 255, 100);

		return point.model;

	}


	mesh cursorArrow() {



		vec3d normalPoint = position + (normal);

		entity normalModel(normalPoint.x, normalPoint.y, normalPoint.z, "cube.obj", 0, 0, 0);
		normalModel.scale = .1;


		return normalModel.model;

	}

	entity cursorPlane() {
		float l = normal.length();
		normal = normal / l;

		entity normalModel(position.x, position.y, position.z, "ground.obj", 0, atan2(normal.z, normal.x) + 3.14, acos(normal.y / normal.length()));
		normalModel.scale = -1;
		return normalModel;
	}
	entity cursorInvertedPlane() {

		entity normalModel(position.x, position.y, position.z, "ground.obj", 0, atan2(normal.z, normal.x) + 3.14, acos(normal.y / normal.length()) + 3.14);
		normalModel.scale = 1;
		return normalModel;
	}



};


struct mat4x4
{
	float m[4][4] = { 0 };
};


struct weapon {
	int damage = 50;
	entity fire(vec3d direction, vec3d location) {
		entity e1(location.x, location.y, location.z, "cube.obj", 0, 0, 0);
		e1.setColour(10, 10, 150);
		e1.scale = .1;
		e1.speed = direction;
		return e1;
	}
	entity weapon;
	int type = 1;//0 for full auto, 1 for semi, n for nburst; 
	float deltaTime = .1;// time between events 
	float lastTime = 0; // time since last event 
};

struct player {
	//camera and movement
	vec3d vCamera = { 0,0,0 };
	vec3d vLookDir = { 0,0,1 };
	vec3d speed = { 0,0,0 };
	float fYaw;
	float fPitch;


	entity inHand;
	weapon wInHand;
	weapon primary;
	weapon secondary;
	int weaponSelected = 0;

	void drawWeapon(int weapon) {
		if (weapon == weaponSelected) {
			inHand = entity();
			weaponSelected = 0;
		}
		else {
			switch (weapon)
			{
			case 1:  	cout << "draw weapon one!" << endl;
				inHand = primary.weapon;
				wInHand = primary;
				weaponSelected = 1;
				break;
			case 2: 	cout << "weapon two!" << endl;
				inHand = secondary.weapon;
				wInHand = secondary;
				weaponSelected = 2;
				break;
			case 3: 	cout << "weapon three!" << endl;
				break;
			default: 	cout << "empty hand" << endl;
			}
		}
	}



	entity fire() {
		switch (weaponSelected)
		{
		case 1:  	cout << "weapon one!" << endl;
			return primary.fire(vLookDir, vCamera);
			break;
		case 2: 	cout << "weapon two!" << endl;
			return secondary.fire(vLookDir, vCamera);
			break;
		case 0: 	cout << "weapon three!" << endl;
			break;
		default: 	cout << "empty hand" << endl;
		}

	}

};

// Override base class with your custom functionality
class FPSengine : public olc::PixelGameEngine
{
public:
	FPSengine()
	{
		// Name you application
		sAppName = "FPSGame";
	}

private: // Private variables
	mesh meshCube;
	mat4x4 matProj;
	float fTheta;
	vector<entity>objects;
	vector<meshTree>trees;
	bool lockMouse;
	float sensitivity = 0.005;
	cursor cursor;
	boolean cursorControl;
	entity selectedObject;
	meshTree topNode;
	meshTree currentNode;
	int lastSelected = 0;
	float timeChecker;
	player player;



	
	vec3d Vector_IntersectPlane(vec3d& plane_p, vec3d& plane_n, vec3d& lineStart, vec3d& lineEnd, float& t)
	{
		plane_n = plane_n.normal();
		float plane_d = - plane_n.dotprod(plane_p) ;
		float ad = lineStart.dotprod(plane_n);
		float bd = lineEnd.dotprod(plane_n);
		t = (-plane_d - ad) / (bd - ad);
		vec3d lineStartToEnd = lineEnd - lineStart;
		vec3d lineToIntersect = lineStartToEnd * t;
		return (lineStart + lineToIntersect);
	}


	//matrix methods
	vec3d MatrixMultiplyVector(vec3d& i, mat4x4& m) {
		vec3d v;
		v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
		v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
		v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
		v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
		return v;
	}
	mat4x4 Matrix_MakeIdentity()
	{
		mat4x4 matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		return matrix;

	}
	mat4x4 Matrix_MakeRotationX(float fAngleRad) {
		mat4x4 matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = cosf(fAngleRad);
		matrix.m[1][2] = sinf(fAngleRad);
		matrix.m[2][1] = -sinf(fAngleRad);
		matrix.m[2][2] = cosf(fAngleRad);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}
	mat4x4 Matrix_MakeRotationY(float fAngleRad) {
		mat4x4 matrix;
		matrix.m[0][0] = cosf(fAngleRad);
		matrix.m[0][2] = sinf(fAngleRad);
		matrix.m[2][0] = -sinf(fAngleRad);
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = cosf(fAngleRad);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}
	mat4x4 Matrix_MakeRotationZ(float fAngleRad) {
		mat4x4 matrix;
		matrix.m[0][0] = cosf(fAngleRad);
		matrix.m[0][1] = sinf(fAngleRad);
		matrix.m[1][0] = -sinf(fAngleRad);
		matrix.m[1][1] = cosf(fAngleRad);
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		return matrix;
	}
	mat4x4 Matrix_MakeScale(float a, float b) {
		mat4x4 matrix;
		matrix.m[0][0] = a / b;
		matrix.m[1][1] = a / b;
		matrix.m[2][2] = a / b;
		matrix.m[3][3] = a / b;
		return matrix;

	}
	mat4x4 Matrix_MakeTranslation(float x, float y, float z) {
		mat4x4 matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		matrix.m[3][0] = x;
		matrix.m[3][1] = y;
		matrix.m[3][2] = z;
		return matrix;

	}
	mat4x4 Matrix_MakeProjection(float fFovd, float fAspectRatio, float fNear, float fFar) {

		float fFovRad = 1.0f / tanf(fFovd * 0.5f / 180.0f * 3.14159f);
		mat4x4 matrix;
		matrix.m[0][0] = fAspectRatio * fFovRad;
		matrix.m[1][1] = fFovRad;
		matrix.m[2][2] = fFar / (fFar - fNear);
		matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
		matrix.m[2][3] = 1.0f;
		matrix.m[3][3] = 0.0f;
		return matrix;
	}
	mat4x4 Matrix_MultiplyMatrix(mat4x4& m1, mat4x4& m2) {
		mat4x4 matrix;
		for (int c = 0; c < 4; c++) {
			for (int r = 0; r < 4; r++) {
				matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
			}
		}
		return matrix;
	}
	mat4x4 Matrix_PointAt(vec3d& pos, vec3d& target, vec3d& up) {
		vec3d newForward = target - pos;
		newForward = newForward.normal();

		vec3d a = newForward * (up.dotprod(newForward));
		vec3d newUp = up - a;

		newUp = newUp.normal();

		vec3d newRight = newUp.cross(newForward);

		//point at MAtrix
		mat4x4 matrix;
		matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
		matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
		matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
		matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
		return matrix;
	}
	mat4x4 Matrix_QuickInverse(mat4x4& m) //only for roataiona nd translation matrix
	{
		mat4x4 matrix;
		matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
		matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
		matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
		matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
		matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
		matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
		matrix.m[3][3] = 1.0f;
		return matrix;
	}
	mat4x4 Matrix_PosAndRot(vec3d pos, float x, float y, float z) {
		mat4x4 matx = Matrix_MakeRotationX(x); mat4x4 maty = Matrix_MakeRotationY(y);
		mat4x4 matz = Matrix_MakeRotationZ(z); mat4x4 matT = Matrix_MakeTranslation(pos.x, pos.y, pos.z);
		mat4x4 matRelocate = Matrix_MultiplyMatrix(matz, matx);
		matRelocate = Matrix_MultiplyMatrix(matRelocate, maty);
		matRelocate = Matrix_MultiplyMatrix(matRelocate, matT);
		return matRelocate;
	}


	//triangle methods
	int Triangle_ClipAgainstPlane(vec3d plane_p, vec3d plane_n, triangle& in_tri, triangle& out_tri1, triangle& out_tri2)
	{//consider optimizing
		// Make sure plane normal is indeed normal
		plane_n.Normilize();

		// Return signed shortest distance from point to plane, plane normal must be normalised
		auto dist = [&](vec3d& p)
		{
			vec3d n = p.normal();
			return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - plane_n.dotprod(plane_p));
		};

		// Create two temporary storage arrays to classify points either side of plane
		// If distance sign is positive, point lies on "inside" of plane
		vec3d* inside_points[3];  int nInsidePointCount = 0;
		vec3d* outside_points[3]; int nOutsidePointCount = 0;
		vec2d* inside_tex[3]; int nInsideTexCount = 0;
		vec2d* outside_tex[3]; int nOutsideTexCount = 0;



		// Get signed distance of each point in triangle to plane
		float d0 = dist(in_tri.p[0]);
		float d1 = dist(in_tri.p[1]);
		float d2 = dist(in_tri.p[2]);

		if (d0 >= 0) { inside_points[nInsidePointCount] = &in_tri.p[0]; inside_tex[nInsidePointCount++] = &in_tri.t[0]; }
		else {
			outside_points[nOutsidePointCount] = &in_tri.p[0]; outside_tex[nOutsidePointCount++] = &in_tri.t[0];
		}
		if (d1 >= 0) {
			inside_points[nInsidePointCount] = &in_tri.p[1]; inside_tex[nInsidePointCount++] = &in_tri.t[1];
		}
		else {
			outside_points[nOutsidePointCount] = &in_tri.p[1]; outside_tex[nOutsidePointCount++] = &in_tri.t[1];
		}
		if (d2 >= 0) {
			inside_points[nInsidePointCount] = &in_tri.p[2]; inside_tex[nInsidePointCount++] = &in_tri.t[2];
		}
		else {
			outside_points[nOutsidePointCount] = &in_tri.p[2]; outside_tex[nOutsidePointCount++] = &in_tri.t[2];
		}
		
		// Now classify triangle points, and break the input triangle into 
		// smaller output triangles if required. There are four possible
		// outcomes...

		if (nInsidePointCount == 0)
		{
			// All points lie on the outside of plane, so clip whole triangle
			// It ceases to exist

			return 0; // No returned triangles are valid
		}

		if (nInsidePointCount == 3)
		{
			// All points lie on the inside of plane, so do nothing
			// and allow the triangle to simply pass through
			out_tri1 = in_tri;
			out_tri1.r = in_tri.r; out_tri1.g = in_tri.g; out_tri1.b = in_tri.b;

			return 1; // Just the one returned original triangle is valid
		}

		if (nInsidePointCount == 1 && nOutsidePointCount == 2)
		{
			// Triangle should be clipped. As two points lie outside
			// the plane, the triangle simply becomes a smaller triangle

			// Copy appearance info to new triangle
			out_tri1.r = in_tri.r; out_tri1.g = in_tri.g; out_tri1.b = in_tri.b;


			// The inside point is valid, so keep that...
			out_tri1.p[0] = *inside_points[0];
			out_tri1.t[0] = *inside_tex[0];

			// but the two new points are at the locations where the 
			// original sides of the triangle (lines) intersect with the plane
			float t;
			out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
			out_tri1.t[1].u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[1].v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1], t);
			out_tri1.t[2].u = t * (outside_tex[1]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[2].v = t * (outside_tex[1]->v - inside_tex[0]->v) + inside_tex[0]->v;

			return 1; // Return the newly formed single triangle
		}

		if (nInsidePointCount == 2 && nOutsidePointCount == 1)
		{
			
			// Triangle should be clipped. As two points lie inside the plane,
			// the clipped triangle becomes a "quad". Fortunately, we can
			// represent a quad with two new triangles

			// Copy appearance info to new triangles
			out_tri1.r = in_tri.r; out_tri1.g = in_tri.g; out_tri1.b = in_tri.b;


			out_tri2.r = in_tri.r; out_tri2.g = in_tri.g; out_tri2.b = in_tri.b;


			// The first triangle consists of the two inside points and a new
			// point determined by the location where one side of the triangle
			// intersects with the plane
			out_tri1.p[0] = *inside_points[0];
			out_tri1.p[1] = *inside_points[1];
			out_tri1.t[0] = *inside_tex[0];
			out_tri1.t[1] = *inside_tex[1];

			float t;
			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0], t);
			out_tri1.t[2].u = t * (outside_tex[0]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[2].v = t * (outside_tex[0]->v - inside_tex[0]->v) + inside_tex[0]->v;


			// The second triangle is composed of one of he inside points, a
			// new point determined by the intersection of the other side of the 
			// triangle and the plane, and the newly created point above
			out_tri2.p[0] = *inside_points[1];
			out_tri2.t[0] = *inside_tex[1];
			out_tri2.t[1] = out_tri1.t[2];
			out_tri2.p[1] = out_tri1.p[2];
			out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0], t);
			out_tri1.t[2].u = t * (outside_tex[0]->u - inside_tex[1]->u) + inside_tex[1]->u;
			out_tri1.t[2].v = t * (outside_tex[0]->v - inside_tex[1]->v) + inside_tex[1]->v;
		
			return 2; // Return two newly formed triangles which form a quad
		}
	}
	boolean Triangle_InPlane(vec3d plane_p, vec3d plane_n, triangle& in_tri)
	{
		// Make sure plane normal is indeed normal
		plane_n.Normilize();

		// Return signed shortest distance from point to plane, plane normal must be normalised
		auto dist = [&](vec3d& p)
		{
			
			return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z -plane_n.dotprod(plane_p));
		};

		// Create two temporary storage arrays to classify points either side of plane
		// If distance sign is positive, point lies on "inside" of plane
		vec3d* inside_points[3];  int nInsidePointCount = 0;
		vec3d* outside_points[3]; int nOutsidePointCount = 0;


		// Get signed distance of each point in triangle to plane
		float d0 = dist(in_tri.p[0]);
		float d1 = dist(in_tri.p[1]);
		float d2 = dist(in_tri.p[2]);

		if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; }
		else {
			outside_points[nOutsidePointCount++] = &in_tri.p[0];
		}
		if (d1 >= 0) {
			inside_points[nInsidePointCount++] = &in_tri.p[1];
		}
		else {
			outside_points[nOutsidePointCount++] = &in_tri.p[1];
		}
		if (d2 >= 0) {
			inside_points[nInsidePointCount++] = &in_tri.p[2];
		}
		else {
			outside_points[nOutsidePointCount++] = &in_tri.p[2];
		}

		// Now classify triangle points, and break the input triangle into 
		// smaller output triangles if required. There are four possible
		// outcomes...



		if (nInsidePointCount == 0)
		{


			return false; // Just the one returned original triangle is valid
		}

		return true;
	}
	int Triangle_InPlaneCritical(vec3d plane_p, vec3d plane_n, triangle& in_tri)
	{
		// Make sure plane normal is indeed normal
		plane_n.Normilize();

		// Return signed shortest distance from point to plane, plane normal must be normalised
		auto dist = [&](vec3d& p)
		{
			
			return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - plane_n.dotprod(plane_p));
		};

		// Create two temporary storage arrays to classify points either side of plane
		// If distance sign is positive, point lies on "inside" of plane
		vec3d* inside_points[3];  int nInsidePointCount = 0;
		vec3d* outside_points[3]; int nOutsidePointCount = 0;

		int p0 = 0;
		int p1 = 0;
		int p2 = 0; // if critical

		// Get signed distance of each point in triangle to plane
		float d0 = dist(in_tri.p[0]);
		float d1 = dist(in_tri.p[1]);
		float d2 = dist(in_tri.p[2]);

		if (d0 >= 0) {
			inside_points[nInsidePointCount++] = &in_tri.p[0];
			p0 = -1;
		}
		else {
			outside_points[nOutsidePointCount++] = &in_tri.p[0];
			p0 = 1;
		}
		if (d1 >= 0) {
			inside_points[nInsidePointCount++] = &in_tri.p[1];
			p1 = -2;
		}
		else {
			outside_points[nOutsidePointCount++] = &in_tri.p[1];
			p1 = 2;
		}
		if (d2 >= 0) {
			inside_points[nInsidePointCount++] = &in_tri.p[2];
			p2 = -3;
		}
		else {
			outside_points[nOutsidePointCount++] = &in_tri.p[2];
			p2 = 3;
		}

		// Now classify triangle points, and break the input triangle into 
		// smaller output triangles if required. There are four possible
		// outcomes...



		if (nInsidePointCount == 3)
		{
			return 0; // no critical
		}
		else if (nOutsidePointCount == 2) {
			return min(p0, min(p1, p2)); // -1 = p1 and p2 being critical -3 is p0 and p1 being critical
		}
		else if (nInsidePointCount == 2) {
			return max(p0, max(p1, p2)); // 1 = p0 critical 
		}
		return -999; // not in clip

	}

	//triangle list methods
	vector<triangle> transformTriangle(vector<triangle> vecTrianglesInWorld, mat4x4 matWorld, mat4x4 matView, mat4x4 matProj, bool inHand, vec3d c1) {
		vector<triangle> vecTri;



		sort(vecTrianglesInWorld.begin(), vecTrianglesInWorld.end(), [c1](triangle& t1, triangle& t2)
			{
				float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f - c1.z;
				float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f - c1.z;
				float x1 = (t1.p[0].x + t1.p[1].x + t1.p[2].x) / 3.0f - c1.x;
				float x2 = (t2.p[0].x + t2.p[1].x + t2.p[2].x) / 3.0f - c1.x;
				float y1 = (t1.p[0].y + t1.p[1].y + t1.p[2].y) / 3.0f - c1.y;
				float y2 = (t2.p[0].y + t2.p[1].y + t2.p[2].y) / 3.0f - c1.y;

				float a1 = z1 * z1 + x1 * x1 + y1 * y1;
				float a2 = z2 * z2 + x2 * x2 + y2 * y2;

				return abs(a1) > abs(a2);
			});

		for (auto& tri : vecTrianglesInWorld) {

			triangle triProjected, triTransformed, triViewed;

			triProjected.r = tri.r; triProjected.g = tri.g; triProjected.b = tri.b; triViewed.t[0] = tri.t[0]; triViewed.t[1] = tri.t[1]; triViewed.t[2] = tri.t[2];
			//cout << tri.r << ":" << tri.g <<":" << tri.b << endl;

			triTransformed = tri;
			triTransformed.p[0] = MatrixMultiplyVector(tri.p[0], matWorld);
			triTransformed.p[1] = MatrixMultiplyVector(tri.p[1], matWorld);
			triTransformed.p[2] = MatrixMultiplyVector(tri.p[2], matWorld);

			vec3d normal, line1, line2;

			line1 = triTransformed.p[1] - triTransformed.p[0];
			line2 = triTransformed.p[2] - triTransformed.p[0];

			// Take cross product of lines to get normal to triangle surface
			normal = line1.cross(line2);

			// You normally need to normalise a normal!
			normal.Normilize();

			vec3d vCameraRay = triTransformed.p[0] - player.vCamera;

			//if (normal.z < 0) 
			if (inHand || (normal.dotprod(vCameraRay) < 0.0f)) // TODO fix as in hand does not work properly (it seems to be working must investigate!)
			{
				vec3d light_direction = { 0,0,1 };
				light_direction = light_direction.normal();



				triViewed.p[0] = MatrixMultiplyVector(triTransformed.p[0], matView);
				triViewed.p[1] = MatrixMultiplyVector(triTransformed.p[1], matView);
				triViewed.p[2] = MatrixMultiplyVector(triTransformed.p[2], matView);

				int nClippedTriangles = 0;
				triangle clipped[2];
				nClippedTriangles = Triangle_ClipAgainstPlane({ 0.0f,0.0f,0.1f }, { 0.0f,0.0f,1.0f }, triViewed, clipped[0], clipped[1]);//loo int how this works exactly 

				//clipped[0] = triViewed;
				//nClippedTriangles = 1;

				for (int n = 0; n < nClippedTriangles; n++) {

					triProjected.t[0] = clipped[n].t[0]; triProjected.t[1] = clipped[n].t[1]; triProjected.t[2] = clipped[n].t[2];
					//3d to 2d
					triProjected.p[0] = MatrixMultiplyVector(clipped[n].p[0], matProj);
					triProjected.p[1] = MatrixMultiplyVector(clipped[n].p[1], matProj);
					triProjected.p[2] = MatrixMultiplyVector(clipped[n].p[2], matProj);

					triProjected.p[0] /= triProjected.p[0].w;
					triProjected.p[1] /= triProjected.p[1].w;
					triProjected.p[2] /= triProjected.p[2].w;

					vec3d vOffsetView = { 1, 1, 0 };

					triProjected.p[0] += vOffsetView;
					triProjected.p[1] += vOffsetView;
					triProjected.p[2] += vOffsetView;

					triProjected.p[0].x *= 0.5f * (float)ScreenWidth();
					triProjected.p[0].y *= 0.5f * (float)ScreenHeight();
					triProjected.p[1].x *= 0.5f * (float)ScreenWidth();
					triProjected.p[1].y *= 0.5f * (float)ScreenHeight();
					triProjected.p[2].x *= 0.5f * (float)ScreenWidth();
					triProjected.p[2].y *= 0.5f * (float)ScreenHeight();



					vecTri.push_back(triProjected);
				}

			}

		}
		return (vecTri);
	}
	vector<triangle> createaabbHitbox(entity e1, vector<triangle> vecTriin) {
		vector<triangle> vecTri;

		float minx = 100000; float miny = 100000; float minz = 100000; float maxx = -100000; float maxy = -100000; float maxz = -100000;

		entity obj = e1;  // load in all enitites 
		mat4x4 matx = Matrix_MakeRotationX(obj.Rotx); mat4x4 maty = Matrix_MakeRotationY(obj.Roty);
		mat4x4 matz = Matrix_MakeRotationZ(obj.Rotz); mat4x4 matT = Matrix_MakeTranslation(obj.position.x, obj.position.y, obj.position.z);
		mat4x4 matRelocate = Matrix_MultiplyMatrix(matz, matx);
		matRelocate = Matrix_MultiplyMatrix(matRelocate, maty);
		matRelocate = Matrix_MultiplyMatrix(matRelocate, matT);

		// loads in triangle of objects
		for (auto& tri : vecTriin) {
			for (int i = 0; i < 3; i++) {
				if (tri.p[i].x > maxx) maxx = tri.p[i].x;
				if (tri.p[i].y > maxy) maxy = tri.p[i].y;
				if (tri.p[i].z > maxz) maxz = tri.p[i].z;
				if (tri.p[i].x < minx) minx = tri.p[i].x;
				if (tri.p[i].y < miny) miny = tri.p[i].y;
				if (tri.p[i].z < minz) minz = tri.p[i].z;
			}
		}

		e1.aabbHitbox = aabb(minx, miny, minz, maxx, maxy, maxz);
		e1.hitboxType = 1; // aabb
		vec3d p000 = { minx,miny,minz };
		vec3d p001 = { minx,miny,maxz };
		vec3d p010 = { minx,maxy,minz };
		vec3d p011 = { minx,maxy,maxz };
		vec3d p100 = { maxx,miny,minz };
		vec3d p101 = { maxx,miny,maxz };
		vec3d p110 = { maxx,maxy,minz };
		vec3d p111 = { maxx,maxy,maxz };

		// return draable hitbox
		vecTri.push_back({ p000,p010,p110 });		vecTri.push_back({ p000, p110, p100 });
		vecTri.push_back({ p100,p110,p111 });	vecTri.push_back({ p100, p111, p101 });
		vecTri.push_back({ p101,p111,p011 });	vecTri.push_back({ p101, p011, p001 });
		vecTri.push_back({ p001,p011,p010 });	vecTri.push_back({ p001, p010, p000 });
		vecTri.push_back({ p010,p011,p111 });	vecTri.push_back({ p010, p111, p110 });
		vecTri.push_back({ p101,p001,p000 });	vecTri.push_back({ p101, p000, p100 });

		return vecTri;
	}
	vector<triangle> getMesh(meshTree top, mat4x4 secondTransform) {
		vector<triangle> output;
		mat4x4 transformation = Matrix_PosAndRot(top.myMeshNode.position, top.myMeshNode.Rotx, top.myMeshNode.Roty, top.myMeshNode.Rotz);
		mat4x4 critTransformation = Matrix_MakeTranslation(top.myMeshNode.position.x, top.myMeshNode.position.y, top.myMeshNode.position.z);
		transformation = Matrix_MultiplyMatrix(transformation, secondTransform);
		critTransformation = Matrix_MultiplyMatrix(critTransformation, secondTransform);

		if (top.hasSubs) {
			for (auto& sub : top.subMeshNodes) {
				vector<triangle> temp = getMesh(sub, transformation);
				for (auto& tri : temp) {
					output.push_back(tri);
				}
			}
		}

		for (auto& tri : top.myMeshNode.tri) {
			triangle newTri;

			for (int i = 0; i < 3; i++) {
				if (find(top.myMeshNode.criticalPoints.begin(), top.myMeshNode.criticalPoints.end(), tri.p[i]) != top.myMeshNode.criticalPoints.end()) {
					newTri.p[i] = MatrixMultiplyVector(tri.p[i], critTransformation);
				}
				else {
					newTri.p[i] = MatrixMultiplyVector(tri.p[i], transformation);
				}
			}

			newTri.r = tri.r; newTri.g = tri.g; newTri.b = tri.b;
			output.push_back(newTri);
		}

		return output;
	}

	//meshtreeMethods
	/*TODO:
	meshTree loadMeshTreeFromFile(file)
	meshTree create() // in struct
	meshTree saveToFile(meshTree)
	*/



	/*TODO figure out animation from file
	*/

	void Raster(vector<triangle> vecTrianglesToRaster, float r, float g, float b, int rastermode) {


		for (auto& triToRaster : vecTrianglesToRaster) {



		


			triangle clipped[2];
			list<triangle> listTriangles;
			listTriangles.push_back(triToRaster);
			int nNewTriangles = 1;


			for (int p = 0; p < 4; p++) {

				int nTrisToAdd = 0;
				while (nNewTriangles > 0) {
					triangle test = listTriangles.front();
					listTriangles.pop_front();
					nNewTriangles--;

					switch (p)
					{
					case 0:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 1:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, (float)ScreenHeight() + 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 2:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 3:	nTrisToAdd = Triangle_ClipAgainstPlane({ (float)ScreenWidth() + 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					}

					for (int w = 0; w < nTrisToAdd; w++)
						listTriangles.push_back(clipped[w]);
				}
				nNewTriangles = listTriangles.size();
			}
		
			switch (rastermode)
			{
			case 0:
				// Draw the transformed, viewed, clipped, projected, sorted, clipped triangles
				for (auto& t : listTriangles)
				{
					FillTriangle(t.p[0].x, t.p[0].y,
						t.p[1].x, t.p[1].y,
						t.p[2].x, t.p[2].y, olc::Pixel(r, g, b));

				}
				break;
			case 1: 	// Draw the transformed, viewed, clipped, projected, sorted, clipped triangles
				for (auto& t : listTriangles)
				{


					DrawTriangle(t.p[0].x, t.p[0].y,
						t.p[1].x, t.p[1].y,
						t.p[2].x, t.p[2].y, olc::Pixel(r, g, b));
				}
				break;
			case 2: 	// Draw the transformed, viewed, clipped, projected, sorted, clipped triangles
				for (auto& t : listTriangles)
				{
					FillTriangle(t.p[0].x, t.p[0].y,
						t.p[1].x, t.p[1].y,
						t.p[2].x, t.�a|� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\winver.h  �~����.P�t�OWsP-"�x#N�:�&���/� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\reason.h  0T�]�yǥ˘�0�H��I��G(������ C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\winreg.h  Zڤ�~����̓s���� ��*	���k��)C�X� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\wnnc.h  �:��2K] �
j�����e�
ʪ�3k�ӭ��W� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\winnetwk.h  ���r��Axa?E$��N�:���i�ٯa5.�P�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\cderr.h  �2�ܗ[�7���{:E��!2W����>�h��|� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\dde.h  ߧ�0�Ħ���~-�X���\P �ggD��� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\ddeml.h  �z��bsF׆E��j]���0�\b*xs*q�����A� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\dlgs.h  �="V�A�D��f�� ����7b�Ե�o1R� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\lzexpand.h  U��{�z��f���Hb�����"�6e`ax� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\mmsyscom.h  ��(`.�mQEo"���"����l��89f���ѕ C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\mciapi.h  z%1mQZ�J��f����H��O8�\+mp;���GO� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\mmiscapi.h  ��v%	��J�Bk��~̝��`��8g���{H� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\mmiscapi2.h  N&�H�J�?�����']:o�\ȑ�B�:;�d�X� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\playsoundapi.h  �pS~C�:�?+[[=vX.�Q�Bꁵ&���g� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\mmeapi.h  �V$9�
:��� �+;�&'��y*�=p�+�;k�R� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\timeapi.h  n�F�*T	�[�'�$9-��)Q���j�>�r��&� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\joystickapi.h  '��|��\��Zv���� ���s�����s�w� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\mmsystem.h  ��8s7�w�Ń&!��qA>����Y~a�k`N3�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\nb30.h  iu�	��&? ����#����u-��b@���F� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\rpcdcep.h  �.ו+�8W(M��6[���Z\��~���[�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\rpcdce.h  ��ˎ�jy*ƱT1A��`��#^$L�
�Ӭ �c�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\rpcnsi.h  ��f�u�Jn\��zx�\�G����q*�V�,��� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\rpcnterr.h  ��������Q�;�[$��ۅ=@����,���� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\rpcasync.h  ^����u���V�y������/�ta��6��ݢ C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\rpc.h  ;�]���)J��5م��s<`Γ|����� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\shellapi.h  �l{`G2`0m�p�{���rg؁��jo��� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\winperf.h  ,�p��)ъ����6e�y26�ޣ6�Xj+X	� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\inaddr.h  ��(�f��J�c��A��/�!c� ;b�r�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\winsock.h  ��]�	��[�yhF�|�a)�{Le6�����R7�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\bcrypt.h  l��;��j���I��}b��f��84�����Y@�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\ncrypt.h  v�E׵V�'<��$cS��T��m�i{}� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\dpapi.h  ���ڀ���1�x�)eF�;��5� ���5��i� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\wincrypt.h  �3j
y:/����`�'2Ős�T_l��ѰCg�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\winefs.h  �����zH"�HLUHFy]��I�uh�J��Cp �� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\rpcnsip.h  |�8�>�$��;�@�[VP���sVf*wۭg� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\rpcsal.h  �YY�c*��}+�;hs�(s��f�i@��n��ܱ��� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\rpcndr.h  ��.b��)���]Ý2��tW���Q�*��� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\wtypesbase.h  @&������?+ڽ粨5�[?D�F�d��GO�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\wtypes.h  UU�P���"@�~~ݤ�=�+�
p�LBu�Z5���� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\winioctl.h  J�R�ڴ���9��,Z):��6�&�7�p �� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\winsmcrd.h  ֒�O�	��]��Md��f G��f1eT�Q~�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\winscard.h  �Cţ/�	� G���n�c��J�.ɽo�T�E� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\prsht.h  Ȕ�J��XDW�v�h�W{�+��2$<F$��B� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\winspool.h  n�c��v����
	��������c`'�|��u� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\unknwnbase.h  ����{�9]�+����	&�r��(�Y���/� � C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\objidlbase.h  @e��sq\���&N)�H�ͬ�@�v�GK�[�Wl� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\cguid.h  �3P��"��[��91Y �\g׻��\�fC�g� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\combaseapi.h  1W�$~�C��%�mJ����~O0��50�yoii\� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\unknwn.h  ����d�0�t[�n\��.��/9��(�׈;u
i3� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\objidl.h  KY��`�	$����`�#)Q�B��+FF�)� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\oaidl.h  RO�b4ۓb��
���y�g>�^?����o���� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\propidlbase.h  x��n��X�^��Uw�B$�S��̜�e%2� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\coml2api.h  ?Z�.�A��q���i.cj�j�~A��3����B� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\oleidl.h  D�i��Ϻ&�Oͼ}3���a�F_�|���$�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\servprov.h  2��Zt��:H��mT":g$z�N�D�K�F� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\msxml.h  �MK*��!�d�i����������$���ղR�*
� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\urlmon.h  # b�����e=�D
f���;�ŒI����@�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\propidl.h  ��O�LC_��x3� ;�eTȅT .F�ǌ�f� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\objbase.h  �e���9LTU��V�R[&��0����C��P2� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\oleauto.h  �7��b�6C�M����u���D. +����9?e� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\ole2.h  ��dW�H�3�v�a�2�r�5����*l��UL� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\commdlg.h  ����V�Rʇ&kZ�gl˦u8�xQWۏz`k+�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\shared\stralign.h  p�& ;Tt<�ѝ0捶��K߿#56�f���� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\winsvc.h  %:���qI0fI�J��ͪ�����<� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\mcx.h  N���NE���\�H3�K4��9�|f�U	�Ӹ�r\� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\ime_cmodes.h  '�N�	8t����A��d�W�G$���l�0"$�&J� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\imm.h  �-��}n�5��L���:I����R��~(�y� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\windows.h  e��i�45���̯Zo�����q�+�<��� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\dpa_dsa.h  ��I��<��GF�pL*n��T=��G����jo_� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\commctrl.h  �[��v�7��S����%߅g^��٤.�٫-��(^� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\uxtheme.h  �����U���J��3K�h�s��4�M�6��)]� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\dwmapi.h  vu�� 1Ky᩶�cm�69�-0�N�M�o��r�B� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GL\gl.h  Z$��{Am�;h��f�8�.��zC�xh�S���� Z� %??_C@_0BD@IINOPBDD@wglSwapIntervalEXT@ � ��  � -   � C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusMem.h  �u�y�}Xfj�5_٘����j�?b!����u�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusBase.h  ����sr˳)N���R�J�i��:�	������� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusEnums.h  �ǿ׻��z��HQ�b��2��
�7-<	vf�}?� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusTypes.h  ��%����=z�Om�=����QC�
+1�s� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusInit.h  rgZ"=�pl޲�k� �[6����m����J� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusPixelFormats.h  ��/�a������ �+��!W���� �+��Ȉo� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusColor.h  ��ۯbj�~I4�r�C�U�����&�gOr� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusMetaHeader.h  ,������aN��AƖ<��Na��O�܅Ŗ� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusImaging.h  �)�Ϭ�Cd��]i��H��Il���0J1�t/u�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusColorMatrix.h  A!�ʶDp�-"x'n�D8�Kc�C�������� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusGpStubs.h  RA&�9V%�ܯ�EB�R򣔌u��O`ŕo���ה� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusHeaders.h  �UP>�9�}�o�5ׁ�8Vz�������2� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusFlat.h  (�����GӬվ�3�v���:!��j��$�M��� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusImageAttributes.h  �MFX��9�՟�w)���e6�\��M����+�V�&� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusMatrix.h  6��$EP��ܠ{y'3�1R��i+��w���ь C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusBrush.h  �d�X���W%�^Y��/�~���kl�4�I��� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusPen.h  �{~��n�fA�s(����U%Eꃽo��]�B� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusStringFormat.h  4����7�G�OViz�)�q��&/!��>9�^թm{�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusPath.h  �	U3���We��X��8���T�vݖ�)6�[� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusLineCaps.h  -��-'8��WFV<�R�^��ԏL���`��{�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusGraphics.h  ��H�����Z�S��x�j:H2-p�@:hE^�=7I� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusMetafile.h  �fح5d�4���Ӡ�,��^	A\RA�,� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusCachedBitmap.h  �)��Kxo�l?��f�!h]�<��d��Ou� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusRegion.h  ~�G���ϖ�9?�)>�.p'�ɉy�ƺs��H�~e�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusFontCollection.h  ��^!d�f���	��_b�߰��[��*�i�NJ� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusFontFamily.h  V!Q�-�}�o�-��@Nܝī�Q�➅�e�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusFont.h  ȥ7����qz���k<ҕ�����%��|g� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusBitmap.h  ��O�ҝ�
dj�~Mʝ���4OS�����Ԟ C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\GdiplusImageCodec.h  &���;�N�@U� ?\$�BJ:^���Q�h��{� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\gdiplus.h  B
��z�\"�J�ߪ�˙��'��&Η�K� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\shtypes.h  8�޲
yl��3�2���Νe0U�T��B1�� C:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\um\Shlwapi.h  �E8]���� 5	%�s��p�W�?U���&� � %??_C@_1CM@DMOOLNON@?$AAO?$AAL?$AAC?$AA_?$AAP?$AAI?$AAX?$AAE?$AAL?$AA_?$AAG?$AAA?$AAM?$AAE?$AA_@ � ,��  ��9    � %??_C@_11LOCGONAA@@ � ��  �:   �D:\2 - game devlopement\2.2 Projects\2.2 A - 3d game engine\Project1\olcPixelGameEngine.h  K3X;*%�wY�ᝒ�?l9�ԨH�n�H��rg� C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.31.31103\include\strstream  ht���$���_0]�C�O6$ԗ�>
��Ȅ�� 1� %??_C@_08BFENOALC@cube?4obj@ � 	��  �|    T� %??_C@_0L@PJFINJCL@ground?4obj@ � ��  ��    ޶ %??_C@_0BB@MBADCGAF@draw?5weapon?5one?$CB@ � ��  �S    � %??_C@_0M@NLLNAAEH@weapon?5two?$CB@ � ��  ��=    � %??_C@_0O@DDNDCEA@weapon?5three?$CB@ � ��  �S    � %??_C@_0L@FOAHNDLL@empty?5hand@ � ��  ��    �� %??_C@_0M@OMALDLBH@weapon?5one?$CB@ � ��  ��=    �� %??_C@_07ONCMKLLM@FPSGame@ � ��  �o    e� %??_C@_09BIJCCNPP@no?5raster@ � 
��  �.)    s� %??_C@_0O@ENOFGLJI@akgunauto?4obj@ � ��  �S    t� %??_C@_0M@NBILOAIF@testArm?4obj@ � ��  ��=     D:\2 - game devlopement\2.2 Projects\2.2 A - 3d game engine\Project1\Dynamic3DEngine.cpp  V�O��q�ӫ(~C!�CUR�?.D�iJj�\e	 A� %??_C@_0CC@PGBJLEPD@invalid?5wchar_t?5filename?5argume@ � "��  �Z!    `� %??_C@_0P@MNLEEHMK@bad?5conversion@ � ��  �h    �� %??_C@_19OKANBMPJ@?$AA?$CB?$AA?$CF?$AAx?$AA?$AA@ � 
��  ��E    .� %??_C@_0BG@ELKJEKEA@?$CFb?5?$CFd?5?$CFH?5?3?5?$CFM?5?3?5?$CFS?5?$CFY@ � ��  ��    4� %??_C@_0N@LIAHCJPJ@?$CFm?5?1?5?$CFd?5?1?5?$CFy@ � ��  �    >� %??_C@_01CLKCMJKC@?5@ � ��  �p    @� %??_C@_0N@LPFKKEBD@?3AM?3am?3PM?3pm@ � ��  �    E� %??_C@_0BA@HGEPHOKI@?$CFI?5?3?5?$CFM?5?3?5?$CFS?5?$CFp@ � ��  �}    G� %??_C@_07GLPBMING@?$CFH?5?3?5?$CFM@ � ��  �o    K� %??_C@_0N@GEFMGHCI@?$CFH?5?3?5?$CFM?5?3?5?$CFS@ � ��  �    P� %??_C@_0N@LADGFGLL@?$CFd?5?1?5?$CFm?5?1?5?$CFy@ � ��  �    � %??_C@_04IHCGILC@?$CB?$CFx?$AA@ � ��  ��    � %??_C@_05NNKBBLJI@?$CF?40Lf@ � ��  ��    7� %??_C@_02BBAHNLBA@?$CFp@ � ��  �    j� %??_C@_02CLHGNPPK@Lu@ � ��  �    s� %??_C@_02HIKPPMOK@Ld@ � ��  �    |� %??_C@_02BDDLJJBK@lu@ � ��  �    �� %??_C@_02EAOCLKAK@ld@ � ��  �    �� %??_C@_0BJ@DHFDPMIM@invalid?5vector?5subscript@ � ��  ��J    t� %??_C@_02MDKMJEGG@eE@ � ��  �    u� %??_C@_02OOPEBDOJ@pP@ � ��  �    x� %??_C@_01LFCBOECM@?4@ � ��  �p    �� %??_C@_0BI@CFPLBAOH@invalid?5string?5position@ � ��  ��    	� %??_C@_1BK@MHIKGOKE@?$AA?3?$AAA?$AAM?$AA?3?$AAa?$AAm?$AA?3?$AAP?$AAM?$AA?3?$AAp?$AAm@ � ��  ��O    �� %??_C@_0O@NKNMEGII@list?5too?5long@ � ��  �S    W� %??_C@_0BA@FOIKENOD@vector?5too?5long@ � ��  �}    �� %??_C@_04FJCIEOIB@$?$CLxv@ � ��  ��    � %??_C@_0HJ@PJOMAIHM@?$CLv$x?$CLv$xv$?$CLxv?$CL$xv$?$CLx?$CL$vx?$CL$vx$v?$CL@ � y��  �mS    � %??_C@_13IMODFHAA@?$AA?9@ � ��  �oS    �� %??_C@_01JOAMLHOP@?9@ � ��  �p   	� .data$rs      	� .rdata$r CONST     	�� .text$yd  CODE     	�� .text$di  CODE     	�� .CRT$XCA DATA     	�� .CRT$XCU DATA     	�� .CRT$XCL DATA     	�� .CRT$XCC DATA     	�� .CRT$XCZ DATA     	�� .CRT$XDU DATA     	�� .CRT$XDL DATA     	�� .CRT$XDC DATA     	�� .CRTMA$XCA DATA     	�� .CRTMA$XCU DATA     	�� .CRTMA$XCL DATA     	�� .CRTMA$XCC DATA     	�� .CRTMA$XCZ DATA     	�� .CRTMP$XCA DATA     	�� .CRTMP$XCU DATA     	�� .CRTMP$XCL DATA     	�� .CRTMP$XCC DATA     	�� .CRTMP$XCZ DATA     	�� .CRTVT$XCA DATA     	� .CRTVT$XCU DATA     	~� .CRTVT$XCL DATA     	}� .CRTVT$XCC DATA     	|� .CRTVT$XCZ DATA     	�4.data$r      	�4.xdata$x CONST     � 	nodiscard � ??2@YAPAXIPAX@Z � ��  	 �lP  �&  ��    ȝ` ��� acosf �  �  �&  �+  �U  ��`  � asinf �  �  �&  ��  �t  ��`  � atan2f �  �  �&  �O	  ��  ��`  � atanf �  �  �&  ��	  ��  ��`  � cosf �  �  �&  �  �  ��`  � coshf �  �  �&  ��  �!  ��`  � expf �  �  �&  �A  �@  ��`  � fabsf �  �  �&  ��  �_  ��`  � fmodf �  �  �&  ��  ��  ��`  � frexpf �  �  �&  ��  ��  Ș` �� hypotf �  �  �&  �:  �  ��`  � ldexpf �  �  �&  ��  �4  Ș` � log10f �  �  �&  �r  �e  ��`  � logf �  �  �&  �  ��  ��`  � modff �  �#  �&  ��  ��  Ș` � powf �  �  �&  �k  ��  ��`  � sinf �  �  �&  �  �/  ��`  � sinhf �  �  �&  ��  �N  ��`    sqrtf �  �  �&  �5  �m  ��`   tanf �  �  �&  ��  ��  ��`   tanhf �  �  �&  �U  ��  ��`   acosl �  �  �&  ��  ��  ��`   asinl �  �  �&  �}  ��  ��`   atan2l �  �  �&  �  �  ��`    atanl �  �  �&  ��  �9  ��`  7 coshl �  �  �&  �  ��  ��`  ; cosl �  �  �&  ��  ��  ��`  C expl �  �  �&  �9  �  ��`  K fabsl �  �  �&  ��  �,  ��`  a fmodl �  �  �&  ��  �j  ��`  g frexpl �  �  �&  ��  ��  Ș` �t hypotl �  �  �&  ��  �  ��`  z ldexpl �  �  �&  �~  �2  Ș` � logl �  �  �&  �   �c  ��`  � log10l �  �  �&  ��   ��  ��`  � modfl �  �  �&  �?!  ��  Ș` � powl �  �  �&  �"  ��  ��`  � sinhl �  �  �&  ��"  �-	  ��`  � sinl �  �  �&  �R#  �L	  ��`  � sqrtl �  �  �&  ��#  �k	  ��`  � tanhl �  �  �&  �r$  ��	  ��`  � tanl �  �  �&  �%  ��	  ��`  z 	nodiscard } ?abs@@YAMM@Z � ��  �  �&  ��(  ��
  ȝ` ��z 	nodiscard  ?acos@@YAMM@Z � ��  �  �&  ��)  ��
  ȝ` ��# 	nodiscard ' ?atan2@@YAMMM@Z � ��  �  �&  �y-  ��  ȝ` ��#9 	nodiscard < ?cos@@YAMM@Z � ��  �  �&  ��/  �a  ȝ` ��9g 	nodiscard j ?floor@@YAMM@Z � ��  �  �&  �K5  ��  ȝ` �g 	nodiscard  ?sin@@YAMM@Z � ��  �  �&  �@I  ��  ȝ` �� 	nodiscard  ?sqrt@@YAMM@Z � ��  �  �&  ��J  ��  ȝ` � ' __local_stdio_printf_options �  �C  �&  �Fs  �  ��` *, __local_stdio_scanf_options �  �C  �&  ��s  �<  ��` /� _vfwprintf_l �  �E  �&  �
t  �g  Ș` '�� _vfwprintf_s_l �  �E  �&  ��u  �@  Ș` '�� _vfwprintf_p_l �  �E  �&  �8w  �  Ș` '�� _vfwscanf_l �  �E  �&  ���  �{&  Ș` ,�� _vfwscanf_s_l �  �E  �&  �4�  �T'  Ș` ,�� _vsnwprintf_l �  �[  �&  ��  �3-  -Ș` '�� _vsnwprintf_s_l �  �Y  �&  �=�  ��-  /Ș` '�� _vswprintf_c_l �  �[  �&  �*�  ��/  +Ș` '�� _vswprintf_l �  �[  �&  � �  �1  Ș` �� __vswprintf_l �  �_  �&  ��  ��1  Ș` �� vswprintf �  �d  �&  �t�  �x2  Ș` � _vswprintf_s_l �  �[  �&  �8�  ��2  +Ș` '� _vswprintf_p_l �  �[  �&  �.�  �+4  +Ș` '�# _vscwprintf_l �  �]  �&  �$�  �_5  )Ș` '�- _vscwprintf_p_l �  �]  �&  ��  �6  )Ș` '�� swprintf_s �` �  �&  ���  ��9  -��` � _scwprintf �` �  �&  ���  �V@  '��` #� _vswscanf_l �  �h  �&  �?�  �lC  Ș` ,�� _vswscanf_s_l �  �h  �&  ���  �ED  Ș` ,� _vsnwscanf_l �  �m  �&  ���  �E  Ș` ,� _vsnwscanf_s_l �  �m  �&  ���  ��E   Ș` ,�� _vfprintf_l �  �t  �&  ���  ��J  Ș` '�� _vfprintf_s_l �  �t  �&  ��  ��K  Ș` '�� _vfprintf_p_l �  �t  �&  ���  ��L  Ș` '�� _vfscanf_l �  �t  �&  ��  �U  Ș` ,�� _vfscanf_s_l �  �t  �&  ���  ��U  Ș` ,�� _vsnprintf_l �  �|  �&  �~  ��[  -Ș` 'w� _vsnprintf �  ��  �&  �� �\  Ș` �� vsnprintf �  ��  �&  �z ��\  +Ș` 'w� _vsprintf_l �  ��  �&  �� ��]  Ș` �� _vsprintf_s_l �  �|  �&  �4 �q^  +Ș` '~� _vsprintf_p_l �  �|  �&  �* ��_  +Ș` '�� _vsnprintf_s_l �  ��  �&  � 	 ��`  /Ș` '�� _vscprintf_l �  ��  �&  � ��b  )Ș` 'w _vscprintf_p_l �  ��  �&  �� ��c  )Ș` '�	 _vscprintf_p �  ��  �&  �� �d  Ș`  _vsnprintf_c_l �  �|  �&  �� �_d  +Ș` 'wV sprintf_s �` �  �&  �K �Qg  -��` �X _scprintf �` �  �&  �5# �_n  '��` �� _vsscanf_l �  ��  �&  ��& ��o  Ș` ,�� _vsscanf_s_l �  ��  �&  �D( �jp  Ș` ,�� vsscanf_s �  ��  �&  �9) ��p  Ș` �K _vcwprintf_l �  �]  �&  ��B ��{  Ș` '=T _vcwprintf_s_l �  �]  �&  ��C �A|  Ș` 'B] _vcwprintf_p_l �  �]  �&  �QE ��|  Ș` 'G� _vcwscanf_l �  �]  �&  ��M ��  Ș` ,�� _vcwscanf_s_l �  �]  �&  �'O ���  Ș` ,��" wmemchr �  ��  �&  ��Z ���  ��`  �" wmemcmp �  �,S  �&  �\ ���  (��`  �" wmemcpy �  ��  �&  �`] �p�  ��`  �" wmemset �  ��N  �&  ��^ �
�  ��`  �% 	nodiscard �% ?min@?$numeric_limits@_N@std@@SA_NXZ � ��   ��  �'  ��` ��  
ȝ`  ��%& 	nodiscard 	& ?min@?$numeric_limits@D@std@@SADXZ � ��  ��  �'  ��e �W�  
ȝ`  �&6& 	nodiscard 7& ?min@?$numeric_limits@C@std@@SACXZ � ��  ��  �'  ��j �̈ 