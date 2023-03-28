/*
Sound : TODO 
Lighting : TODO
Animation - keyframes : TODO 
Animations - states : TODO
Animations - State Transition Function : TODO 
Animation - grabbing : WBN
Physics - hitboxes : TODO 
Physics - kinamatics : TODO


*/


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
	float w = 1;
	vec2d operator/=(const float rhs) {
		this->u /= rhs;
		this->v /= rhs;
		return *this;
	}
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
	vector<vec2d> texs;

	float r = 0;
	float g = 0;
	float b = 0;



	bool loadFromObjectFile(string sFilename, bool hasTexture = false)
	{
		ifstream f(sFilename);
		if (!f.is_open())
			return false;

		int count = 0;

		while (!f.eof()) {
			char line[128];
			f.getline(line, 128);

			strstream s;
			s << line;
			char junk;

			if (line[0] == 'v') {
				if (line[1] == 't') {
					vec2d v;
					s >> junk >> junk >> v.u >> v.v;
					//std::cout << "text" << v.u << v.v << endl;
					texs.push_back(v);
				}
				else {
					vec3d v;
					s >> junk >> v.x >> v.y >> v.z;
					verts.push_back(v);
				}
			}
			if (!hasTexture) {
				if (line[0] == 'f') {
					int f[3];
					s >> junk >> f[0] >> f[1] >> f[2];
					tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
				}
			}
			else {
				if (line[0] == 'f' ) {
					
					s >> junk;

					string tokens[6];
					int nTokenCount = -1;
					
					while (!s.eof()) {
						char c = s.get();
						if (c == ' ' || c == '/')
							nTokenCount++;
						else
							tokens[nTokenCount].append(1, c);
					
					}
					tokens[nTokenCount].pop_back();
					
					triangle tri; 
					tri.p[0] = verts[stoi(tokens[0]) - 1];
					tri.p[1] = verts[stoi(tokens[2]) - 1];
					tri.p[2] = verts[stoi(tokens[4]) - 1];

					tri.t[0] = texs[stoi(tokens[1]) - 1];
					tri.t[1] = texs[stoi(tokens[3]) - 1];
					tri.t[2] = texs[stoi(tokens[5]) - 1];
					//std::cout << tri.p[0].x << ":" << tri.t[0].u << "end";

					/*tris.push_back({ verts[stoi(tokens[0]) - 1], verts[stoi(tokens[2]) - 1], verts[stoi(tokens[4]) - 1],
						texs[stoi(tokens[1]) - 1], texs[stoi(tokens[3]) - 1], texs[stoi(tokens[5]) - 1] });
					*/
					
					tris.push_back(tri);
					count++;
					
				}
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
	float Scale = 1;

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
		int x, y, z, count;
		x = y = z = count = 0;
		for (auto& obj : cP) {
			criticalPoints.push_back(obj);
			x = obj.x; y = obj.y; x = obj.z; count++;
		}
		toOrigin = {(float) ( - 1.0 * x / count),(float)(-1.0 * y / count), (float)(-1.0 * z / count)};
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
	void clearModel() {
		myMeshNode.tri.clear();
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
	entity(float a, float b, float c, string filename, bool hasTexture = false) {
		position.x = a; position.y = b; position.z = c; model.loadFromObjectFile(filename, hasTexture);
	}

	entity(float a, float b, float c, string filename, float rx, float ry, float rz, bool hasTexture = false) {
		position.x = a; position.y = b; position.z = c; model.loadFromObjectFile(filename, hasTexture);
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
	meshTree topTree;
	meshTree *currentNode;
	int lastSelected = 0;
	float timeChecker;
	player player;
	std::unique_ptr<olc::Sprite> tempSprite;
	float* m_DepthBuffer = nullptr;
	entity e1;
	double pi = 2 * acos(0.0);

	vec3d Vector_IntersectPlane(vec3d& plane_p, vec3d& plane_n, vec3d& lineStart, vec3d& lineEnd, float& t)
	{
		plane_n = plane_n.normal();
		float plane_d = -plane_n.dotprod(plane_p);
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
	mat4x4 Matrix_MakeTranslation(vec3d v) {
		mat4x4 matrix;
		matrix.m[0][0] = 1.0f;
		matrix.m[1][1] = 1.0f;
		matrix.m[2][2] = 1.0f;
		matrix.m[3][3] = 1.0f;
		matrix.m[3][0] = v.x;
		matrix.m[3][1] = v.y;
		matrix.m[3][2] = v.z;
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
	mat4x4 Matrix_MakeRelocate(vec3d pos, float x, float y, float z) {
		mat4x4 matx = Matrix_MakeRotationX(x); mat4x4 maty = Matrix_MakeRotationY(y);
		mat4x4 matz = Matrix_MakeRotationZ(z); mat4x4 matT = Matrix_MakeTranslation(pos.x, pos.y, pos.z);
		mat4x4 matRelocate = Matrix_MultiplyMatrix(matz, matx);
		matRelocate = Matrix_MultiplyMatrix(matRelocate, maty);
		matRelocate = Matrix_MultiplyMatrix(matRelocate, matT);
		return matRelocate;
	}
	mat4x4 Matrix_MakeRelocate(mat4x4& matRx, mat4x4& matRy, mat4x4& matRz, mat4x4& matT) {
		mat4x4 matRelocate = Matrix_MultiplyMatrix(matRz, matRx);
		matRelocate = Matrix_MultiplyMatrix(matRelocate, matRy);
		matRelocate = Matrix_MultiplyMatrix(matRelocate, matT);

		return matRelocate;
	}
	mat4x4 Matrix_Relocate(mat4x4& matRelocate, vec3d pos, float x, float y, float z) {
		mat4x4 matx = Matrix_MakeRotationX(x); mat4x4 maty = Matrix_MakeRotationY(y);
		mat4x4 matz = Matrix_MakeRotationZ(z); mat4x4 matT = Matrix_MakeTranslation(pos.x, pos.y, pos.z);
		 matRelocate = Matrix_MultiplyMatrix(matRelocate, matz);
		 matRelocate = Matrix_MultiplyMatrix(matRelocate, matx);
		matRelocate = Matrix_MultiplyMatrix(matRelocate, maty);
		matRelocate = Matrix_MultiplyMatrix(matRelocate, matT);
		return matRelocate;
	}
	mat4x4 Matrix_Relocate(mat4x4& matRelocate, mat4x4& matRx, mat4x4& matRy, mat4x4& matRz, mat4x4& matT) {
		 matRelocate = Matrix_MultiplyMatrix(matRelocate,matRz);
		 matRelocate = Matrix_MultiplyMatrix(matRelocate, matRx);
		matRelocate = Matrix_MultiplyMatrix(matRelocate, matRy);
		matRelocate = Matrix_MultiplyMatrix(matRelocate, matT);

		return matRelocate;
	}

	//triangle methods
	int Triangle_ClipAgainstPlane(vec3d plane_p, vec3d plane_n, triangle& in_tri, triangle& out_tri1, triangle& out_tri2)
	{//consider optimizing
		// Make sure plane normal is indeed normal
		plane_n.Normilize();

		out_tri1.t[0] = in_tri.t[0];
		out_tri2.t[0] = in_tri.t[0];
		out_tri1.t[1] = in_tri.t[1];
		out_tri2.t[1] = in_tri.t[1];
		out_tri1.t[2] = in_tri.t[2];
		out_tri2.t[2] = in_tri.t[2];

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
			out_tri1.r = in_tri.r; out_tri1.g = in_tri.g; out_tri1.b = in_tri.b; out_tri1.t[0] = in_tri.t[0]; out_tri1.t[1] = in_tri.t[1]; out_tri1.t[2] = in_tri.t[2];

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
			out_tri1.t[1].w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;

			out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1], t);
			out_tri1.t[2].u = t * (outside_tex[1]->u - inside_tex[0]->u) + inside_tex[0]->u;
			out_tri1.t[2].v = t * (outside_tex[1]->v - inside_tex[0]->v) + inside_tex[0]->v;
			out_tri1.t[2].w = t * (outside_tex[1]->w - inside_tex[0]->w) + inside_tex[0]->w;


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
			out_tri1.t[2].w = t * (outside_tex[0]->w - inside_tex[0]->w) + inside_tex[0]->w;

			// The second triangle is composed of one of he inside points, a
			// new point determined by the intersection of the other side of the 
			// triangle and the plane, and the newly created point above
			out_tri2.p[0] = *inside_points[1];
			out_tri2.t[0] = *inside_tex[1];
			out_tri2.p[1] = out_tri1.p[2];
			out_tri2.t[1] = out_tri1.t[2];
			out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0], t);
			out_tri2.t[2].u = t * (outside_tex[0]->u - inside_tex[1]->u) + inside_tex[1]->u;
			out_tri2.t[2].v = t * (outside_tex[0]->v - inside_tex[1]->v) + inside_tex[1]->v;
			out_tri2.t[2].w = t * (outside_tex[0]->w - inside_tex[1]->w) + inside_tex[1]->w;
			
			return 2;
		}
	}
	boolean Triangle_InPlane(vec3d plane_p, vec3d plane_n, triangle& in_tri)
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
	triangle transformTriangleIntoWorld(float scale, mat4x4& relocateMatrix, triangle& tri) {
		triangle scaledTri, newTri;
		newTri.t[0] = tri.t[0];
		newTri.t[1] = tri.t[1];
		newTri.t[2] = tri.t[2];
		scaledTri.p[0] = tri.p[0] * scale;
		scaledTri.p[1] = tri.p[1] * scale;
		scaledTri.p[2] = tri.p[2] * scale;


		newTri.p[0] = MatrixMultiplyVector(scaledTri.p[0], relocateMatrix);
		newTri.p[1] = MatrixMultiplyVector(scaledTri.p[1], relocateMatrix);
		newTri.p[2] = MatrixMultiplyVector(scaledTri.p[2], relocateMatrix);
		
		return newTri;
	}
	triangle transformTriangleIntoWorldRespectCrits(float scale, mat4x4& relocateMatrix, triangle& tri) {
		triangle scaledTri, newTri;
		scaledTri.p[0] = tri.p[0] * scale;
		scaledTri.p[1] = tri.p[1] * scale;
		scaledTri.p[2] = tri.p[2] * scale;
		switch (tri.criticalPoints) {

		case 1:
			newTri.p[0] = scaledTri.p[0];
			newTri.p[1] = MatrixMultiplyVector(scaledTri.p[1], relocateMatrix);
			newTri.p[2] = MatrixMultiplyVector(scaledTri.p[2], relocateMatrix);
			break;
		case 2:
			newTri.p[0] = MatrixMultiplyVector(scaledTri.p[0], relocateMatrix);
			newTri.p[1] = scaledTri.p[1];
			newTri.p[2] = MatrixMultiplyVector(scaledTri.p[2], relocateMatrix);
			break;
		case 3:
			newTri.p[0] = MatrixMultiplyVector(scaledTri.p[0], relocateMatrix);
			newTri.p[1] = MatrixMultiplyVector(scaledTri.p[1], relocateMatrix);
			newTri.p[2] = scaledTri.p[2];
			break;
		case -1:
			newTri.p[0] = MatrixMultiplyVector(scaledTri.p[0], relocateMatrix);
			newTri.p[1] = scaledTri.p[1];
			newTri.p[2] = scaledTri.p[2];
			break;
		case -2:
			newTri.p[0] = scaledTri.p[0];
			newTri.p[1] = MatrixMultiplyVector(scaledTri.p[1], relocateMatrix);
			newTri.p[2] = scaledTri.p[2];
			break;
		case -3:
			newTri.p[0] = scaledTri.p[0];
			newTri.p[1] = scaledTri.p[1];
			newTri.p[2] = MatrixMultiplyVector(scaledTri.p[2], relocateMatrix);
			break;


		default:
			newTri.p[0] = MatrixMultiplyVector(scaledTri.p[0], relocateMatrix);
			newTri.p[1] = MatrixMultiplyVector(scaledTri.p[1], relocateMatrix);
			newTri.p[2] = MatrixMultiplyVector(scaledTri.p[2], relocateMatrix);
			break;
		}

		newTri.t[0] = tri.t[0];
		newTri.t[1] = tri.t[1];
		newTri.t[2] = tri.t[2];
		return newTri;
	}
	
	//triangle list methods
	vector<triangle> triangleRelocate(vector<triangle> model, mat4x4 matRelocate, int scale) {
		vector<triangle> triToWorld;
		for (auto& tri : model) {
			triangle newTri = transformTriangleIntoWorld(scale, matRelocate, tri);
			newTri.r = 155; newTri.g = 40; newTri.b = 40;
			newTri.t[0] = tri.t[0]; newTri.t[1] = tri.t[1]; newTri.t[2] = tri.t[2];
			triToWorld.push_back(newTri);// send it to world 	
		}
		return triToWorld;
	}
	vector<triangle> triangleRelocate(vector<triangle> model, vec3d pos, float x, float y, float z, int scale) {
		vector<triangle> triToWorld;
		mat4x4 matRelocate = Matrix_MakeRelocate(pos, x, y, z);
		for (auto& tri : model) {
			triangle newTri = transformTriangleIntoWorld(scale, matRelocate, tri);
			newTri.r = 155; newTri.g = 40; newTri.b = 40;
			newTri.t[0] = tri.t[0]; newTri.t[1] = tri.t[1]; newTri.t[2] = tri.t[2];
			triToWorld.push_back(newTri);// send it to world 	
		}
		return triToWorld;
	}
	void triangleRelocate(vector<triangle> model, mat4x4 matRelocate, int scale, vector<triangle>& triList) {
		for (auto& tri : model) {
			triangle newTri = transformTriangleIntoWorld(scale, matRelocate, tri);
			newTri.r = 155; newTri.g = 40; newTri.b = 40;
			newTri.t[0] = tri.t[0]; newTri.t[1] = tri.t[1]; newTri.t[2] = tri.t[2];
			triList.push_back(newTri);// send it to world 	
		}
		
	}
	vector<triangle> transformTriangleIntoProjection(vector<triangle> vecTrianglesInWorld, mat4x4 matWorld, mat4x4 matView, mat4x4 matProj, bool inHand, vec3d c1)
	{
		vector<triangle> vecTri;



		/*sort(vecTrianglesInWorld.begin(), vecTrianglesInWorld.end(), [c1](triangle& t1, triangle& t2)
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
			});*/

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
					triProjected.t[0] /= triProjected.p[0].w;
					triProjected.t[1] /= triProjected.p[1].w;
					triProjected.t[2] /= triProjected.p[2].w;

					triProjected.t[0].w = 1.0f / triProjected.p[0].w;
					triProjected.t[1].w = 1.0f / triProjected.p[1].w;
					triProjected.t[2].w = 1.0f / triProjected.p[2].w;

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
	vector<triangle> getMesh(meshTree top, mat4x4 secondTransform) { // will probably used remade version below 
		vector<triangle> output;
		mat4x4 transformation = Matrix_MakeRelocate(top.myMeshNode.position, top.myMeshNode.Rotx, top.myMeshNode.Roty, top.myMeshNode.Rotz);
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
	
	void triangleRelocateKeepCritPoints(vector<triangle> model, mat4x4 matRelocate, int scale, vector<triangle>& triList) {
		for (auto& tri : model) {
			triangle newTri = transformTriangleIntoWorldRespectCrits(scale, matRelocate, tri);
			newTri.r = 155; newTri.g = 40; newTri.b = 40;
			triList.push_back(newTri);// send it to world 	
		}
	}

	//meshtreeMethods
	/*TODO:
	meshTree loadMeshTreeFromFile(file)
	meshTree create() // in struct
	meshTree saveToFile(meshTree)
	*/



	/*TODO figure out animation from file
	*/


	
	vector<triangle> MeshTreeToWorld(meshTree TopBranch) {
		mat4x4 matRelocate = Matrix_MakeRelocate({ TopBranch.myMeshNode.position.x, TopBranch.myMeshNode.position.y, TopBranch.myMeshNode.position.z }, TopBranch.myMeshNode.Rotx, TopBranch.myMeshNode.Roty, TopBranch.myMeshNode.Rotz);
		vector<triangle> triList = triangleRelocate(TopBranch.myMeshNode.tri,matRelocate,TopBranch.myMeshNode.Scale);
		

		for (auto& leaf : TopBranch.subMeshNodes) {
			MeshTreeSubtoWorld(leaf, matRelocate, triList);
		}
		
		return triList;
	
	}

	void MeshTreeSubtoWorld(meshTree Branch, mat4x4 matRelocateSuper, vector<triangle>& triList) {
		vector<triangle> branchModel; 
		mat4x4 matRelocate;
		
		//move to rotation point 
		matRelocate = Matrix_MakeTranslation(Branch.myMeshNode.toOrigin);

		//rotate and offset Position position
		Matrix_Relocate(matRelocate, Branch.myMeshNode.position, Branch.myMeshNode.Rotx, Branch.myMeshNode.Roty, Branch.myMeshNode.Rotz);
		triangleRelocateKeepCritPoints(Branch.myMeshNode.tri, matRelocate, Branch.myMeshNode.Scale, branchModel);
		//undo the move to rotation point 
		Matrix_Relocate(matRelocate, Branch.myMeshNode.toOrigin*-1, 0,0,0);
		//apply super Relocation Matrix 
		matRelocate = Matrix_MultiplyMatrix(matRelocate, matRelocateSuper);
		branchModel = triangleRelocate(branchModel, Branch.myMeshNode.toOrigin * -1,0,0,0,1);

		triangleRelocateKeepCritPoints(branchModel, matRelocateSuper, 1, triList);//apply transformation to triangles and add to triList


		//do subs 
		for (auto& leaf : Branch.subMeshNodes) {
			MeshTreeSubtoWorld(leaf, matRelocate, triList);
		}
	}

	
	void DrawTexture(int x1, int y1, float u1, float v1, float w1,
		int x2, int y2, float u2, float v2, float w2,
		int x3, int y3, float u3, float v3, float w3, olc::Sprite* spr)

	{
			if (y2 < y1)
			{
				swap(y1, y2);
				std::swap(x1, x2);
				std::swap(u1, u2);
				std::swap(v1, v2);
				std::swap(w1, w2);
			}

			if (y3 < y1)
			{
				std::swap(y1, y3);
				std::swap(x1, x3);
				std::swap(u1, u3);
				std::swap(v1, v3);
				std::swap(w1, w3);
			}

			if (y3 < y2)
			{
				std::swap(y2, y3);
				std::swap(x2, x3);
				std::swap(u2, u3);
				std::swap(v2, v3);
				std::swap(w2, w3);
			}

			int dy1 = y2 - y1;
			int dx1 = x2 - x1;
			float dv1 = v2 - v1;
			float du1 = u2 - u1;
			float dw1 = w2 - w1;

			int dy2 = y3 - y1;
			int dx2 = x3 - x1;
			float dv2 = v3 - v1;
			float du2 = u3 - u1;
			float dw2 = w3 - w1;

			float tex_u, tex_v, tex_w;

			float dax_step = 0, dbx_step = 0,
				du1_step = 0, dv1_step = 0,
				du2_step = 0, dv2_step = 0,
				dw1_step = 0, dw2_step = 0;

			if (dy1) dax_step = dx1 / (float)abs(dy1);
			if (dy2) dbx_step = dx2 / (float)abs(dy2);

			if (dy1) du1_step = du1 / (float)abs(dy1);
			if (dy1) dv1_step = dv1 / (float)abs(dy1);
			if (dy1) dw1_step = dw1 / (float)abs(dy1);

			if (dy2) du2_step = du2 / (float)abs(dy2);
			if (dy2) dv2_step = dv2 / (float)abs(dy2);
			if (dy2) dw2_step = dw2 / (float)abs(dy2);

			if (dy1)
			{
				for (int i = y1; i <= y2; i++)
				{
					int ax = x1 + (float)(i - y1) * dax_step;
					int bx = x1 + (float)(i - y1) * dbx_step;

					float tex_su = u1 + (float)(i - y1) * du1_step;
					float tex_sv = v1 + (float)(i - y1) * dv1_step;
					float tex_sw = w1 + (float)(i - y1) * dw1_step;

					float tex_eu = u1 + (float)(i - y1) * du2_step;
					float tex_ev = v1 + (float)(i - y1) * dv2_step;
					float tex_ew = w1 + (float)(i - y1) * dw2_step;

					if (ax > bx)
					{
						std::swap(ax, bx);
						std::swap(tex_su, tex_eu);
						std::swap(tex_sv, tex_ev);
						std::swap(tex_sw, tex_ew);
					}

					tex_u = tex_su;
					tex_v = tex_sv;
					tex_w = tex_sw;

					float tstep = 1.0f / ((float)(bx - ax));
					float t = 0.0f;

					for (int j = ax; j < bx; j++)
					{
						tex_u = (1.0f - t) * tex_su + t * tex_eu;
						tex_v = (1.0f - t) * tex_sv + t * tex_ev;
						tex_w = (1.0f - t) * tex_sw + t * tex_ew;
						if (tex_w > m_DepthBuffer[i * ScreenWidth() + j] )
						{
							Draw(j, i, spr->Sample(tex_u / tex_w, tex_v / tex_w));
							m_DepthBuffer[i * ScreenWidth() + j] = tex_w;
						}
						t += tstep;
					}

				}
			}

			dy1 = y3 - y2;
			dx1 = x3 - x2;
			dv1 = v3 - v2;
			du1 = u3 - u2;
			dw1 = w3 - w2;

			if (dy1) dax_step = dx1 / (float)abs(dy1);
			if (dy2) dbx_step = dx2 / (float)abs(dy2);

			du1_step = 0, dv1_step = 0;
			if (dy1) du1_step = du1 / (float)abs(dy1);
			if (dy1) dv1_step = dv1 / (float)abs(dy1);
			if (dy1) dw1_step = dw1 / (float)abs(dy1);

			if (dy1)
			{
				for (int i = y2; i <= y3; i++)
				{
					int ax = x2 + (float)(i - y2) * dax_step;
					int bx = x1 + (float)(i - y1) * dbx_step;

					float tex_su = u2 + (float)(i - y2) * du1_step;
					float tex_sv = v2 + (float)(i - y2) * dv1_step;
					float tex_sw = w2 + (float)(i - y2) * dw1_step;

					float tex_eu = u1 + (float)(i - y1) * du2_step;
					float tex_ev = v1 + (float)(i - y1) * dv2_step;
					float tex_ew = w1 + (float)(i - y1) * dw2_step;

					if (ax > bx)
					{
						std::swap(ax, bx);
						std::swap(tex_su, tex_eu);
						std::swap(tex_sv, tex_ev);
						std::swap(tex_sw, tex_ew);
					}

					tex_u = tex_su;
					tex_v = tex_sv;
					tex_w = tex_sw;

					float tstep = 1.0f / ((float)(bx - ax));
					float t = 0.0f;

					for (int j = ax; j < bx; j++)
					{
						tex_u = (1.0f - t) * tex_su + t * tex_eu;
						tex_v = (1.0f - t) * tex_sv + t * tex_ev;
						tex_w = (1.0f - t) * tex_sw + t * tex_ew;




						if (tex_w > m_DepthBuffer[i * ScreenWidth() + j])
						{
							Draw(j, i, spr->Sample(tex_u / tex_w, tex_v / tex_w));
							m_DepthBuffer[i * ScreenWidth() + j] = tex_w;
						}
						t += tstep;
					}
				}
			}
	
	}
	
	void Raster(vector<triangle> vecTrianglesToRaster, float r, float g, float b, int rastermode , olc::Sprite *spr) {
		//TODO documentation for rastermodes

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
					case 1:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, (float)ScreenHeight() - 1, 0.0f }, { 0.0f, -1.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 2:	nTrisToAdd = Triangle_ClipAgainstPlane({ 0.0f, 0.0f, 0.0f }, { 1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
					case 3:	nTrisToAdd = Triangle_ClipAgainstPlane({ (float)ScreenWidth() - 1, 0.0f, 0.0f }, { -1.0f, 0.0f, 0.0f }, test, clipped[0], clipped[1]); break;
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
						t.p[2].x, t.p[2].y, olc::Pixel(r, g, b));

					DrawTriangle(t.p[0].x, t.p[0].y,
						t.p[1].x, t.p[1].y,
						t.p[2].x, t.p[2].y, olc::GREY);
					/*float z1 = (triToRaster.p[0].z + triToRaster.p[1].z + triToRaster.p[2].z) / 3.0f;
					float x1 = (triToRaster.p[0].x + triToRaster.p[1].x + triToRaster.p[2].x) / 3.0f;
					float y1 = (triToRaster.p[0].y + triToRaster.p[1].y + triToRaster.p[2].y) / 3.0f;*/
				}
				break;
			case 3: 	// Draw the transformed, viewed, clipped, projected, sorted, clipped triangles
				for (auto& t : listTriangles)
				{
					FillTriangle(t.p[0].x, t.p[0].y,
						t.p[1].x, t.p[1].y,
						t.p[2].x, t.p[2].y, olc::Pixel(t.r, t.g, t.b));

					DrawTriangle(t.p[0].x, t.p[0].y,
						t.p[1].x, t.p[1].y,
						t.p[2].x, t.p[2].y, olc::GREY);
					/*float z1 = (triToRaster.p[0].z + triToRaster.p[1].z + triToRaster.p[2].z) / 3.0f;
					float x1 = (triToRaster.p[0].x + triToRaster.p[1].x + triToRaster.p[2].x) / 3.0f;
					float y1 = (triToRaster.p[0].y + triToRaster.p[1].y + triToRaster.p[2].y) / 3.0f;*/
				}
				break;
			case 4: //colour and text 
				for (auto& t : listTriangles)
				{
					//FillTriangle(t.p[0].x, t.p[0].y,t.p[1].x, t.p[1].y,t.p[2].x, t.p[2].y, olc::Pixel(t.r, t.g, t.b));
					DrawTriangle(t.p[0].x, t.p[0].y,
						t.p[1].x, t.p[1].y,
						t.p[2].x, t.p[2].y, olc::Pixel(t.r, t.g, t.b));


				}
				break;
			case 5: // text 
				for (auto& t : listTriangles)
				{
					//FillTriangle(t.p[0].x, t.p[0].y,t.p[1].x, t.p[1].y,t.p[2].x, t.p[2].y, olc::Pixel(t.r, t.g, t.b));
					DrawTexture(t.p[0].x, t.p[0].y, t.t[0].u, t.t[0].v, t.t[0].w,
						t.p[1].x, t.p[1].y, t.t[1].u, t.t[1].v, t.t[1].w,
						t.p[2].x, t.p[2].y, t.t[2].u, t.t[2].v, t.t[2].w, spr) ;


				} break;
			default: 	cout << "no raster" << endl;
			}


		}
	}


public:
	bool OnUserCreate() override
	{
		
	m_DepthBuffer = new float[ScreenWidth() * ScreenHeight()]; 

		//meshCube.loadFromObjectFile("cube.obj");

		e1 = entity(10,-2,10,"CompTable.obj",0,0,0, true);
		e1.setColour(10, 10, 150);
		e1.scale = 1;
		//objects.push_back(e1);
		
		entity e6 = entity(10, 0, 10, "MONKEYreturn.obj", 0, 2, 0, false);
		e6.setColour(10, 10, 150);
		e6.scale = 1;
		objects.push_back(e6);
		

		entity e2(0, -2, 0, "ground.obj", 0, 0, 0);
		e2.rasterFirst = true;
		e2.scale = 15;
		e2.setColour(10, 150, 10);
		e2.model.setColour(10, 150, 10);
		objects.push_back(e2);

		entity e4(-0.15, -0.15, -0, "doesnotexsist.obj", 0, 3.2/2, 0);
		e4.setColour(10, 150, 10);
		e4.scale = 1;
		e4.model.setColour(10, 10, 10);
		player.primary.weapon = e4;
		player.primary.type = 0;// full auto
		
		
		/*entity e5(10, -2, 0, "ground.obj", 0, 0, 0, true);
		e5.setColour(150, 150, 150);
		e5.scale = .04;
		e5.model.setColour(150, 150, 150);*/
		//player.secondary.weapon = e5;	
		player.secondary.deltaTime = .5;
		player.secondary.type = 1;// semi auto
		//selectedObject = e5;

		//objects.push_back(e5);

		matProj = Matrix_MakeProjection(90.0f, (float)ScreenHeight() / (float)ScreenWidth(), 0.1f, 1000.0f);
		
		tempSprite = std::make_unique<olc::Sprite>("./DSteelPlate.jpg");
		
		return true;
	}
	// game loop 
	bool OnUserUpdate(float fElapsedTime) override
	{

		vector<triangle>vecPoints;
		vector<triangle>vecSuper;
		vector<triangle>vecCursor;
		//input // 
		vec3d vForward = { 0,0,0 };
		vForward.x = 4.0f * cos(player.fYaw + 3.14 / 2) * fElapsedTime;	// Travel Along X-Axis
		vForward.z = 4.0f * sin(player.fYaw + 3.14 / 2) * fElapsedTime;	// Travel Along Z-Axis;;
		vec3d vSideways = { -vForward.z,0,vForward.x };
		{
			//testing key t for testing
			if (GetKey(olc::Key::T).bHeld) {
				currentNode->myMeshNode.Rotx+=0.1;
				
			}
			

			if (GetKey(olc::Key::TAB).bPressed) {
				//change selected object
				if (objects.size() > 0) {
					lastSelected = (lastSelected + 1) % objects.size();
					selectedObject = objects.at(lastSelected);
				}

			}
			if (GetKey(olc::Key::X).bPressed) {
				//change selected object
				vector<triangle> temp;
				for (auto& tri : selectedObject.model.tris) { // do triangle clip


					mat4x4 matRelocate = Matrix_MakeRelocate({ selectedObject.position.x, selectedObject.position.y, selectedObject.position.z }, selectedObject.Rotx, selectedObject.Roty, selectedObject.Rotz); //Matrix_MultiplyMatrix(matz, matx);


					triangle newTri;
					newTri.p[0] = MatrixMultiplyVector(tri.p[0], matRelocate);
					newTri.p[1] = MatrixMultiplyVector(tri.p[1], matRelocate);
					newTri.p[2] = MatrixMultiplyVector(tri.p[2], matRelocate);

					if (Triangle_InPlane(cursor.position, cursor.normal, newTri)) { // sub triangle 
						temp.push_back(tri);
					}
				}
				//selectedObject.model.tris = temp;
				selectedObject.model.tris.clear();
				for (auto& tri : temp) {

					selectedObject.model.tris.push_back(tri);
				}
			}
			if (GetKey(olc::Key::B).bPressed) {
				//set TOP mesh

				meshNode temp(selectedObject.model.tris, selectedObject.Rotx, selectedObject.Roty, selectedObject.Rotz, selectedObject.position);
				topTree = meshTree(temp);
				currentNode = &topTree;
				trees.push_back(topTree);
			}
			if (GetKey(olc::Key::V).bPressed) {
				currentNode = &topTree;
			}
			if (GetKey(olc::Key::N).bPressed && currentNode!= NULL) {
				//change selected object
				vector<triangle> topTemp;
				vector<triangle> temp;
				vector<vec3d>critPoints;
				for (auto& tri : currentNode->myMeshNode.tri) { // do triangle clip

					mat4x4 matRelocate = Matrix_MakeRelocate({ currentNode->myMeshNode.position.x, currentNode->myMeshNode.position.y, currentNode->myMeshNode.position.z }, currentNode->myMeshNode.Rotx, currentNode->myMeshNode.Roty, currentNode->myMeshNode.Rotz); //Matrix_MultiplyMatrix(matz, matx);


					triangle newTri;
					newTri.p[0] = MatrixMultiplyVector(tri.p[0], matRelocate);
					newTri.p[1] = MatrixMultiplyVector(tri.p[1], matRelocate);
					newTri.p[2] = MatrixMultiplyVector(tri.p[2], matRelocate);

					int critCheck = Triangle_InPlaneCritical(cursor.position, cursor.normal, newTri); // sub triangle 
					switch (critCheck) {

					case 1:
						critPoints.push_back(tri.p[0]);
						break;
					case 2:
						critPoints.push_back(tri.p[1]);
						break;
					case -1:
						critPoints.push_back(tri.p[1]);
						critPoints.push_back(tri.p[2]);
						break;
					case -2:
						critPoints.push_back(tri.p[0]);
						critPoints.push_back(tri.p[2]);
						break;
					case -3:
						critPoints.push_back(tri.p[0]);
						critPoints.push_back(tri.p[1]);

						break;
					case 3:
						critPoints.push_back(tri.p[2]);
						break;
					}
					if (critCheck != -999) { //if sub mesh
						tri.criticalPoints = critCheck;
						temp.push_back(tri);
					}
					else {//top mesh 
						topTemp.push_back(tri);
					}


				}
				//Create sub node  Triangle_InPlaneCritical
				meshNode subNode(temp, 0, 0, 0, {0,0,0}, critPoints);
				meshTree subTree = meshTree(subNode);
				currentNode->addSub(subTree);
				//currentNode->clearModel();
				currentNode->myMeshNode.tri = topTemp;
				currentNode = &(currentNode->subMeshNodes.back());
			}
			if (GetKey(olc::Key::M).bHeld && currentNode != NULL) {
				for (auto& tri : currentNode->myMeshNode.tri) { // do triangle clip

					mat4x4 matRelocate = Matrix_MakeRelocate({ selectedObject.position.x, selectedObject.position.y, selectedObject.position.z }, selectedObject.Rotx, selectedObject.Roty, selectedObject.Rotz); //Matrix_MultiplyMatrix(matz, matx);


					triangle newTri;
					newTri.p[0] = MatrixMultiplyVector(tri.p[0], matRelocate);
					newTri.p[1] = MatrixMultiplyVector(tri.p[1], matRelocate);
					newTri.p[2] = MatrixMultiplyVector(tri.p[2], matRelocate);

					if (Triangle_InPlaneCritical(cursor.position, cursor.normal, newTri) != -999) { // sub triangle 
						vecPoints.push_back(newTri);
					}
				}
			}

			if (GetKey(olc::Key::L).bHeld) {
				//rotate and move cursor
				//move
				if (GetKey(olc::Key::SHIFT).bHeld) {
					cursor.position += {0, 0, -1 * fElapsedTime};
				}
				else {
					//rotate
					cursor.normal += {0, 0, -1 * fElapsedTime};
				}
			}
			if (GetKey(olc::Key::K).bHeld) {
				//rotate and move cursor
				if (GetKey(olc::Key::SHIFT).bHeld) {
					cursor.position += {0, -1 * fElapsedTime, 0};
				}
				else {
					//rotate
					cursor.normal += {0, -1 * fElapsedTime, 0};
				}
			}
			if (GetKey(olc::Key::J).bHeld) {
				if (GetKey(olc::Key::SHIFT).bHeld) {
					cursor.position += { -1 * fElapsedTime, 0, 0};
				}
				else {
					//rotate
					cursor.normal += {-1 * fElapsedTime, 0, 0};
				}
			}
			if (GetKey(olc::Key::O).bHeld) {
				//rotate and move cursor
				//move
				if (GetKey(olc::Key::SHIFT).bHeld) {
					cursor.position += {0, 0, 1 * fElapsedTime};
				}
				else {
					//rotate
					cursor.normal += {0, 0, 1 * fElapsedTime};
				}
			}
			if (GetKey(olc::Key::I).bHeld) {
				//rotate and move cursor
				if (GetKey(olc::Key::SHIFT).bHeld) {
					cursor.position += {0, 1 * fElapsedTime, 0};
				}
				else {
					//rotate
					cursor.normal += {0, 1 * fElapsedTime, 0};
				}
			}
			if (GetKey(olc::Key::U).bHeld) {
				if (GetKey(olc::Key::SHIFT).bHeld) {
					cursor.position += { 1 * fElapsedTime, 0, 0};
				}
				else {
					//rotate
					cursor.normal += {1 * fElapsedTime, 0, 0};
				}
			}
			
			if (GetKey(olc::Key::A).bHeld) {
				player.vCamera -= vSideways;
			}
			if (GetKey(olc::Key::D).bHeld) {
				player.vCamera += vSideways;
			}
			if (GetKey(olc::Key::K2).bPressed) {
				//draw gun
				player.drawWeapon(2);
			}
			if (GetKey(olc::Key::K1).bPressed) {
				//draw gun
				player.drawWeapon(1);
			}
			if (GetKey(olc::Key::W).bHeld) player.vCamera += vForward;
			if (GetKey(olc::Key::S).bHeld) player.vCamera -= vForward;
			if (GetKey(olc::Key::SPACE).bPressed) {
				if (player.vCamera.y < 1.0f) {
					player.speed.y = 0;
					player.speed += {0.0f, 70.0f * fElapsedTime, 0.0f};
				}
			}
			if (GetKey(olc::Key::ESCAPE).bPressed) return (0);
			//mouse control
			lockMouse = IsFocused();

			if (lockMouse) {
				SetCursor(NULL);
				SetCursor(LoadCursor(NULL, NULL));
				ShowCursor(false);

				POINT point;
				GetCursorPos(&point);
				player.fYaw += (point.x - 200) * sensitivity;
				player.fPitch += (point.y - 200) * sensitivity;
				SetCursorPos(200, 200);
			}
			else
			{
				ShowCursor(true);
			}

			player.wInHand.lastTime += fElapsedTime;

			if (GetMouse(0).bPressed && player.weaponSelected != 0 && player.wInHand.deltaTime <= fElapsedTime + player.wInHand.lastTime) {
				objects.push_back(player.fire());
				player.wInHand.lastTime = 0;
			}
			if (GetMouse(0).bHeld && player.weaponSelected != 0 && player.wInHand.type == 0) {

				if (player.wInHand.deltaTime <= fElapsedTime + player.wInHand.lastTime) {
					player.wInHand.lastTime = 0;
					objects.push_back(player.fire());
				}

			}
		}

		// world changes 
		player.speed -= {0.0f, 9.8f * fElapsedTime, 0.0f}; // gravity 

		player.vCamera += player.speed / 4;

		if (player.vCamera.y < 0.25f) { // player in floor
			player.speed.y = 0;
			player.vCamera.y = 0.25f;
		}


		if (player.fPitch > 3.14 / 2) {
			player.fPitch = 3.14 / 2;
		}
		else if (player.fPitch < -3.14 / 2) {
			player.fPitch = -3.14 / 2;
		}

		//clear depth buffer
		memset(m_DepthBuffer, 0,ScreenWidth()* ScreenHeight() * sizeof(float));

		/// position and draw objects//
		// drawing update 
		FillRect(0, 0, ScreenWidth(), ScreenHeight(), olc::Pixel(140, 140, 255));

		// get camera matrix // learn more about later 
		vec3d vUp = { 0,-1,0 };
		vec3d vTarget = { 0,0,1 };
		mat4x4 matCamerRoty = Matrix_MakeRotationY(player.fYaw);
		mat4x4 matCamerRotx = Matrix_MakeRotationX(player.fPitch);
		mat4x4 matCamerRot = Matrix_MultiplyMatrix(matCamerRotx, matCamerRoty);
		player.vLookDir = MatrixMultiplyVector(vTarget, matCamerRot);

		vTarget = player.vCamera + player.vLookDir;
		mat4x4 matCamera = Matrix_PointAt(player.vCamera, vTarget, vUp);
		mat4x4 matView = Matrix_QuickInverse(matCamera);
		mat4x4 matWorld;
		matWorld = Matrix_MakeTranslation(0.0f, 0.0f, 0.0f);
	



		vector<triangle>vecFirst;
		vector<triangle>vecTrianglesInWorld; // holds triangles in world 
		vector<triangle>vecHitbox; // holds triangles in world 

		vector<entity>points;

		
		//creates the points
		for (auto& obj : objects) { //object interactions for every object create a box to repersent it as a point

			obj.position += obj.speed * fElapsedTime;

			//create origin points 
			entity point(obj.position.x, obj.position.y, obj.position.z, "cube.obj", 0, 0, 0);
			point.setColour(150, 10, 10);
			point.scale = .1;
			point.model.setColour(150, 10, 10);
			points.push_back(point);
		}

		//loads in the points
		for (auto& obj : points) { // load in all points
		mat4x4 matRelocate = Matrix_MakeRelocate({ obj.position.x, obj.position.y, obj.position.z }, obj.Rotx, obj.Roty, obj.Rotz); //Matrix_MultiplyMatrix(matz, matx);
			// loads in triangle of objects
			for (auto& tri : obj.model.tris) {
				triangle newTri = transformTriangleIntoWorld(obj.scale,matRelocate, tri);
				newTri.r = obj.model.r; newTri.g = obj.model.g; newTri.b = obj.model.b;
				vecPoints.push_back(newTri);// send it to world 	
			}	
		}
		
		//load in all entities
		for (auto& obj : objects ) { // load in all enitites 
			vector<triangle>vecTemp;
			mat4x4 matRelocate = Matrix_MakeRelocate({ obj.position.x, obj.position.y, obj.position.z }, obj.Rotx, obj.Roty, obj.Rotz); //Matrix_MultiplyMatrix(matz, matx);

			// loads in triangle of objects
			for (auto& tri : obj.model.tris) { 
				triangle newTri = transformTriangleIntoWorld(obj.scale, matRelocate, tri);
				newTri.r = obj.model.r; newTri.g = obj.model.g; newTri.b = obj.model.b;
				
				if (obj.rasterFirst) {
					vecFirst.push_back(newTri);
				}
				else {
					vecTrianglesInWorld.push_back(newTri);// send it to world 
					vecTemp.push_back(newTri);
				}
			}
			

		}
		
		//load in the cursor
		mesh tempMesh = cursor.cursorPoint();
		 mat4x4 matT = Matrix_MakeTranslation(cursor.position.x, cursor.position.y, cursor.position.z);
		 mat4x4 matRelocate = matT;

		for (auto& tri : tempMesh.tris) {
			
			triangle newTri = transformTriangleIntoWorld(.01, matRelocate, tri);
			newTri.r = 255; newTri.g = 255; newTri.b = 255;
			vecCursor.push_back(newTri);// send it to world 	
		}
		tempMesh = cursor.cursorArrow();
		matT = Matrix_MakeTranslation(cursor.position.x + cursor.normal.x / 2, cursor.position.y + cursor.normal.y / 2, cursor.position.z + cursor.normal.z / 2);
		matRelocate = matT;

		for (auto& tri : tempMesh.tris) {
			triangle newTri = transformTriangleIntoWorld(.05, matRelocate, tri);
			newTri.r = 0; newTri.g = 255; newTri.b = 0;
			vecCursor.push_back(newTri);// send it to world 	
		}

		
		

		entity tempEntity = cursor.cursorPlane();
		tempMesh = tempEntity.model;
		matRelocate = Matrix_MakeRelocate({ tempEntity.position.x, tempEntity.position.y, tempEntity.position.z }, tempEntity.Rotx, tempEntity.Roty, tempEntity.Rotz); //Matrix_MultiplyMatrix(matz, matx);
	

		for (auto& tri : tempMesh.tris) {
			triangle newTri = transformTriangleIntoWorld(tempEntity.scale, matRelocate, tri);
			newTri.r = 40; newTri.g = 155; newTri.b = 40;
			vecCursor.push_back(newTri);// send it to world 	
		}


		tempEntity = cursor.cursorInvertedPlane();
		tempMesh = tempEntity.model;
		matRelocate = Matrix_MakeRelocate({ tempEntity.position.x, tempEntity.position.y, tempEntity.position.z }, tempEntity.Rotx, tempEntity.Roty, tempEntity.Rotz);
		
		for (auto& tri : tempMesh.tris) {
			triangle newTri = transformTriangleIntoWorld(tempEntity.scale, matRelocate, tri);
			newTri.r = 155; newTri.g = 40; newTri.b = 40;
			vecCursor.push_back(newTri);// send it to world 	
		}

		//load in in hand object
		entity obj = player.inHand;

		matRelocate = Matrix_MakeRelocate({ obj.position.x, obj.position.y, obj.position.z }, obj.Rotx, obj.Roty, obj.Rotz);
		matRelocate = Matrix_Relocate(matRelocate, { player.vCamera.x, player.vCamera.y, player.vCamera.z }, player.fPitch, player.fYaw, 0);
		for (auto& tri : obj.model.tris) {

			triangle newTri = transformTriangleIntoWorld(obj.scale, matRelocate, tri);
			newTri.r = obj.model.r; newTri.g = obj.model.g; newTri.b = obj.model.b;

			if (obj.rasterFirst) {
				vecFirst.push_back(newTri);
			}
			else {
				vecTrianglesInWorld.push_back(newTri);// send it to world 

			}
		}
		
		
	
		vector<triangle>vecTrianglesToRaster;// holds visable triangles in world to be drawn 
		vector<triangle>vecTrianglesToRasterFirst;// holds visable triangles in world to be drawn first 
		vector<triangle>vecHitboxRaster;// holds hitbox to be drawn 
		vector<triangle>vecPointsRaster;// holds points to be drawn 
		vector<triangle>vecCursorRaster;// holds points to be drawn 

		// transform world objects to projected space 
		vecTrianglesToRaster = transformTriangleIntoProjection(vecTrianglesInWorld, matWorld, matView, matProj,false,player.vCamera);
		vecTrianglesToRasterFirst = transformTriangleIntoProjection(vecFirst, matWorld, matView, matProj,false, player.vCamera);
		vecPointsRaster = transformTriangleIntoProjection(vecPoints, matWorld, matView, matProj, false, player.vCamera);
		vecCursorRaster = transformTriangleIntoProjection(vecCursor, matWorld, matView, matProj, false, player.vCamera);

		

		// send it to world 
		Raster(vecTrianglesToRasterFirst, 10,150,10,0, NULL);
		Raster(vecTrianglesToRaster, 10,10,150,3, NULL);
		//Raster(vecPointsRaster, 150, 10, 10, 2, NULL);
		Raster(vecCursorRaster, 255, 10, 10, 3, NULL);
		
		if (currentNode != NULL) {
			vector<triangle>TopNodeTriList;

			matRelocate = Matrix_MakeRelocate({ topTree.myMeshNode.position.x, topTree.myMeshNode.position.y, topTree.myMeshNode.position.z }, topTree.myMeshNode.Rotx, topTree.myMeshNode.Roty, topTree.myMeshNode.Rotz);
			for (auto& tri : topTree.myMeshNode.tri) {
				triangle newTri = transformTriangleIntoWorld(topTree.myMeshNode.Scale, matRelocate, tri);
				newTri.r = 155; newTri.g = 40; newTri.b = 40;
				TopNodeTriList.push_back(newTri);// send it to world 	
			}

			matRelocate = Matrix_MakeRelocate({ currentNode->myMeshNode.position.x, currentNode->myMeshNode.position.y, currentNode->myMeshNode.position.z }, currentNode->myMeshNode.Rotx, currentNode->myMeshNode.Roty, currentNode->myMeshNode.Rotz);
			vector<triangle>CurrentNodeTriList = triangleRelocate(currentNode->myMeshNode.tri, matRelocate, currentNode->myMeshNode.Scale);
			vector<triangle>TopNodeTriListRaster = transformTriangleIntoProjection(TopNodeTriList, matWorld, matView, matProj, false, player.vCamera);
			vector<triangle>CurrentNodeTriListRaster = transformTriangleIntoProjection(CurrentNodeTriList, matWorld, matView, matProj, false, player.vCamera);
			Raster(TopNodeTriListRaster, 0, 0, 255, 1, NULL);
			Raster(CurrentNodeTriListRaster, 0, 255, 0, 1, NULL);
			
			

			matWorld = Matrix_MakeTranslation(4.0f, -1.0f, 0.0f);
			vector<triangle> treeTri = MeshTreeToWorld(topTree);
			vector<triangle> treeTriRaster = transformTriangleIntoProjection(treeTri, matWorld, matView, matProj, false, player.vCamera);
			Raster(treeTriRaster, 0, 255, 0, 2, NULL);
			
		}

		matWorld = Matrix_MakeTranslation(0.0f, 0.0f, 0.0f);

		vector<triangle> tex = e1.model.tris;
		vector<triangle>testTriList = triangleRelocate(tex, {0,-2,0},0,0,0, 1);
		vector<triangle> texRast = transformTriangleIntoProjection(testTriList, matWorld, matView, matProj, false, player.vCamera);
		Raster(texRast, 0, 255, 0, 5, tempSprite.get());
		
	}
};

int main()
{
	FPSengine demo;
	if (demo.Construct(1440, 900, 1, 1, 1))
		
		demo.Start();
	return 0;
}
