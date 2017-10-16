#include <iostream>
#include <sstream> 
#include <bitset>
#include <functional>
#include <set>
#include <map>
#include <unordered_map>
#include "myMath.h"
#include "BoundBox.h"

#ifdef MATHFUNCSDLL_EXPORTS
#define MATHFUNCSDLL_API __declspec(dllexport) 
#else
#define MATHFUNCSDLL_API __declspec(dllimport) 
#endif
using namespace myMath;
using namespace std;
class ModelPCA
{
public:
	Point3D O, basis0, basis1, basis2;
	double eigValues[3];

	ModelPCA() {};
	~ModelPCA() {};

	char* ToString()
	{
		stringstream strStream;

		strStream << "O:" << O[0] << "," << O[1] << "," << O[2] << "\n";
		strStream << "eigveter:\n";
		strStream << basis0.x << " " << basis0.y << " " << basis0.z << " eigvalue: " << eigValues[0] << "\n";
		strStream << basis1.x << " " << basis1.y << " " << basis1.z << " eigvalue: " << eigValues[1] << "\n";
		strStream << basis2.x << " " << basis2.y << " " << basis2.z << " eigvalue: " << eigValues[2] << "\n";

		string str = strStream.str();
		std::vector<char> cstr(str.c_str(), str.c_str() + str.size() + 1);
		char* ca = &cstr[0];
		return ca;
	}
};

class VBOobj
{
public:
	unsigned int* VBO;
	unsigned int VBOsize;

	VBOobj();
	VBOobj(const VBOobj& V);
	virtual ~VBOobj();

	VBOobj(unsigned int s);

	bool releaseVBO();
};

class Model :public VBOobj
{
public:
	
	TransformationMatrix tm;

	
	bool showmesh;
	bool dynamicshowmesh;

	ModelPCA PAxis;

	spvector<ModelTriangle> tarr;
	spvector<triangleinfo> tarrinfo;

private:
	vector<Point3D> sourceP; //點資料

public:
	Model();
	Model(const Model& M);
	Model(Model&& M);
	Model& operator=(const Model& M) = delete;
	virtual ~Model();

	explicit Model(const spvector<baseTriangle>& V);
	explicit Model(const string& fileName);

	static Model LoadObjASCI(const string& fileName);
private:
	void SetModel(float* xyz, size_t numofTri);
public:
	void SaveStlBinary(const string& fileName) const;
	void SavePointCloudASCII(const string& fileName) const;
	void SavePointCloudBinary(const string& fileName) const;
	const ModelTriangle& operator[](unsigned int ind) const
	{
		if (ind < tarr.size())
		{
			return *(tarr[ind]);
		}
		return *(tarr[0]);
	}
	ModelTriangle& operator[](unsigned int ind)
	{
		if (ind < tarr.size())
		{
			return *(tarr[ind]);
		}
		return *(tarr[0]);
	}

	void Add_tri(spvector<baseTriangle> bt);

	void ApplyTMat();
	void RebuildOctCloud();
	void Triangle_update();
	void SetMinTriSize();
	void BuildAABB();
	void CalModelPCA();


	//檢測退化三角形相關
	tuple<bool, char, char, double> checkDegenerationTriangle(const ModelTriangle& T, double llimit) const;
	spvector<baseTriangle> searchDegenerationTriangle(double llimit = 1e-6) const;
	void removeDegenerationTriangle(double llimit = 1e-6);

	//三角形相交
	static pair<bool, QSegment3D*> checkTriangleCross(const ModelTriangle& T1, const ModelTriangle& T2);
	//三角形自相交
	static pair<bool, QSegment3D*> checkTriangleCross_onesame(const ModelTriangle& T1, const ModelTriangle& T2);

	TriCrossResult searchTriangleCross();
	TriCrossResult searchTriangleCross(const Model* rhs, const vectors<unsigned int>& testind) const;
	TriSeparateResult separateTriangle(const TriCrossResult& x);

	static spvector<Segment3D> searchTriangleCross_fordisplay(const TriCrossResult& x);
	static spvectors<Segment3D> searchTriangleCross_fordisplay2(const TriCrossResult& x);
	static spvector<baseTriangle> subtri_fordisplay(const TriSeparateResult& x);
	spvectors<Segment3D> searchEdgecontour_fordisplay() const;
	static spvectors<Segment3D> searchEdgecontour_fordisplay(const Model& m, const vectors<Edgeind>& vei);
	static void repairtsr(vector<reference_wrapper<spvector<ModelTriangle>>>& modelarr, vector<Model::TriSeparateResult>& tsr);

	spvector<baseTriangle> searchEdgeTriangle() const;
	spvector<Segment3D> searchEdge() const;
	vector<Edgeind> searchEdgeind() const;
	vectors<Edgeind> searchEdgecontour() const;
	vectors<Edgeind> searchEdgecontour(vector<Edgeind>& sourceEdge) const;

	sp<baseTriangle> searchTriangle(const Point3D& P) const;
	unsigned int searchTriangleind(const Point3D& P) const;
	spvector<baseTriangle> getTriangles(const vector<unsigned int>& ind) const;
	spvectors<baseTriangle> getTriangles(const vectors<unsigned int>& ind) const;

	double calcMeanCurvature(const Point3DNode* pn) const;
	double calcGaussCurvature(const Point3DNode* pn) const;

	std::pair<double, double> calcCurvature(const Point3DNode* pn) const;
	bool signCurvature(const Point3DNode* pn) const;
	double calcCurvature(const Point3D& P) const;

	Point3D PointinModelCoor(const Point3D& P) const;
	Vector3D VectorinModelCoor(const Vector3D& V) const;
	Point3D PointinWorldCoor(const Point3D& P) const;
	Vector3D VectorinWorldCoor(const Vector3D& V) const;

	spvector<baseTriangle> shiftTriangle(double D, unsigned int dir) const;
	spvector<baseTriangle> shiftTriangle2(double D, unsigned int dir, Vector3D n) const;
	spvector<baseTriangle> shiftTriangle3(const map<Point3DNode*, double>& mpnd, double shiftd, unsigned int dir) const;
	spvector<baseTriangle> shiftTriangle4(const map<Point3DNode*, double>& mpnd, double D, unsigned int dir, Vector3D n) const;




	pair<Point3D, Point3D> Intersection(const Line3D& i) const;
	sp<Model> filpmodel() const;
};