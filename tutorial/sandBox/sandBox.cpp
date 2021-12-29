#include "sandBox.h"
#include "igl/edge_flaps.h"
#include "igl/collapse_edge.h"
#include "Eigen/dense"
#include <functional>
#include <Eigen/Core>
#include "igl/opengl/ViewerCore.h"
#include "igl/opengl/glfw/renderer.h"
#include "igl/decimate.h"
#include "igl/writeOBJ.h"
#include <igl/circulation.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/parallel_for.h>
#include <igl/opengl/glfw/Viewer.h>
#include <set>
#include "simplifier.h"
#include "igl/swept_volume_bounding_box.h"


SandBox::SandBox() 
	: 
	objectsData{ new std::vector<ObjectData*>{} }, 
	linksCounter{ 0 }, destIndex{ 0 }, 
	firstLinkIndex{ 1 }, 
	lastLinkIndex{ 1 }, 
	isCCD{ false },
	xVelocity{-0.005},
	yVelocity{ 0 }
{

}

SandBox::~SandBox()
{
    delete objectsData;
}

void SandBox::OnNewMeshLoad()
{
	if (linksCounter > 0)
	{
		data().MyTranslate(Eigen::Vector3d(0, 1.6, 0), true);
		parents.push_back(linksCounter);
	}
	else
	{
		parents.push_back(-1);
	}

	data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(0, 0, 1));
	data().show_overlay_depth = false;
	data().point_size = 10;
	data().line_width = 2;
	data().set_visible(true, 2);
	data().show_faces = 3;
	data().show_texture = 2;

	data().dirty = 157; //this line prevents texture coordinates
	data().MyScale(Eigen::Vector3d(1, 1, 1));
	data().SetCenterOfRotation(Eigen::Vector3d(0, -0.8, 0));

	Eigen::MatrixXd axis(7, 3);
	axis << 0, -0.8, 0, 0, -0.8, 1.6, 0, -0.8, -1.6, 0, 0.8, 0, 0, -2.4, 0, -1.6, -0.8, 0, 1.6, -0.8, 0;
	data().add_points(axis, Eigen::RowVector3d(0, 0, 1));
 	data().add_edges(axis.row(0), axis.row(1), Eigen::RowVector3d(0, 0, 1));
	data().add_edges(axis.row(0), axis.row(2), Eigen::RowVector3d(0, 0, 1));
	data().add_edges(axis.row(0), axis.row(3), Eigen::RowVector3d(0, 1, 0));
	data().add_edges(axis.row(0), axis.row(4), Eigen::RowVector3d(0, 1, 0));
	data().add_edges(axis.row(0), axis.row(5), Eigen::RowVector3d(1, 0, 0));
	data().add_edges(axis.row(0), axis.row(6), Eigen::RowVector3d(1, 0, 0));

	ObjectData* od = new ObjectData();
	InitObjectData(*od, data().V, data().F);
	objectsData->push_back(od);
	data().AddBoundingBox(od->tree->m_box, Eigen::RowVector3d(0, 1, 0));

 	linksCounter++;
	lastLinkIndex = linksCounter;
}

void SandBox::Init(const std::string &config)
{
	std::string item_name;
	std::ifstream nameFileout;
	nameFileout.open(config);

	if (!nameFileout.is_open())
	{
		std::cout << "Can't open file " << config << std::endl;
	}
	else
	{
		while (nameFileout >> item_name)
		{
			std::cout << "openning " << item_name << std::endl;
			load_mesh_from_file(item_name);

			if (item_name.find("sphere") != std::string::npos)
			{
				parents.push_back(-1);
				data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(0, 0, 1));
				data().show_overlay_depth = false;
				data().point_size = 10;
				data().line_width = 2;
				data().set_visible(false, 1);
			}
			else
			{
				if (linksCounter > 0)
				{
					data().MyTranslate(Eigen::Vector3d(0, 1.6, 0), true);
					parents.push_back(linksCounter);
				}
				else
				{
					parents.push_back(-1);
				}

				data().add_points(Eigen::RowVector3d(0, 0, 0), Eigen::RowVector3d(0, 0, 1));
				data().show_overlay_depth = false;
				data().point_size = 10;
				data().line_width = 2;
				data().set_visible(false, 1);

				data().dirty = 157; //this line prevents texture coordinates
				data().MyScale(Eigen::Vector3d(1, 1, 1));
				data().SetCenterOfRotation(Eigen::Vector3d(0, -0.8, 0));

				Eigen::MatrixXd axis(7, 3);
				axis << 
					   0, -0.8,    0, 
					   0, -0.8,  1.6, 
					   0, -0.8,	-1.6, 
					   0,  0.8,	   0, 
					   0, -2.4,    0, 
					-1.6, -0.8,    0, 
					 1.6, -0.8,	   0;

				data().add_points(axis, Eigen::RowVector3d(0, 0, 1));
				data().add_edges(axis.row(0), axis.row(1), Eigen::RowVector3d(0, 0, 1));
				data().add_edges(axis.row(0), axis.row(2), Eigen::RowVector3d(0, 0, 1));
				data().add_edges(axis.row(0), axis.row(3), Eigen::RowVector3d(0, 1, 0));
				data().add_edges(axis.row(0), axis.row(4), Eigen::RowVector3d(0, 1, 0));
				data().add_edges(axis.row(0), axis.row(5), Eigen::RowVector3d(1, 0, 0));
				data().add_edges(axis.row(0), axis.row(6), Eigen::RowVector3d(1, 0, 0));

				ObjectData* od = new ObjectData();
				InitObjectData(*od, data().V, data().F);
				objectsData->push_back(od);
				data().AddBoundingBox(od->tree->m_box, Eigen::RowVector3d(0, 1, 0));

				linksCounter++;
			}
		}
		lastLinkIndex = linksCounter;
		nameFileout.close();
	}
	data().set_colors(Eigen::RowVector3d(0.9, 0.1, 0.1));
	isActive = false;
}

void SandBox::InitObjectData(ObjectData& od, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    
    od.V = new Eigen::MatrixXd(V);
    od.F = new Eigen::MatrixXi(F);
    od.E = new Eigen::MatrixXi();
    od.EF = new Eigen::MatrixXi();
    od.EI = new Eigen::MatrixXi();
    od.EMAP = new Eigen::VectorXi();
    od.Q = new PriorityQueue();
    od.num_collapsed = 0;
    od.tree = new igl::AABB<Eigen::MatrixXd, 3>{};

    igl::edge_flaps(*od.F, *od.E, *od.EMAP, *od.EF, *od.EI);

    od.tree->init(*od.V, *od.F);

    od.C = new Eigen::MatrixXd(od.E->rows(), od.V->cols());

    ComputeNormals(od, data());

    ComputeQMatrices(od);

    ComputePriorityQueue(od);
}

void SandBox::ReInitObjectData(ObjectData& od)
{
    Eigen::MatrixXd V = *od.V; 
    Eigen::MatrixXi F = *od.F;

    // Remove old object data
    ObjectData* odToRemove = objectsData->at(selected_data_index);
    ClearObjectData(*odToRemove);
    delete odToRemove;

    // Add new object data after collapsing edges 
    ObjectData* odToAdd = new ObjectData();
    InitObjectData(*odToAdd, V, F);
    objectsData->at(selected_data_index) = odToAdd;
}

void SandBox::ClearObjectData(ObjectData& od)
{
    delete od.V;
    delete od.F;
    delete od.EMAP;
    delete od.E;
    delete od.EF;
    delete od.EI;
    delete od.Q;
    delete od.C;
    delete od.F_NORMALS;

    std::for_each(od.QMATRICES.begin(), od.QMATRICES.end(), [](Eigen::Matrix4d* m) -> void { delete m; });
    od.QMATRICES.clear();
}

void SandBox::Simplify()
{
    ObjectData* objectToSimplify = objectsData->at(selected_data_index);
    int num_to_collapse = std::ceil(0.05 * objectToSimplify->Q->size());
    Simplify(num_to_collapse, *objectToSimplify);
}

void SandBox::Simplify(int num_to_collapse, ObjectData& od)
{

    if (!od.Q->empty())
    {
        bool something_collapsed = false;
        for (int j = 0; j < num_to_collapse; j++)
        {
              if (!collapse_edge(od))
            {
                break;
            }
            something_collapsed = true;
            od.num_collapsed++;
        }
         
        if (something_collapsed)
        {
            data().clear();
            data().set_mesh(*od.V, *od.F);
            data().set_face_based(true);
            data().dirty = 157; //this line prevents texture coordinates
            ReInitObjectData(od);
        }
    }
}

void SandBox::MoveTo(double x, double y)
{
    data().TranslateInSystem(GetRotation(), Eigen::Vector3d(x, 0, 0));
    data().TranslateInSystem(GetRotation(), Eigen::Vector3d(0, y, 0));
    WhenTranslate();
}

void SandBox::Animate()
{
	if (isActive)
	{
		
		if (isCCD)
		{
			IK_CCD();
		}
		else
		{
			IK_FABRIK();
		}
	}
}

bool SandBox::ObjectsCollide(igl::AABB<Eigen::MatrixXd, 3>* firstTree, igl::AABB<Eigen::MatrixXd, 3>* secondTree)
{
	if (firstTree == nullptr || secondTree == nullptr)
	{
		return false;
	}
 	if (BoxesIntersect(firstTree->m_box, secondTree->m_box)) {
		if (firstTree->is_leaf() && secondTree->is_leaf())
		{
			data_list.at(0).AddBoundingBox(firstTree->m_box, Eigen::RowVector3d(0, 1, 0));
			data_list.at(1).AddBoundingBox(secondTree->m_box, Eigen::RowVector3d(0, 1, 0));
			return true;
		}
		else if (!firstTree->is_leaf() && secondTree->is_leaf())
		{
			return	ObjectsCollide(firstTree->m_left, secondTree) ||
					ObjectsCollide(firstTree->m_right, secondTree);
		}
		else if (firstTree->is_leaf() && !secondTree->is_leaf())
		{

			return  ObjectsCollide(firstTree, secondTree->m_left) ||
					ObjectsCollide(firstTree, secondTree->m_right);
		}
		else if (!firstTree->is_leaf() && !secondTree->is_leaf())
		{
			return  ObjectsCollide(firstTree->m_left, secondTree->m_left) ||
					ObjectsCollide(firstTree->m_left, secondTree->m_right) ||
					ObjectsCollide(firstTree->m_right, secondTree->m_left) ||
					ObjectsCollide(firstTree->m_right, secondTree->m_right);
		}
	}
	return false;
}

bool SandBox::BoxesIntersect(Eigen::AlignedBox <double, 3>& firstBox, Eigen::AlignedBox <double, 3>& secondBox)
{
	Eigen::Matrix4d firstTrans = data_list.at(0).MakeTransd();
	Eigen::Matrix4d secondTrans = data_list.at(1).MakeTransd();
	 
	Eigen::Vector3d firstBoxCenter = firstBox.center();
	Eigen::Vector3d secondBoxCenter = secondBox.center();

	Eigen::Vector4d firstCenter = firstTrans * Eigen::Vector4d(firstBoxCenter(0), firstBoxCenter(1), firstBoxCenter(2), 1);
	Eigen::Vector4d secondCenter = secondTrans * Eigen::Vector4d(secondBoxCenter(0), secondBoxCenter(1), secondBoxCenter(2), 1);

	Eigen::Vector3d C0(firstCenter(0), firstCenter(1), firstCenter(2));
	Eigen::Vector3d C1(secondCenter(0), secondCenter(1), secondCenter(2));

	Eigen::Vector3d D = C1 - C0;

	Eigen::Matrix3d A = data_list.at(0).GetRotation();
	Eigen::Matrix3d B = data_list.at(1).GetRotation();

	Eigen::Matrix3d A_matrix;
	A_matrix << A(0, 0), A(1, 0), A(2, 0),
				A(0, 1), A(1, 1), A(2, 1),
				A(0, 2), A(1, 2), A(2, 2);

	Eigen::Matrix3d B_matrix;
	B_matrix << B(0, 0), B(1, 0), B(2, 0),
				B(0, 1), B(1, 1), B(2, 1),
				B(0, 2), B(1, 2), B(2, 2);

	Eigen::RowVector3d a;
	a << firstBox.sizes()(0) / 2, firstBox.sizes()(1) / 2, firstBox.sizes()(2) / 2;

	Eigen::RowVector3d b;
	b << secondBox.sizes()(0) / 2, secondBox.sizes()(1) / 2, secondBox.sizes()(2) / 2;

	Eigen::Matrix3d C = A.transpose() * B;

	double R0, R1, R;

	// L = A0
	R0 = a(0);
	R1 = b(0) * abs(C(0, 0)) + b(1) * abs(C(0, 1)) + b(2) * abs(C(0, 2));
	R = (A_matrix.row(0) * D).norm();
	if (R > R0 + R1) return false;

	// L = A1
	R0 = a(1);
	R1 = b(0) * abs(C(1, 0)) + b(1) * abs(C(1, 1)) + b(2) * abs(C(1, 2));
	R = (A_matrix.row(1) * D).norm();
	if (R > R0 + R1) return false;

	// L = A2
	R0 = a(2);
	R1 = b(0) * abs(C(2, 0)) + b(1) * abs(C(2, 1)) + b(2) * abs(C(2, 2));
	R = (A_matrix.row(2) * D).norm();
	if (R > R0 + R1) return false;

	// L = B0
	R0 = a(0) * abs(C(0, 0)) + a(1) * abs(C(1, 0)) + a(2) * abs(C(2, 0));
	R1 = b(0);
	R = (B_matrix.row(0) * D).norm();
	if (R > R0 + R1) return false;

	// L = B1
	R0 = a(0) * abs(C(0, 1)) + a(1) * abs(C(1, 1)) + a(2) * abs(C(2, 1));
	R1 = b(1);
	R = (B_matrix.row(1) * D).norm();
	if (R > R0 + R1) return false;

	// L = B2
	R0 = a(0) * abs(C(0, 2)) + a(1) * abs(C(1, 2)) + a(2) * abs(C(2, 2));
	R1 = b(2);
	R = (B_matrix.row(2) * D).norm();
	if (R > R0 + R1) return false;

	// L = A0 x B0
	R0 = a(1) * abs(C(2, 0)) + a(2) * abs(C(1, 0));
	R1 = b(1) * abs(C(0, 2)) + b(2) * abs(C(0, 1));
	R = (C(1, 0) * A_matrix.row(2) * D - C(2, 0) * A_matrix.row(1) * D).norm();
	if (R > R0 + R1) return false;

	// L = A0 x B1
	R0 = a(1) * abs(C(2, 1)) + (a(2))*abs(C(1, 1));
	R1 = b(0) * abs(C(0, 2)) + b(2) * abs(C(0, 0));
	R = (C(1, 1) * A_matrix.row(2) * D - C(2, 1) * A_matrix.row(1) * D).norm();
	if (R > R0 + R1) return false;

	// L = A0 x B2
	R0 = a(1) * abs(C(2, 2)) + a(2) * abs(C(1, 2));
	R1 = b(0) * abs(C(0, 1)) + b(1) * abs(C(0, 0));
	R = (C(1, 2) * A_matrix.row(2) * D - C(2, 2) * A_matrix.row(1) * D).norm();
	if (R > R0 + R1) return false;

	// L = A1 x B0
	R0 = a(0) * abs(C(2, 0)) + a(2) * abs(C(0, 0));
	R1 = b(1) * abs(C(1, 2)) + b(2) * abs(C(1, 1));
	R = (C(2, 0) * A_matrix.row(0) * D - C(0, 0) * A_matrix.row(2) * D).norm();
	if (R > R0 + R1) return false;

	// L = A1 x B1
	R0 = a(0) * abs(C(2, 1)) + a(2) * abs(C(0, 1));
	R1 = b(0) * abs(C(1, 2)) + b(2) * abs(C(1, 0));
	R = (C(2, 1) * A_matrix.row(0) * D - C(0, 1) * A_matrix.row(2) * D).norm();
	if (R > R0 + R1) return false;

	/// L = A1 x B2
	R0 = a(0) * abs(C(2, 2)) + a(2) * abs(C(0, 2));
	R1 = b(0) * abs(C(1, 1)) + b(1) * abs(C(1, 0));
	R = (C(2, 2) * A_matrix.row(0) * D - C(0, 2) * A_matrix.row(2) * D).norm();
	if (R > R0 + R1) return false;

	// L = A2 x B0
	R0 = a(0) * abs(C(1, 0)) + a(1) * abs(C(0, 0));
	R1 = b(1) * abs(C(2, 2)) + b(2) * abs(C(2, 1));
	R = (C(0, 0) * A_matrix.row(1) * D - C(1, 0) * A_matrix.row(0) * D).norm();
	if (R > R0 + R1) return false;

	// L = A2 x B1
	R0 = a(0) * abs(C(1, 1)) + a(1) * abs(C(0, 1));
	R1 = b(0) * abs(C(2, 2)) + b(2) * abs(C(2, 0));
	R = (C(0, 1) * A_matrix.row(1) * D - C(1, 1) * A_matrix.row(0) * D).norm();
	if (R > R0 + R1) return false;

	// L = A2 x B2
	R0 = a(0) * abs(C(1, 2)) + a(1) * abs(C(0, 2));
	R1 = b(0) * abs(C(2, 1)) + b(1) * abs(C(2, 0));
	R = (C(0, 2) * A_matrix.row(1) * D - C(1, 2) * A_matrix.row(0) * D).norm();
	if (R > R0 + R1) return false;
	return true;
}

void SandBox::SetVelocity(double x, double y)
{
	xVelocity = x;
	yVelocity = y;
}

void  SandBox::IK_CCD()
{
	Eigen::Matrix4d destTrans = data_list[destIndex].MakeTransd();
	Eigen::Matrix4d firstTrans = data_list[firstLinkIndex].MakeTransd();
	Eigen::Matrix4d lastTrans = data_list[lastLinkIndex].MakeTransd();
	Eigen::Matrix4d lastParentsTrans = CalcParentsTrans(lastLinkIndex);

	Eigen::Vector4d lowerTip(0, -0.8, 0, 1);
	Eigen::Vector4d upperTip(0, 0.8, 0, 1);

	Eigen::Vector3d dst = ExtractPosition(destTrans * Eigen::Vector4d(0, 0, 0, 1));
	Eigen::Vector3d src = ExtractPosition(firstTrans * lowerTip);
	Eigen::Vector3d tip = ExtractPosition(lastParentsTrans * lastTrans * upperTip).transpose();

	double distToDest = (dst - tip).norm();
	double absDistance = (dst - src).norm();
	double maxReach = linksCounter * 1.6;
	double threshold = 0.1;

	std::cout << distToDest << std::endl;

	if (maxReach < absDistance) 
	{
		std::cout << "Cant reach" << std::endl;
		isActive = false;
		return;
	}

	if (distToDest < threshold)
	{
		std::cout << "Reached" << std::endl;
		isActive = false;
		return;
	}

	for (int i = firstLinkIndex + lastLinkIndex - 1; i >= firstLinkIndex; i--) 
	{
		Eigen::Matrix4d iParentsTrans = CalcParentsTrans(i);
		Eigen::Matrix4d iTrans = data_list[i].MakeTransd();
		Eigen::Vector3d D = dst;
		Eigen::Vector3d R = ExtractPosition(iParentsTrans * iTrans * lowerTip);
		Eigen::Vector3d E = ExtractPosition(lastParentsTrans * lastTrans * upperTip);
		Eigen::Vector3d RD = D - R;
		Eigen::Vector3d RE = E - R;
		float RD_dot_RE = RD.normalized().dot(RE.normalized());
		if (abs(RD_dot_RE) > 1)
		{
			std::cout << "Dot product greater then |1|, rounding value" << std::endl;
			RD_dot_RE = RD_dot_RE > 0 ? 1 : -1;
		}
		double theta = acosf(RD_dot_RE);
		Eigen::Matrix3d m = iParentsTrans.inverse().block(0, 0, 3, 3);
		Eigen::Vector3d rotAxis = m * RD.cross(RE);
		theta = theta / 10; // Smooth the rotation
		data_list[i].MyRotate(rotAxis, -theta);
	}
}

void  SandBox::IK_FABRIK()
{
	int n = linksCounter + 1;
	float di = 1.6;

	Eigen::Matrix4d destTrans = data_list[destIndex].MakeTransd();
	Eigen::Matrix4d firstTrans = data_list[firstLinkIndex].MakeTransd();
	Eigen::Matrix4d lastTrans = data_list[lastLinkIndex].MakeTransd();
	Eigen::Matrix4d lastParentsTrans = CalcParentsTrans(lastLinkIndex);

	Eigen::Vector4d lowerTip(0, -0.8, 0, 1);
	Eigen::Vector4d upperTip(0, 0.8, 0, 1);

	Eigen::Vector3d t = ExtractPosition(destTrans * Eigen::Vector4d(0, 0, 0, 1));
	Eigen::Vector3d p1 = ExtractPosition(firstTrans * lowerTip);
	Eigen::Vector3d pn = ExtractPosition(lastParentsTrans * lastTrans * upperTip).transpose();

	float dist = (p1 - t).norm();

	if (dist > di * (n - 1))
	{
		std::cout << "The target is unreachable" << std::endl;
		isActive = false;
		return;
	}

	float difA = (pn - t).norm();
	float tol = 0.1;

	std::cout << difA << std::endl;

	if (difA <= tol)
	{
		std::cout << "Reached target" << std::endl;
		isActive = false;
		return;
	}

	std::vector<Eigen::Vector3d> p;
	p.push_back(Eigen::Vector3d::Identity()); // Dummy
	p.push_back(p1);
	for (int i = 1; i <= linksCounter; i++)
	{
		p.push_back(ExtractPosition(data_list.at(i).MakeTransd() * CalcParentsTrans(i) * upperTip));
	}

	Eigen::Vector3d b = p.at(1);

	// FORWARD REACHING
	p.at(n) = t;

	for (int i = n - 1; i >= 1; i--)
	{
		Eigen::Vector3d pi = p.at(i);
		Eigen::Vector3d pi1 = p.at(i + 1);
		float ri = (pi1 - pi).norm();
		float li = di / ri;
		p.at(i) = ((1 - li) * pi1) + (li * pi);
	}

	// BACKWARD REACHING
	p.at(1) = b;

	for (int i = 1; i <= n - 1; i++)
	{
		Eigen::Vector3d pi = p.at(i);
		Eigen::Vector3d pi1 = p.at(i + 1);
		float ri = (pi1 - pi).norm();
		float li = di / ri;
		p.at(i + 1) = ((1 - li) * pi) + (li * pi1);
	}

	// UPDATE JOINTS POSITIONS
	for (int i = 1; i <= linksCounter; i++)
	{
		Eigen::Matrix4d iParentsTrans = CalcParentsTrans(i);
		Eigen::Matrix4d iTrans = data_list.at(i).MakeTransd();
		Eigen::Vector3d D = p.at(i + 1);
		Eigen::Vector3d R = ExtractPosition(iParentsTrans * iTrans * lowerTip);
		Eigen::Vector3d E = ExtractPosition(iParentsTrans * iTrans * upperTip);
		Eigen::Vector3d RD = D - R;
		Eigen::Vector3d RE = E - R;

		Eigen::Vector3d axis = RD.cross(RE).normalized();
		float RD_dot_RE = RD.normalized().dot(RE.normalized());
		if (abs(RD_dot_RE) > 1)
		{
			std::cout << "Dot product greater then |1|, rounding value" << std::endl;
			RD_dot_RE = RD_dot_RE > 0 ? 1 : -1;
		}
		float theta = acosf(RD_dot_RE);
		Eigen::Matrix4d m = iParentsTrans.inverse();
		Eigen::Matrix3d m1;
		m1 <<
			m(0, 0), m(0, 1), m(0, 2),
			m(1, 0), m(1, 1), m(1, 2),
			m(2, 0), m(2, 1), m(2, 2);

		Eigen::Vector3d rotAxis = m1 * RD.cross(RE);
		theta = theta / 10; // Smooth the rotation
		data_list.at(i).MyRotate(rotAxis, -theta);
	}
}

Eigen::Vector3d SandBox::ExtractPosition(Eigen::Vector4d m)
{
	return Eigen::Vector3d(m(0), m(1), m(2));
}
	