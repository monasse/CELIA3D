#include "solide.hpp"
#include "intersections.hpp"
#ifndef SOLIDE_CPP
#define SOLIDE_CPP

//const double eps_relat = numeric_limits<double>::epsilon();
const double eps_relat =0.000001;

Vertex::Vertex()
{
  pos = Point_3(0.,0.,0.);
  num = 0;
}


Vertex::Vertex(const Point_3 p, std::vector<int> & parts)
{
  pos = p;
  for(int i=0; i<parts.size(); i++){
	particules.push_back(parts[i]);
  }
}

Face::Face()
{
  centre = Point_3(0.,0.,0.);
  normale = Vector_3(1.,0.,0.);
  voisin = -1;
  D0 = 1.;
}


Face::Face(std::vector<Vertex> & v, int part)
{
  std::vector<Point_3> points;
  for(int i=0; i<v.size(); i++){
	vertex.push_back(v[i]);
	points.push_back(v[i].pos);
  }
  centre = centroid(points.begin(),points.end());
  normale = orthogonal_vector(points[0],points[1],points[2]);
  double norm = sqrt(CGAL::to_double(normale.squared_length()));
  normale = normale*1./norm;
  voisin = part;
  D0 = 1.;
}

Face::Face(std::vector<Vertex> & v, int part, double dist)
{
  std::vector<Point_3> points;
  for(int i=0; i<v.size(); i++){
	vertex.push_back(v[i]);
	points.push_back(v[i].pos);
  }
  centre = centroid(points.begin(),points.end());
  normale = orthogonal_vector(points[0],points[1],points[2]);
  double norm = sqrt(CGAL::to_double(normale.squared_length()));
  normale = normale*1./norm;
  voisin = part;
  D0 = dist;
}


Particule::Particule()
{   
	min_x = 0.; 
	min_y = 0.;
	min_z = 0.;
	max_x = 1. ;
	max_y = 1.;
	max_z = 1.;

	x0 = Point_3(0.5,0.5,0.5);
	
	cube = true;
	
	const Point_3 s1(min_x,min_y,min_z);
	const Point_3 r1(max_x, min_y, min_z);
	const Point_3 t1(max_x, max_y, min_z);
	const Point_3 v1(min_x, max_y, min_z);
	
	
	const Point_3 s2(min_x,min_y,max_z);
	const Point_3 r2(max_x, min_y, max_z);
	const Point_3 t2(max_x, max_y, max_z);
	const Point_3 v2(min_x, max_y, max_z);

	//Face 1
	std::vector<Vertex> vert1(4);
	vert1[0].pos = s1;
	vert1[0].particules.push_back(-1);
	vert1[1].pos = v1;
	vert1[1].particules.push_back(-1);
	vert1[2].pos = t1;
	vert1[2].particules.push_back(-1);
	vert1[3].pos = r1;
	vert1[3].particules.push_back(-1);
	vector<Point_3> points1(4);
	for(int i=0;i<4;i++){
	  points1[i] = vert1[i].pos;
	}
	Face face1 = Face(vert1,-1,1.);
	face1.centre = centroid(points1.begin(),points1.end());
	face1.normale = orthogonal_vector(points1[0],points1[1],points1[2]);
	double norm1 = sqrt(CGAL::to_double(face1.normale.squared_length()));
	face1.normale = face1.normale*1./norm1;

	//Face 2
	std::vector<Vertex> vert2(4);
	vert2[0].pos = s1;
	vert2[0].particules.push_back(-1);
	vert2[1].pos = r1;
	vert2[1].particules.push_back(-1);
	vert2[2].pos = r2;
	vert2[2].particules.push_back(-1);
	vert2[3].pos = s2;
	vert2[3].particules.push_back(-1);
	vector<Point_3> points2(4);
	for(int i=0;i<4;i++){
	  points2[i] = vert2[i].pos;
	}
	Face face2 = Face(vert2,-1,1.);
	face2.centre = centroid(points2.begin(),points2.end());
	face2.normale = orthogonal_vector(points2[0],points2[1],points2[2]);
	double norm2 = sqrt(CGAL::to_double(face2.normale.squared_length()));
	face2.normale = face2.normale*1./norm2;

	//Face 3
	std::vector<Vertex> vert3(4);
	vert3[0].pos = r1;
	vert3[0].particules.push_back(-1);
	vert3[1].pos = t1;
	vert3[1].particules.push_back(-1);
	vert3[2].pos = t2;
	vert3[2].particules.push_back(-1);
	vert3[3].pos = r2;
	vert3[3].particules.push_back(-1);
	vector<Point_3> points3(4);
	for(int i=0;i<4;i++){
	  points3[i] = vert3[i].pos;
	}
	Face face3 = Face(vert3,-1,1.);
	face3.centre = centroid(points3.begin(),points3.end());
	face3.normale = orthogonal_vector(points3[0],points3[1],points3[2]);
	double norm3 = sqrt(CGAL::to_double(face3.normale.squared_length()));
	face3.normale = face3.normale*1./norm3;

	//Face 4
	std::vector<Vertex> vert4(4);
	vert4[0].pos = t1;
	vert4[0].particules.push_back(-1);
	vert4[1].pos = v1;
	vert4[1].particules.push_back(-1);
	vert4[2].pos = v2;
	vert4[2].particules.push_back(-1);
	vert4[3].pos = t2;
	vert4[3].particules.push_back(-1);
	vector<Point_3> points4(4);
	for(int i=0;i<4;i++){
	  points4[i] = vert4[i].pos;
	}
	Face face4 = Face(vert4,-1,1.);
	face4.centre = centroid(points4.begin(),points4.end());
	face4.normale = orthogonal_vector(points4[0],points4[1],points4[2]);
	double norm4 = sqrt(CGAL::to_double(face4.normale.squared_length()));
	face4.normale = face4.normale*1./norm4;
	
	//Face 5
	std::vector<Vertex> vert5(4);
	vert5[0].pos = s1;
	vert5[0].particules.push_back(-1);
	vert5[1].pos = r1;
	vert5[1].particules.push_back(-1);
	vert5[2].pos = r2;
	vert5[2].particules.push_back(-1);
	vert5[3].pos = s2;
	vert5[3].particules.push_back(-1);
	vector<Point_3> points5(4);
	for(int i=0;i<4;i++){
	  points5[i] = vert5[i].pos;
	}
	Face face5 = Face(vert5,-1,1.);
	face5.centre = centroid(points5.begin(),points5.end());
	face5.normale = orthogonal_vector(points5[0],points5[1],points5[2]);
	double norm5 = sqrt(CGAL::to_double(face5.normale.squared_length()));
	face5.normale = face5.normale*1./norm5;
  
	//Face 6
	std::vector<Vertex> vert6(4);
	vert6[0].pos = s1;
	vert6[0].particules.push_back(-1);
	vert6[1].pos = r1;
	vert6[1].particules.push_back(-1);
	vert6[2].pos = r2;
	vert6[2].particules.push_back(-1);
	vert6[3].pos = s2;
	vert6[3].particules.push_back(-1);
	vector<Point_3> points6(4);
	for(int i=0;i<4;i++){
	  points6[i] = vert6[i].pos;
	}
	Face face6 = Face(vert6,-1,1.);
	face6.centre = centroid(points6.begin(),points6.end());
	face6.normale = orthogonal_vector(points6[0],points6[1],points6[2]);
	double norm6 = sqrt(CGAL::to_double(face6.normale.squared_length()));
	face6.normale = face6.normale*1./norm6;

	std::vector<Face> f(6);
	f[0] = face1;
	f[1] = face2;
	f[2] = face3;
	f[3] = face4;
	f[4] = face5;
	f[5] = face6;

	faces = f;
	
	//triangles par face	
	Triangle_3 Tri1(s1,r1,v1);
	Triangle_3 Tri2(t1,r1,v1);	
	triangles.push_back(Tri1);
	triangles.push_back(Tri2);
	normales.push_back(face1.normale);
	normales.push_back(face1.normale);
	fluide.push_back(true);
	fluide.push_back(true);
	//face2
	Triangle_3 Tri5(s2,r2,v2);
	Triangle_3 Tri6(t2,r2,v2);
	triangles.push_back(Tri5);
	triangles.push_back(Tri6);
	normales.push_back(face2.normale);
	normales.push_back(face2.normale);
	fluide.push_back(true);
	fluide.push_back(true);
	//face3
	Triangle_3 Tri9(s2,s1,v2);
	Triangle_3 Tri10(v1,s1,v2);
	triangles.push_back(Tri9);
	triangles.push_back(Tri10);
	normales.push_back(face3.normale);
	normales.push_back(face3.normale);
	fluide.push_back(true);
	fluide.push_back(true);
	//face4	
	Triangle_3 Tri13(r2,r1,t2);
	Triangle_3 Tri14(t1,r1,t2);
	triangles.push_back(Tri13);
	triangles.push_back(Tri14);
	normales.push_back(face4.normale);
	normales.push_back(face4.normale);
	fluide.push_back(true);
	fluide.push_back(true);
	//face5	
	Triangle_3 Tri17(v2,v1,t2);
	Triangle_3 Tri18(t1,v1,t2);
	triangles.push_back(Tri17);
	triangles.push_back(Tri18);
	normales.push_back(face5.normale);
	normales.push_back(face5.normale);
	fluide.push_back(true);
	fluide.push_back(true);
	//face6
	Triangle_3 Tri21(s2,s1,r2);
	Triangle_3 Tri22(r1,s1,r2);
	triangles.push_back(Tri21);
	triangles.push_back(Tri22); 
	normales.push_back(face6.normale);
	normales.push_back(face6.normale);
	fluide.push_back(true);
	fluide.push_back(true);
	Points_interface.resize(triangles.size(), std::vector<Point_3>(0));
	Triangles_interface.resize(triangles.size(), std::vector<Triangle_3>(0));
}

Particule::Particule(const double x_min, const double y_min, const double z_min, 
							 const double x_max, const double y_max,const double z_max)
{   
	min_x = x_min; 
	min_y = y_min;
	min_z = z_min;
	max_x = x_max ;
	max_y = y_max;
	max_z = z_max;

	x0 = Point_3((x_min+x_max)/2.,(y_min+y_max)/2.,(z_min+z_max)/2.);
	
	cube = true;
	
	const Point_3 s1(min_x,min_y,min_z);
	const Point_3 r1(max_x, min_y, min_z);
	const Point_3 t1(max_x, max_y, min_z);
	const Point_3 v1(min_x, max_y, min_z);
	
	
	const Point_3 s2(min_x,min_y,max_z);
	const Point_3 r2(max_x, min_y, max_z);
	const Point_3 t2(max_x, max_y, max_z);
	const Point_3 v2(min_x, max_y, max_z);
	
	//Face 1
	std::vector<Vertex> vert1(4);
	vert1[0].pos = s1;
	vert1[0].particules.push_back(-1);
	vert1[1].pos = v1;
	vert1[1].particules.push_back(-1);
	vert1[2].pos = t1;
	vert1[2].particules.push_back(-1);
	vert1[3].pos = r1;
	vert1[3].particules.push_back(-1);
	vector<Point_3> points1(4);
	for(int i=0;i<4;i++){
	  points1[i] = vert1[i].pos;
	}
	Face face1 = Face(vert1,-1,1.);
	face1.centre = centroid(points1.begin(),points1.end());
	face1.normale = orthogonal_vector(points1[0],points1[1],points1[2]);
	double norm1 = sqrt(CGAL::to_double(face1.normale.squared_length()));
	face1.normale = face1.normale*1./norm1;
	
	//Face 2
	std::vector<Vertex> vert2(4);
	vert2[0].pos = s1;
	vert2[0].particules.push_back(-1);
	vert2[1].pos = r1;
	vert2[1].particules.push_back(-1);
	vert2[2].pos = r2;
	vert2[2].particules.push_back(-1);
	vert2[3].pos = s2;
	vert2[3].particules.push_back(-1);
	vector<Point_3> points2(4);
	for(int i=0;i<4;i++){
	  points2[i] = vert2[i].pos;
	}
	Face face2 = Face(vert2,-1,1.);
	face2.centre = centroid(points2.begin(),points2.end());
	face2.normale = orthogonal_vector(points2[0],points2[1],points2[2]);
	double norm2 = sqrt(CGAL::to_double(face2.normale.squared_length()));
	face2.normale = face2.normale*1./norm2;
	
	//Face 3
	std::vector<Vertex> vert3(4);
	vert3[0].pos = r1;
	vert3[0].particules.push_back(-1);
	vert3[1].pos = t1;
	vert3[1].particules.push_back(-1);
	vert3[2].pos = t2;
	vert3[2].particules.push_back(-1);
	vert3[3].pos = r2;
	vert3[3].particules.push_back(-1);
	vector<Point_3> points3(4);
	for(int i=0;i<4;i++){
	  points3[i] = vert3[i].pos;
	}
	Face face3 = Face(vert3,-1,1.);
	face3.centre = centroid(points3.begin(),points3.end());
	face3.normale = orthogonal_vector(points3[0],points3[1],points3[2]);
	double norm3 = sqrt(CGAL::to_double(face3.normale.squared_length()));
	face3.normale = face3.normale*1./norm3;

	//Face 4
	std::vector<Vertex> vert4(4);
	vert4[0].pos = t1;
	vert4[0].particules.push_back(-1);
	vert4[1].pos = v1;
	vert4[1].particules.push_back(-1);
	vert4[2].pos = v2;
	vert4[2].particules.push_back(-1);
	vert4[3].pos = t2;
	vert4[3].particules.push_back(-1);
	vector<Point_3> points4(4);
	for(int i=0;i<4;i++){
	  points4[i] = vert4[i].pos;
	}
	Face face4 = Face(vert4,-1,1.);
	face4.centre = centroid(points4.begin(),points4.end());
	face4.normale = orthogonal_vector(points4[0],points4[1],points4[2]);
	double norm4 = sqrt(CGAL::to_double(face4.normale.squared_length()));
	face4.normale = face4.normale*1./norm4;

	//Face 5
	std::vector<Vertex> vert5(4);
	vert5[0].pos = s1;
	vert5[0].particules.push_back(-1);
	vert5[1].pos = r1;
	vert5[1].particules.push_back(-1);
	vert5[2].pos = r2;
	vert5[2].particules.push_back(-1);
	vert5[3].pos = s2;
	vert5[3].particules.push_back(-1);
	vector<Point_3> points5(4);
	for(int i=0;i<4;i++){
	  points5[i] = vert5[i].pos;
	}
	Face face5 = Face(vert5,-1,1.);
	face5.centre = centroid(points5.begin(),points5.end());
	face5.normale = orthogonal_vector(points5[0],points5[1],points5[2]);
	double norm5 = sqrt(CGAL::to_double(face5.normale.squared_length()));
	face5.normale = face5.normale*1./norm5;

	//Face 6
	std::vector<Vertex> vert6(4);
	vert6[0].pos = s1;
	vert6[0].particules.push_back(-1);
	vert6[1].pos = r1;
	vert6[1].particules.push_back(-1);
	vert6[2].pos = r2;
	vert6[2].particules.push_back(-1);
	vert6[3].pos = s2;
	vert6[3].particules.push_back(-1);
	vector<Point_3> points6(4);
	for(int i=0;i<4;i++){
	  points6[i] = vert6[i].pos;
	}
	Face face6 = Face(vert6,-1,1.);
	face6.centre = centroid(points6.begin(),points6.end());
	face6.normale = orthogonal_vector(points6[0],points6[1],points6[2]);
	double norm6 = sqrt(CGAL::to_double(face6.normale.squared_length()));
	face6.normale = face6.normale*1./norm6;

	std::vector<Face> f(6);
	f[0] = face1;
	f[1] = face2;
	f[2] = face3;
	f[3] = face4;
	f[4] = face5;
	f[5] = face6;

	faces = f;
	
	//triangles par face	
	Triangle_3 Tri1(s1,r1,v1);
	Triangle_3 Tri2(t1,r1,v1);	
	triangles.push_back(Tri1);
	triangles.push_back(Tri2);
	normales.push_back(face1.normale);
	normales.push_back(face1.normale);
	fluide.push_back(true);
	fluide.push_back(true);
	//face2
	Triangle_3 Tri5(s2,r2,v2);
	Triangle_3 Tri6(t2,r2,v2);
	triangles.push_back(Tri5);
	triangles.push_back(Tri6);
	normales.push_back(face2.normale);
	normales.push_back(face2.normale);
	fluide.push_back(true);
	fluide.push_back(true);
	//face3
	Triangle_3 Tri9(s2,s1,v2);
	Triangle_3 Tri10(v1,s1,v2);
	triangles.push_back(Tri9);
	triangles.push_back(Tri10);
	normales.push_back(face3.normale);
	normales.push_back(face3.normale);
	fluide.push_back(true);
	fluide.push_back(true);
	//face4	
	Triangle_3 Tri13(r2,r1,t2);
	Triangle_3 Tri14(t1,r1,t2);
	triangles.push_back(Tri13);
	triangles.push_back(Tri14);
	normales.push_back(face4.normale);
	normales.push_back(face4.normale);
	fluide.push_back(true);
	fluide.push_back(true);
	//face5	
	Triangle_3 Tri17(v2,v1,t2);
	Triangle_3 Tri18(t1,v1,t2);
	triangles.push_back(Tri17);
	triangles.push_back(Tri18);
	normales.push_back(face5.normale);
	normales.push_back(face5.normale);
	fluide.push_back(true);
	fluide.push_back(true);
	//face6
	Triangle_3 Tri21(s2,s1,r2);
	Triangle_3 Tri22(r1,s1,r2);
	triangles.push_back(Tri21);
	triangles.push_back(Tri22); 
	normales.push_back(face6.normale);
	normales.push_back(face6.normale);
	fluide.push_back(true);
	fluide.push_back(true);
	Points_interface.resize(triangles.size(), std::vector<Point_3>(0));
	Triangles_interface.resize(triangles.size(), std::vector<Triangle_3>(0));
}


Particule::Particule(Point_3 c, const double x_min, const double y_min, const double z_min, 
							 const double x_max, const double y_max,const double z_max, 
					 std::vector<Face> & F)
{   
	min_x = x_min; 
	min_y = y_min;
	min_z = z_min;
	max_x = x_max ;
	max_y = y_max;
	max_z = z_max;

	x0 = c;
	
	cube = false;

	faces = F;

	for(int i=0;i<faces.size();i++){
	  Point_3 s,r,v;
	  s = faces[i].vertex[0].pos;
	  for(int k=1;k<faces[i].size()-1;k++){
		r = faces[i].vertex[k].pos;
		v = faces[i].vertex[k+1].pos;
		Vector_3 vect0(s,r);
		Vector_3 vect1(s,v);
		//Verification que les faces ne sont pas alignees
		for(int j=k+2;(j<faces[i].size()) && (CGAL::to_double(vect0*vect1) == 0.);j++){
		  v = faces[i].vertex[j].pos;
		  vect1 = Vector_3(s,v);
		  k++;
		}
		Triangle_3 Tri(s,r,v);
		triangles.push_back(Tri);
		normales.push_back(faces[i].normale);
		if(faces[i].voisin == -1){
		  fluide.push_back(true);
		} else {
		  fluide.push_back(false);
		}
	  }
	}
	Points_interface.resize(triangles.size(), std::vector<Point_3>(0));
	Triangles_interface.resize(triangles.size(), std::vector<Triangle_3>(0));
}
//Destructeur
Particule::~Particule(){
	
}



void Particule::Affiche(){
	
//	std::cout<<" volume of solide := "<<volume()<<std::endl;
// 	std::cout<<" Point min := "<< min_x<<std::endl;
// 	std::cout<<" Point max := "<< max_x<<std::endl;

// std::cout<<" triangles size := "<< triangles.size()<<std::endl;
//std::cout<<" P interface size := "<< Points_interface.size()<<std::endl;
// std::cout<<" T interface size := "<< Triangles_interface.size()<<std::endl;

//  for(int i=0; i<Triangles_interface.size(); i++){
// 	if(!fluide[i])
// 	//std::cout<<"nr triangle "<<j<<" contact fluide "<<fluide[j]<<"  points interface:= "<<interface[j][i]<<std::endl;
// 	std::cout<<"triangles size"<<Triangles_interface[i].size()<<std::endl;
// }

}




double Particule::volume(){
	
	double vol = 0.;
	std::vector<Point_3> Points_poly; 
	
	for(int l= 0; l<triangles.size(); l++)
	{
		Points_poly.push_back(triangles[l].operator[](0));
		Points_poly.push_back(triangles[l].operator[](1));
		Points_poly.push_back(triangles[l].operator[](2));
	}	
	Finite_cells_iterator cit;
	Triangulation T(Points_poly.begin(), Points_poly.end());
	
	for (cit = T.finite_cells_begin(); cit!= T.finite_cells_end(); cit++){
		vol+= CGAL::to_double(T.tetrahedron( cit).volume());
	}
	
	return vol;
}


void Solide::impression(int n){ //Sortie au format vtk
int nb_part = solide.size();

int nb_triangles = 0.;
for(int it=0; it<nb_part; it++){
	nb_triangles += solide[it].triangles.size();
}

const char* solidevtk;
{
	std::ostringstream oss;
	oss << "resultats/solide" << n << ".vtk";
	string s = oss.str();
	//cout << s << endl;
	solidevtk = s.c_str();
}

//Ouverture des flux en donne en ecriture
std::ofstream vtk(solidevtk,ios::out);
if(vtk)
{
	// cout <<"ouverture de xt.vtk reussie" << endl;
} else {
	cout <<"ouverture de solide" << n << ".vtk rate" << endl;
}
//Initialisation du fichier vtk
vtk << "# vtk DataFile Version 3.0" << endl;
vtk << "#Simulation Euler" << endl;
vtk << "ASCII" << endl;
vtk<<"\n";
vtk << "DATASET UNSTRUCTURED_GRID" << endl;
vtk << "POINTS " << 3*nb_triangles << " DOUBLE" << endl;

for(int it=0; it<nb_part; it++){
	for(int l= 0; l<solide[it].triangles.size(); l++){
		vtk << solide[it].triangles[l].operator[](0).operator[](0) << " " << solide[it].triangles[l].operator[](0).operator[](1) << " " << solide[it].triangles[l].operator[](0).operator[](2) << endl;
		vtk << solide[it].triangles[l].operator[](1).operator[](0) << " " << solide[it].triangles[l].operator[](1).operator[](1) << " " << solide[it].triangles[l].operator[](1).operator[](2) << endl;
		vtk << solide[it].triangles[l].operator[](2).operator[](0) << " " << solide[it].triangles[l].operator[](2).operator[](1) << " " << solide[it].triangles[l].operator[](2).operator[](2) << endl;
	}
}
vtk<<"\n";
vtk << "CELLS " << nb_triangles << " " << 4*nb_triangles<< endl;
int num=0;
for(int it=0; it<nb_part; it++){
	for(int l= 0; l<solide[it].triangles.size(); l++){
		vtk << 3 << " " << 3*num << " " << 3*num+1 << " " << 3*num+2 << endl;
		num++;
	}
}
vtk << "\n";
vtk << "CELL_TYPES " << nb_triangles << endl;
for(int l= 0; l<nb_triangles; l++)
{
	vtk << 5 << endl;
}
vtk << "\n";
vtk << "CELL_DATA " << nb_triangles << endl;
//Deplacement
vtk << "SCALARS pression double 1" << endl;
vtk << "LOOKUP_TABLE default" << endl;
for(int it=0; it<nb_part; it++){
	for(int l= 0; l<solide[it].triangles.size(); l++)
	{
		vtk << 0. << endl;
	}
}
}
	
bool inside_convex_polygon(const Particule& S, const Point_3& P){
	
	bool in = false;
	
	if((S.min_x - P.x())<= eps_relat && (S.min_y - P.y())<= eps_relat && 
		 (S.min_z - P.z())<= eps_relat && (S.max_x - P.x())>=-eps_relat && 
		 (S.max_y - P.y())>=-eps_relat && (S.max_z - P.z())>=-eps_relat )
	{
		if(S.cube) {in = true;}
		else{
			in = true;
			for(int l= 0; l<S.triangles.size(); l++){
				Point_3 vertex = S.triangles[l].operator[](0);
				Vector_3 vect(P,vertex);
				if((CGAL::to_double(vect*S.normales[l])) < 0.){in = false;}
			}
		}
	}
	
	return in;
}	

bool box_inside_convex_polygon(const Particule& S, const Bbox& cell){
	
	bool in = false;
	
	if ((S.min_x - cell.xmin()) <= eps_relat && (S.min_y - cell.ymin() <= eps_relat) && 
		  (S.min_z - cell.zmin()) <= eps_relat && (S.max_x - cell.xmax() >=-eps_relat) && 
		  (S.max_y - cell.ymax()) >=-eps_relat && (S.max_z - cell.zmax()>=-eps_relat) ) 
	{
		if(S.cube) { return S.cube;}
		
		else{
			
	  in = true;
		
		Point_3 p1(cell.xmin(),cell.ymin(),cell.zmin());
		in = inside_convex_polygon(S,p1);
		if(!in) {return in;}
		
		Point_3 p2(cell.xmax(),cell.ymax(),cell.zmax());
		in = inside_convex_polygon(S,p2);
		if(!in) {return in;}
		
		Point_3 p3(cell.xmax(),cell.ymin(),cell.zmin());
		in = inside_convex_polygon(S,p3);
		if(!in) {return in;}
		
		Point_3 p4(cell.xmax(),cell.ymin(),cell.zmax());
		in = inside_convex_polygon(S,p4);
		if(!in) {return in;}
		
		Point_3 p5(cell.xmax(),cell.ymax(),cell.zmin());
		in = inside_convex_polygon(S,p5);
		if(!in) {return in;}
		
		Point_3 p6(cell.xmin(),cell.ymax(),cell.zmin());
		in = inside_convex_polygon(S,p6);
		if(!in) {return in;}
		
		Point_3 p7(cell.xmin(),cell.ymax(),cell.zmax());
		in = inside_convex_polygon(S,p7);
		if(!in) {return in;}
		
		Point_3 p8(cell.xmin(),cell.ymin(),cell.zmax());
		in = inside_convex_polygon(S,p8);
		if(!in) {return in;}
		
		}
	}
	
	return in;
}


bool inside_box(const Bbox& cell, const Point_3& P){
	
	bool in = false;
	
	if((cell.xmin() - P.x())<= eps_relat && (cell.ymin() - P.y())<= eps_relat &&
		 (cell.zmin() - P.z())<= eps_relat && (cell.xmax() - P.x())>=-eps_relat &&
		 (cell.ymax() - P.y())>=-eps_relat && (cell.zmax() - P.z())>=-eps_relat )
	  { in = true; }
	
	return in;
}

Solide::Solide(){
}

Solide::Solide(std::vector<Particule> & Part){
  for(int i=0; i<Part.size(); i++){
	solide.push_back(Part[i]);
  }
	
}

Solide::~Solide(){   
}


void Solide::Affiche(){
	
	for(int i=0; i<solide.size(); i++){
		solide[i].Affiche();
	}

}

void Solide::init(const char* s){
  std::ifstream maillage(s,ios::in);
  if(maillage){
	// cout <<"ouverture de xt.vtk reussie" << endl;
  } else {
	cout <<"ouverture de " << s << " ratee" << endl;
  }

  //Recuperation du maillage solide
  int Npoint;
  string sp;
  maillage >> sp >> Npoint;
  const int nb_points = Npoint;
  
  vector<Point_3> Points(nb_points);
  
  for(int i=0;i<nb_points;i++){
	double x,y,z;
	maillage >> x >> y >> z;
	Points[i] = Point_3(x,y,z);
  }
  
  int Npart;
  string sP;
  maillage >> sP >> Npart;
  const int nb_particule = Npart;
  
  vector<Particule> P(nb_particule);
  
  bool points_particules[nb_points][nb_particule];
  
  for(int i=0;i<nb_particule;i++){
	for(int j=0;j<nb_points;j++){
	  //Remise a zero du tableau
	  points_particules[j][i] = false;
	}
  }
  for(int i=0;i<nb_particule;i++){
	int Nfaces;
	bool fixe;
	double X,Y,Z,u,v,w,theta,phi,psi,xmin,ymin,zmin,xmax,ymax,zmax;
	string s;
	maillage >> s >> Nfaces >> fixe;
	maillage >> s >> X >> Y >> Z;
	Point_3 centre(X,Y,Z);
	maillage >> s >> u >> v >> w;
	maillage >> s >> theta >> phi >> psi;
	xmin = X;
	ymin = Y;
	zmin = Z;
	const int nb_faces = Nfaces;
	std::vector<Face> Faces(nb_faces);
	for(int j=0;j<nb_faces;j++){
	  int Nvertex;
	  maillage >> Nvertex;
	  const int nb_vertex = Nvertex;
	  std::vector<Vertex> Vertex(nb_vertex);
	  for(int k=0;k<nb_vertex;k++){
		int p;
		maillage >> p;
		Vertex[k].pos = Points[p];
		Vertex[k].num = p;
		points_particules[p][i] = true;
		double x = CGAL::to_double(Points[p].operator[](0));
		double y = CGAL::to_double(Points[p].operator[](1));
		double z = CGAL::to_double(Points[p].operator[](2));
		xmin = min(x,xmin);
		xmax = max(x,xmax);
		ymin = min(y,ymin);
		ymax = max(y,ymax);
		zmin = min(z,zmin);
		zmax = max(z,zmax);
	  }
	  int voisin;
	  maillage >> voisin;
	  Faces[j] = Face::Face(Vertex, voisin);
	}
	P[i] = Particule::Particule(centre, xmin, ymin, zmin, xmax, ymax, zmax, Faces);
  }
  //Boucle de mise a jour des particules sur les sommets du maillage
  //Mise a jour des distances a l'equilibre entre particules en meme temps
  for(int i=0;i<P.size();i++){
	for(int j=0;j<P[i].size();j++){
	  for(int k=0;k<P[i].faces[j].size();k++){
		for(int l=0;l<P.size();l++){
		  if(points_particules[P[i].faces[j].vertex[k].num][l]){
			P[i].faces[j].vertex[k].particules.push_back(l);
		  }
		}
		if(P[i].faces[j].voisin>=0){
		  P[i].faces[j].D0 = sqrt(CGAL::to_double(CGAL::squared_distance(P[i].x0,P[P[i].faces[j].voisin].x0)));
		}
	  }
	}
  }

  for(int i=0; i<P.size(); i++){
	solide.push_back(P[i]);
  }
  
  //Initialisation de la position du solide
  for(int i=0; i<solide.size(); i++){
	solide[i].Dx = Point_3(0.,0.,0.);
	solide[i].Dxprev = Point_3(0.,0.,0.);
	solide[i].Fi = Vector_3(0.,0.,0.);
	solide[i].Ff = Vector_3(0.,0.,0.);
	solide[i].Mi = Vector_3(0.,0.,0.);
	solide[i].Mf = Vector_3(0.,0.,0.);
	solide[i].rot[0][0] = 1.;
	solide[i].rot[1][1] = 1.;
	solide[i].rot[2][2] = 1.;
	solide[i].rot[0][1] = 0.;
	solide[i].rot[0][2] = 0.;
	solide[i].rot[1][0] = 0.;
	solide[i].rot[1][2] = 0.;
	solide[i].rot[2][0] = 0.;
	solide[i].rot[2][1] = 0.;
  }
}


#endif
