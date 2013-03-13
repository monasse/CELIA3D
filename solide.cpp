#include "solide.hpp"
#include "intersections.hpp"
#include<iostream>
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

// Vertex & Vertex:: operator=(const Vertex &V){
// 	
// 	assert(this != &V);
// 	pos = V.pos;
// 	num = V.num;
// 	for(int i=0; i<V.particules.size(); i++){
// 	particules[i]= V.particules[i];
// 	}
// }

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

void Face::Inertie(){
  //Choix initial d'un repere orthonorme de la face
  if(normale.operator[](0)!=0. || normale.operator[](1)!=0.){
    s = Vector_3(-normale.operator[](1),normale.operator[](0),0.);
    s = s/(sqrt(CGAL::to_double(s.squared_length())));
    t = cross_product(normale,s);
  } else {
    s = Vector_3(0.,0.,1.);
    s = s/(sqrt(CGAL::to_double(s.squared_length())));
    t = cross_product(normale,s);
  }
  //Calcul du centre de la face
  double T1 = 0.;
  double Ts = 0.;
  double Tt = 0.;
  for(int i=0;i<size();i++){
    int ip = (i+1)%(size());
    Vector_3 v1(centre,vertex[i].pos);
    Vector_3 v2(centre,vertex[ip].pos);
    T1 += 1./2.*CGAL::to_double(cross_product(v1,v2)*normale);
    Ts += 1./6.*CGAL::to_double(cross_product(v1,v2)*normale)*CGAL::to_double((v1+v2)*s);
    Tt += 1./6.*CGAL::to_double(cross_product(v1,v2)*normale)*CGAL::to_double((v1+v2)*t);
  }
  centre = centre + (Ts/T1)*s + (Tt/T1)*t;
  //Calcul de la matrice d'inertie de la face dans les deux axes a l'origine centre
  double Tss = 0.;
  double Ttt = 0.;
  double Tst = 0.;
  for(int i=0;i<size();i++){
    int ip = (i+1)%(size());
    Vector_3 v1(centre,vertex[i].pos);
    Vector_3 v2(centre,vertex[ip].pos);
    double As = CGAL::to_double(v1*s);
    double At = CGAL::to_double(v1*t);
    double Bs = CGAL::to_double(v2*s);
    double Bt = CGAL::to_double(v2*t);
    Tss += 1./12.*(As*As+As*Bs+Bs*Bs);
    Ttt += 1./12.*(At*At+At*Bt+Bt*Bt);
    Tst += 1./24.*(2.*As*At+As*Bt+At*Bs+2.*Bs*Bt);
  }
  //Calcul des moments d'inertie
  double Delta = pow(Tss-Ttt,2)+4.*Tst*Tst;
  Is = (Tss+Ttt+sqrt(Delta))/2.;
  It = (Tss+Ttt-sqrt(Delta))/2.;
  //Diagonalisation
  if(abs(Tss-Ttt)>eps){
    if(abs(Tss-Is)>eps){
      Vector_3 stemp = -Tst*s+(Tss-Is)*t;
      s = stemp/(sqrt(CGAL::to_double(stemp.squared_length())));
      t = cross_product(normale,s);
    } else {
      Vector_3 stemp = -Tst*t+(Ttt-Is)*s;
      s = stemp/(sqrt(CGAL::to_double(stemp.squared_length())));
      t = cross_product(normale,s);
    }
  } else {
    if(abs(Tst)>eps){
      Vector_3 stemp = s+t;
      Vector_3 ttemp = -s+t;
      s = stemp/(sqrt(CGAL::to_double(stemp.squared_length())));
      t = stemp/(sqrt(CGAL::to_double(ttemp.squared_length())));
    }
  }
}


// Face & Face:: operator=(const Face &F){
// 	
// 	assert(this != &F);
//   centre = F.centre;
//   normale = F.normale;
//   voisin = F.voisin;
//   D0  = F.D0; 
// 	
// 	for(int i= 0; i<F.vertex.size(); i++){
// 	vertex[i] = F.vertex[i];
// 	}
// }

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
	Position_Triangles_interface.resize(triangles.size(), std::vector< std::vector<int> >(0));
	Points_interface_prev.resize(triangles.size(), std::vector<Point_3>(0));
	Triangles_interface_prev.resize(triangles.size(), std::vector<Triangle_3>(0));
	Position_Triangles_interface_prev.resize(triangles.size(), std::vector<std::vector<int> >(0));
	Ff = Vector_3(0.,0.,0.); Ffprev = Vector_3(0.,0.,0.); 
	Mf = Vector_3(0.,0.,0.); Mfprev = Vector_3(0.,0.,0.);
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
	Position_Triangles_interface.resize(triangles.size(), std::vector< std::vector<int> >(0));
	Points_interface_prev.resize(triangles.size(), std::vector<Point_3>(0));
	Triangles_interface_prev.resize(triangles.size(), std::vector<Triangle_3>(0));
	Position_Triangles_interface_prev.resize(triangles.size(), std::vector< std::vector<int> >(0));
	Ff = Vector_3(0.,0.,0.); Ffprev = Vector_3(0.,0.,0.); 
	Mf = Vector_3(0.,0.,0.); Mfprev = Vector_3(0.,0.,0.);
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
    s = faces[i].centre;
    for(int k=0;k<faces[i].size();k++){
      int kp = (k+1)%(faces[i].size());
      r = faces[i].vertex[k].pos;
      v = faces[i].vertex[kp].pos;
      Vector_3 vect0(s,r);
      Vector_3 vect1(s,v);
      /*Verification que les faces ne sont pas alignees
      for(int j=kp+1;(j<faces[i].size()) && (CGAL::to_double(vect0*vect1) == 0.);j++){
	v = faces[i].vertex[j].pos;
	vect1 = Vector_3(s,v);
	k++;
	}*/
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
	Position_Triangles_interface.resize(triangles.size(), std::vector< std::vector<int> >(0));
	Points_interface_prev.resize(triangles.size(), std::vector<Point_3>(0));
	Triangles_interface_prev.resize(triangles.size(), std::vector<Triangle_3>(0));
	Position_Triangles_interface_prev.resize(triangles.size(), std::vector< std::vector<int> >(0));
	Ff = Vector_3(0.,0.,0.); Ffprev = Vector_3(0.,0.,0.); 
	Mf = Vector_3(0.,0.,0.); Mfprev = Vector_3(0.,0.,0.);
}
//Destructeur
Particule::~Particule(){
	
}

// Particule & Particule:: operator=(const Particule &P){
// 	
// 	assert(this != &P);
// // 	vector< vector<int> >::iterator iter_ii;
// // 	vector<int>::iterator                 iter_jj;
// 	
// 	min_x = P.min_x;
// 	min_y = P.min_y;
// 	min_z = P.min_z;
// 	max_x = P.max_x;
// 	max_y = P.max_y;
// 	max_z = P.max_z;
// 	cube  = P.cube;
// 	
// 	//faces.resize(P.faces.size());
// 	faces = P.faces;
// // 	triangles.resize(P.triangles);
// // 	normales.resize(P.normales);
// // 	fluide.resize(P.fluide);
// // 	assert(triangles.size()==normales.size());
// 	for(int i=0; i<triangles.size(); i++){
// 		triangles[i] = P.triangles[i];
// 		normales[i] = P.normales[i];
// 		fluide[i] = P.fluide[i];
// 	}
// 	
// 	x0 = P.x0; Dx = P.Dx; Dxprev = P.Dxprev; Fi = P.Fi; Ff=P.Ff; Mi= P.Mi; Mf=P.Mf;
// 	
// 	for(int i= 0; i<3;i++){
// 		for(int j= 0; j<3;j++){
// 			rot[i][j] = P.rot[i][j];
// 		}
// 	}
// 	
// // 	Points_interface.resize(P.Points_interface.size(), std::vector<Point_3>(0));
// // 	Triangles_interface.resize(P.Triangles_interface.size(), std::vector<Triangle_3>(0));
// 	
// 	for(int i=0; i<Points_interface.size(); i++ ){
// 		for(int j=0; j<Points_interface[i].size(); j++ ){
// 		Points_interface[i][j]= P.Points_interface[i][j];
// 		}
// 	}
// 	
// 	for(int i=0; i<Triangles_interface.size(); i++ ){
// 		for(int j=0; j<Triangles_interface[i].size(); j++ ){
// 		Triangles_interface[i][j]= P.Triangles_interface[i][j];
// 		}
// 	}
// 
// }
void Particule::Affiche(){
	
//	std::cout<<" volume of solide := "<<volume()<<std::endl;
 	std::cout<<" Point min x:= "<< min_x<<std::endl;
  std::cout<<" Point max x:= "<< max_x<<std::endl;
	std::cout<<" Point min y:= "<< min_y<<std::endl;
	std::cout<<" Point max y:= "<< max_y<<std::endl;
	std::cout<<" Point min z:= "<< min_z<<std::endl;
	std::cout<<" Point max z:= "<< max_z<<std::endl;
// 	for(int i=0; i<triangles.size(); i++){
// 		std::cout<<" triangles:= "<< triangles[i]<<std::endl;
// 	}
//std::cout<<" P interface size := "<< Points_interface.size()<<std::endl;
// std::cout<<" T interface size := "<< Triangles_interface.size()<<std::endl;

//  for(int i=0; i<Triangles_interface.size(); i++){
// 	if(!fluide[i])
// 	//std::cout<<"nr triangle "<<j<<" contact fluide "<<fluide[j]<<"  points interface:= "<<interface[j][i]<<std::endl;
// 	std::cout<<"triangles size"<<Triangles_interface[i].size()<<std::endl;
// }

}

inline double signe(const double x)
{
  return (x < 0.) ? -1. : 1. ;
}

void Particule::solve_position(double dt){
  double rot[3][3];
  if(fixe){
    Dx = Vector_3(0.,0.,0.);
    Dxprev = Vector_3(0.,0.,0.);
    u = Vector_3(0.,0.,0.);
    u_half = u;
    e = Vector_3(0.,0.,0.);
    rot[0][0] = rot[1][1]= rot[2][2] =1.;
    rot[0][1]=rot[0][2]=rot[1][0]=rot[1][2]=rot[2][0]=rot[2][1]=0.;
    eprev = Vector_3(0.,0.,0.);
    omega = Vector_3(0.,0.,0.);
    omega_half = omega;
  } else {
    Dxprev = Dx;
    u = u+(Fi+Ff)/2.*(dt/m);
    u_half = u;
    Dx = Dx+u*dt;
		//Tests pour v�rifier qu'on a toujours une matrice de rotation
		for(int i=0;i<3;i++){
			double norm = rotref[i][0]*rotref[i][0]+rotref[i][1]*rotref[i][1]+rotref[i][2]*rotref[i][2];
			if(abs(norm-1.)>eps){
				cout << "Matrice de rotation rotref renomalisee" <<endl;
				cout << "Ligne " << i << " norme=" << norm << endl;
			}
		}
		double vectrot1 = rotref[0][2]-(rotref[1][0]*rotref[2][1]-rotref[2][0]*rotref[1][1]);
		double vectrot2 = rotref[1][2]-(rotref[2][0]*rotref[0][1]-rotref[0][0]*rotref[2][1]);
		double vectrot3 = rotref[2][2]-(rotref[0][0]*rotref[1][1]-rotref[1][0]*rotref[0][1]);
		if(vectrot1*vectrot1+vectrot2*vectrot2+vectrot3*vectrot3>eps){
			cout << "Erreur rotation rotref " << vectrot1 << " " << vectrot2 << " " << vectrot3 << endl;
			//getchar();
		}
    //Calcul de la matrice de rotation totale depuis le rep�re inertiel jusqu'au temps t et stockage de eprev
    double Q[3][3];
    double e0 = sqrt(1.-CGAL::to_double(e.squared_length()));
    //Recuperation de la matrice de rotation
    rot[0][0] = 1.-2.*CGAL::to_double(e.operator[](1)*e.operator[](1)+e.operator[](2)*e.operator[](2));
    rot[0][1] = 2.*CGAL::to_double(-e0*e.operator[](2)+e.operator[](0)*e.operator[](1));
    rot[0][2] = 2.*CGAL::to_double(e0*e.operator[](1)+e.operator[](0)*e.operator[](2));
    rot[1][0] = 2.*CGAL::to_double(e0*e.operator[](2)+e.operator[](1)*e.operator[](0));
    rot[1][1] = 1.-2.*CGAL::to_double(e.operator[](0)*e.operator[](0)+e.operator[](2)*e.operator[](2));
    rot[1][2] = 2.*CGAL::to_double(-e0*e.operator[](0)+e.operator[](1)*e.operator[](2));
    rot[2][0] = 2.*CGAL::to_double(-e0*e.operator[](1)+e.operator[](2)*e.operator[](0));
    rot[2][1] = 2.*CGAL::to_double(e0*e.operator[](0)+e.operator[](2)*e.operator[](1));
    rot[2][2] = 1.-2.*CGAL::to_double(e.operator[](0)*e.operator[](0)+e.operator[](1)*e.operator[](1));
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	Q[i][j] = rot[i][0]*rotref[0][j];
	Q[i][j] += rot[i][1]*rotref[1][j];
	Q[i][j] += rot[i][2]*rotref[2][j];
      }
    }    
    eprev = e;
    //Recuperation du e global a partir de omega
    double Omega[3];
    Omega[0] = Omega[1] = Omega[2] = 0.;
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	Omega[j] += CGAL::to_double(omega.operator[](k)*Q[k][j]);
      }
    }
    //cout << "debut " << Omega[0] << " " << Omega[1] << " " << Omega[2] << endl;
    //getchar();
    double norm2 = dt*dt*(Omega[0]*Omega[0]+Omega[1]*Omega[1]+Omega[2]*Omega[2]);
    if(norm2>1.){
      cout << "pas de temps trop grand : dt=" << dt << " Omega=" << Omega[0] << " " << Omega[1] << " " << Omega[2] << endl;
      getchar();
    }
    double eglob0 = 1.;//sqrt((1+sqrt(1-norm2))/2.);
    double eglob[3];
    for(int j=0;j<3;j++){
      eglob[j] = 0.;//dt*Omega[j]/2./eglob0;
    }
    //Recuperation de la matrice Zn
    double z[3][3];
    z[0][0] = 0.;
    z[0][1] = -Omega[2];
    z[0][2] = Omega[1];
    z[1][0] = Omega[2];
    z[1][1] = 0.;
    z[1][2] = -Omega[0];
    z[2][0] = -Omega[1];
    z[2][1] = Omega[0];
    z[2][2] = 0.;
    //Calcul de la matrice A
    double a[3];
    double d1 = (I[0]+I[1]+I[2])/2.-I[0];
    double d2 = (I[0]+I[1]+I[2])/2.-I[1];
    double d3 = (I[0]+I[1]+I[2])/2.-I[2];
    //Calcul du moment dans le rep�re inertiel
    double Mx = CGAL::to_double(Q[0][0]*((Mi+Mf).operator[](0))+Q[1][0]*((Mi+Mf).operator[](1))+Q[2][0]*((Mi+Mf).operator[](2)));
    double My = CGAL::to_double(Q[0][1]*((Mi+Mf).operator[](0))+Q[1][1]*((Mi+Mf).operator[](1))+Q[2][1]*((Mi+Mf).operator[](2)));
    double Mz = CGAL::to_double(Q[0][2]*((Mi+Mf).operator[](0))+Q[1][2]*((Mi+Mf).operator[](1))+Q[2][2]*((Mi+Mf).operator[](2)));
    a[0] = -(I[0]*z[1][2]-dt/2.*Mx);
    a[1] = (I[1]*z[0][2]+dt/2.*My);
    a[2] = -(I[2]*z[0][1]-dt/2.*Mz);
    //R�solution du probl�me non lin�aire
    double etemp0 = 1.;
    double etemp1 = 0.;
    double etemp2 = 0.;
    double etemp3 = 0.;
    double err1 = 1.;
    double err2 = 1.;
    double err3 = 1.;
    double epsilon = 1.e-15;
    int k=0;
    for(k=0; k<1000 && (err1>epsilon || err2>epsilon || err3>epsilon); k++){
      double x1 = (dt*a[0]-2.*(d2-d3)*etemp2*etemp3)/(2.*(d2+d3)*etemp0);
      double x2 = (dt*a[1]-2.*(d3-d1)*etemp1*etemp3)/(2.*(d1+d3)*etemp0);
      double x3 = (dt*a[2]-2.*(d1-d2)*etemp1*etemp2)/(2.*(d1+d2)*etemp0);
      etemp1 = x1;
      etemp2 = x2;
      etemp3 = x3;
// 	  //Test : on fixe la rotation en x et y
// 			etemp1 = 0.;
// 			etemp2 = 0.;
// 			//fin test 
      if(etemp1*etemp1+etemp2*etemp2+etemp3*etemp3>0.5){
	    etemp1 /=2.;
	    etemp2 /=2.;
	    etemp3 /=2.;
	    etemp0 = sqrt(1.-etemp1*etemp1-etemp2*etemp2-etemp3*etemp3);
      }
      else{etemp0 = sqrt(1.-etemp1*etemp1-etemp2*etemp2-etemp3*etemp3);}
      err1 = fabs((dt*a[0]-2.*(d2-d3)*etemp2*etemp3)/(2.*(d2+d3)*etemp0)-etemp1);
      err2 = fabs((dt*a[1]-2.*(d3-d1)*etemp1*etemp3)/(2.*(d1+d3)*etemp0)-etemp2);
      err3 = fabs((dt*a[2]-2.*(d1-d2)*etemp1*etemp2)/(2.*(d1+d2)*etemp0)-etemp3);
    }
    //cout << k << endl;
    eglob[0] = etemp1;
    eglob[1] = etemp2;
    eglob[2] = etemp3;
    eglob0 = etemp0;
    //Reconstruction de Z^n+1/2
    z[0][0] = (-2.*(eglob[1]*eglob[1]+eglob[2]*eglob[2]))/dt;
    z[0][1] = (-2.*eglob0*eglob[2]+2.*eglob[0]*eglob[1])/dt;
    z[0][2] = (2.*eglob0*eglob[1]+2.*eglob[0]*eglob[2])/dt;
    z[1][0] = (2.*eglob0*eglob[2]+2.*eglob[0]*eglob[1])/dt;
    z[1][1] = (-2.*(eglob[0]*eglob[0]+eglob[2]*eglob[2]))/dt;
    z[1][2] = (-2.*eglob0*eglob[0]+2.*eglob[1]*eglob[2])/dt;
    z[2][0] = (-2.*eglob0*eglob[1]+2.*eglob[0]*eglob[2])/dt;
    z[2][1] = (2.*eglob0*eglob[0]+2.*eglob[1]*eglob[2])/dt;
    z[2][2] = (-2.*(eglob[0]*eglob[0]+eglob[1]*eglob[1]))/dt;
    //Update de la matrice Q
    double Qprev[3][3];
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	       Qprev[i][j] = Q[i][j];
      }
    }
    Q[0][0] = Qprev[0][0]*(1.+dt*z[0][0])+Qprev[0][1]*dt*z[1][0]+Qprev[0][2]*dt*z[2][0];
    Q[0][1] = Qprev[0][0]*dt*z[0][1]+Qprev[0][1]*(1.+dt*z[1][1])+Qprev[0][2]*dt*z[2][1];
    Q[0][2] = Qprev[0][0]*dt*z[0][2]+Qprev[0][1]*dt*z[1][2]+Qprev[0][2]*(1.+dt*z[2][2]);
    Q[1][0] = Qprev[1][0]*(1.+dt*z[0][0])+Qprev[1][1]*dt*z[1][0]+Qprev[1][2]*dt*z[2][0];
    Q[1][1] = Qprev[1][0]*dt*z[0][1]+Qprev[1][1]*(1.+dt*z[1][1])+Qprev[1][2]*dt*z[2][1];
    Q[1][2] = Qprev[1][0]*dt*z[0][2]+Qprev[1][1]*dt*z[1][2]+Qprev[1][2]*(1.+dt*z[2][2]);
    Q[2][0] = Qprev[2][0]*(1.+dt*z[0][0])+Qprev[2][1]*dt*z[1][0]+Qprev[2][2]*dt*z[2][0];
    Q[2][1] = Qprev[2][0]*dt*z[0][1]+Qprev[2][1]*(1.+dt*z[1][1])+Qprev[2][2]*dt*z[2][1];
    Q[2][2] = Qprev[2][0]*dt*z[0][2]+Qprev[2][1]*dt*z[1][2]+Qprev[2][2]*(1.+dt*z[2][2]);
    //Tests pour v�rifier qu'on a toujours une matrice de rotation
    for(int i=0;i<3;i++){
      double norm = Q[i][0]*Q[i][0]+Q[i][1]*Q[i][1]+Q[i][2]*Q[i][2];
      Q[i][0] /= norm;
      Q[i][1] /= norm;
      Q[i][2] /= norm;
      if(abs(norm-1.)>eps){
	cout << "Matrice de rotation renomalisee" <<endl;
	cout << "Ligne " << i << " norme=" << norm << endl;
      }
    }
    double vect1 = Q[0][2]-(Q[1][0]*Q[2][1]-Q[2][0]*Q[1][1]);
    double vect2 = Q[1][2]-(Q[2][0]*Q[0][1]-Q[0][0]*Q[2][1]);
    double vect3 = Q[2][2]-(Q[0][0]*Q[1][1]-Q[1][0]*Q[0][1]);
    if(vect1*vect1+vect2*vect2+vect3*vect3>eps){
      cout << "Erreur rotation " << vect1 << " " << vect2 << " " << vect3 << endl;
      //getchar();
    }
    //R�cup�ration de la matrice de rotation de la particule
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	rot[i][j] = Q[i][0]*rotref[j][0];
	rot[i][j] += Q[i][1]*rotref[j][1];
	rot[i][j] += Q[i][2]*rotref[j][2];
      }
    }
    //Calcul de e a partir de rot
    double q1 = rot[2][1]-rot[1][2];
    double q2 = rot[0][2]-rot[2][0];
    double q3 = rot[1][0]-rot[0][1];
    double e1,e2,e3;
    e1 = signe(q1)*sqrt(max((1.+rot[0][0]-rot[1][1]-rot[2][2])/4.,0.));
    e2 = signe(q2)*sqrt(max((1.+rot[1][1]-rot[0][0]-rot[2][2])/4.,0.));
    e3 = signe(q3)*sqrt(max((1.+rot[2][2]-rot[0][0]-rot[1][1])/4.,0.));
    e = Vector_3(e1,e2,e3);
    
    //Calcul de Omega^n+1/2
    double omega1 = 0.;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	omega1 -= 1./2.*Qprev[1][i]*z[i][j]*(Qprev[2][j]+Q[2][j]);
      }
    }
    double omega2 = 0.;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	omega2 += 1./2.*Qprev[0][i]*z[i][j]*(Qprev[2][j]+Q[2][j]);
      }
    }
    double omega3 = 0.;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	omega3 -= 1./2.*Qprev[0][i]*z[i][j]*(Qprev[1][j]+Q[1][j]);
      }
    }
    omega = Vector_3(omega1,omega2,omega3);
    omega_half = omega;
  }//Fin du calcul dans le cas d'une particule libre
  /*//Test de fixer la rotation
  rot[0][0]= rot[1][1] = rot[2][2] =1.;
	rot[0][1] = rot[0][2] =rot[1][0] = rot[1][2] = rot[2][0] = rot[2][1] = 0.;
	rotprev[0][0]= rotprev[1][1] = rotprev[2][2] =1.;
	rotprev[0][1] = rotprev[0][2] =rotprev[1][0] = rotprev[1][2] = rotprev[2][0] = rotprev[2][1] = 0.;
	omega = Vector_3(0.,0.,0.);
	omega_half = omega;
	//fin test */
  //Mise � jour de la transformation donnant le mouvement de la particule
  mvt_tprev = mvt_t;
  Aff_transformation_3 rotation(rot[0][0],rot[0][1],rot[0][2],rot[1][0],rot[1][1],rot[1][2],rot[2][0],rot[2][1],rot[2][2]);
  Aff_transformation_3 translation(CGAL::TRANSLATION,Vector_3(Point_3(0.,0.,0.),x0)+Dx);
  Aff_transformation_3 translation_inv(CGAL::TRANSLATION,Vector_3(x0,Point_3(0.,0.,0.))-Dx);
  mvt_t = translation*(rotation*translation_inv);
}

void Particule::solve_vitesse(double dt){
  if(fixe){
    u = Vector_3(0.,0.,0.);
    omega = Vector_3(0.,0.,0.);
  } else {
    u = u+(Fi+Ff)/2.*(dt/m);
    //Calcul de la matrice de rotation totale depuis le rep�re inertiel jusqu'au temps t
    double Q[3][3];
    double e0 = sqrt(1.-CGAL::to_double(e.squared_length()));
    //Recuperation de la matrice de rotation
    double rot[3][3];
    rot[0][0] = 1.-2.*CGAL::to_double(e.operator[](1)*e.operator[](1)+e.operator[](2)*e.operator[](2));
    rot[0][1] = 2.*CGAL::to_double(-e0*e.operator[](2)+e.operator[](0)*e.operator[](1));
    rot[0][2] = 2.*CGAL::to_double(e0*e.operator[](1)+e.operator[](0)*e.operator[](2));
    rot[1][0] = 2.*CGAL::to_double(e0*e.operator[](2)+e.operator[](1)*e.operator[](0));
    rot[1][1] = 1.-2.*CGAL::to_double(e.operator[](0)*e.operator[](0)+e.operator[](2)*e.operator[](2));
    rot[1][2] = 2.*CGAL::to_double(-e0*e.operator[](0)+e.operator[](1)*e.operator[](2));
    rot[2][0] = 2.*CGAL::to_double(-e0*e.operator[](1)+e.operator[](2)*e.operator[](0));
    rot[2][1] = 2.*CGAL::to_double(e0*e.operator[](0)+e.operator[](2)*e.operator[](1));
    rot[2][2] = 1.-2.*CGAL::to_double(e.operator[](0)*e.operator[](0)+e.operator[](1)*e.operator[](1));
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	Q[i][j] = rot[i][0]*rotref[0][j];
	Q[i][j] += rot[i][1]*rotref[1][j];
	Q[i][j] += rot[i][2]*rotref[2][j];
      }
    }    
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	Q[i][j] = rot[i][0]*rotref[0][j];
	Q[i][j] += rot[i][1]*rotref[1][j];
	Q[i][j] += rot[i][2]*rotref[2][j];
      }
    }    
    //Recuperation de Zn+1/2 � partir de omega
    double Omega[3];
    Omega[0] = Omega[1] = Omega[2] = 0.;
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	Omega[j] += CGAL::to_double(omega.operator[](k)*Q[k][j]);
      }
    }
    //cout << "debut " << Omega[0] << " " << Omega[1] << " " << Omega[2] << endl;
    //getchar();
    double norm2 = dt*dt*(Omega[0]*Omega[0]+Omega[1]*Omega[1]+Omega[2]*Omega[2]);
    if(norm2>1.){
      cout << "pas de temps trop grand : dt=" << dt << " Omega=" << Omega[0] << " " << Omega[1] << " " << Omega[2] << endl;
      getchar();
    }
    double eglob0 = sqrt((1.+sqrt(1.-norm2))/2.);
    double eglob[3];
    for(int j=0;j<3;j++){
      eglob[j] = dt*Omega[j]/2./eglob0;
    }
    //Recuperation de la matrice Zn+1/2
    //double eglob0 = sqrt(1.-eglob[0]*eglob[0]-eglob[1]*eglob[1]-eglob[2]*eglob[2]);
    double z[3][3];
    z[0][0] = (-2.*(eglob[1]*eglob[1]+eglob[2]*eglob[2]))/dt;
    z[0][1] = (-2.*eglob0*eglob[2]+2.*eglob[0]*eglob[1])/dt;
    z[0][2] = (2.*eglob0*eglob[1]+2.*eglob[0]*eglob[2])/dt;
    z[1][0] = (2.*eglob0*eglob[2]+2.*eglob[0]*eglob[1])/dt;
    z[1][1] = (-2.*(eglob[0]*eglob[0]+eglob[2]*eglob[2]))/dt;
    z[1][2] = (-2.*eglob0*eglob[0]+2.*eglob[1]*eglob[2])/dt;
    z[2][0] = (-2.*eglob0*eglob[1]+2.*eglob[0]*eglob[2])/dt;
    z[2][1] = (2.*eglob0*eglob[0]+2.*eglob[1]*eglob[2])/dt;
    z[2][2] = (-2.*(eglob[0]*eglob[0]+eglob[1]*eglob[1]))/dt;
    //Calcul de la matrice A
    double a[3];
    double d1 = (I[0]+I[1]+I[2])/2.-I[0];
    double d2 = (I[0]+I[1]+I[2])/2.-I[1];
    double d3 = (I[0]+I[1]+I[2])/2.-I[2];
    //Calcul du moment dans le rep�re inertiel
    double Mx = CGAL::to_double(Q[0][0]*((Mi+Mf).operator[](0))+Q[1][0]*((Mi+Mf).operator[](1))+Q[2][0]*((Mi+Mf).operator[](2)));
    double My = CGAL::to_double(Q[0][1]*((Mi+Mf).operator[](0))+Q[1][1]*((Mi+Mf).operator[](1))+Q[2][1]*((Mi+Mf).operator[](2)));
    double Mz = CGAL::to_double(Q[0][2]*((Mi+Mf).operator[](0))+Q[1][2]*((Mi+Mf).operator[](1))+Q[2][2]*((Mi+Mf).operator[](2)));
    a[0] = -(d2*z[1][2]-d3*z[2][1]-dt/2.*Mx);
    a[1] = (d1*z[0][2]-d3*z[2][0]+dt/2.*My);
    a[2] = -(d1*z[0][1]-d2*z[1][0]-dt/2.*Mz);
    //R�solution du probl�me lin�aire sur Zn+1
    z[0][0] = 0.;
    z[0][1] = -a[2]/I[2];
    z[0][2] = a[1]/I[1];
    z[1][0] = -z[0][1];
    z[1][1] = 0.;
    z[1][2] = -a[0]/I[0];
    z[2][0] = -z[0][2];
    z[2][1] = -z[1][2];
    z[2][2]= 0.;
    //Calcul de Omega^n+1
    double omega1 = 0.;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	omega1 -= Q[1][i]*z[i][j]*Q[2][j];
      }
    }
    double omega2 = 0.;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	omega2 += Q[0][i]*z[i][j]*Q[2][j];
      }
    }
    double omega3 = 0.;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	omega3 -= Q[0][i]*z[i][j]*Q[1][j];
      }
    }
//     //Test pour fixer les composantes x et y de la rotation
// 		omega1 = 0.;
// 		omega2 = 0.;
// 		//fin test 
    omega = Vector_3(omega1,omega2,omega3);
		/*//Test de fixer la rotation
		rot[0][0]= rot[1][1] = rot[2][2] =1.;
		rot[0][1] = rot[0][2] =rot[1][0] = rot[1][2] = rot[2][0] = rot[2][1] = 0.;
		rotprev[0][0]= rotprev[1][1] = rotprev[2][2] =1.;
		rotprev[0][1] = rotprev[0][2] =rotprev[1][0] = rotprev[1][2] = rotprev[2][0] = rotprev[2][1] = 0.;
		omega = Vector_3(0.,0.,0.);
		omega_half = omega;
		//fin test */
  }//Fin du calcul dans le cas d'une particule libre
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


Vector_3 Particule::vitesse_parois(Point_3& X_f){
		
	Vector_3 V_f = u_half + cross_product(omega_half, Vector_3(Point_3(x0.operator[](0) + Dx.operator[](0), x0.operator[](1) + Dx.operator[](1),x0.operator[](2) + Dx.operator[](2)),X_f));

	return V_f;
}	

Vector_3 Particule::vitesse_parois_prev(Point_3& X_f){
	
	Vector_3 V_f = u_half + cross_product(omega_half, Vector_3(Point_3(x0.operator[](0) + Dxprev.operator[](0), x0.operator[](1) + Dxprev.operator[](1),x0.operator[](2) + Dxprev.operator[](2)),X_f));
	
	return V_f;
}	
void Face::compProjectionIntegrals(double &P1, double &Pa, double &Pb, double &Paa, double &Pab, double &Pbb, double &Paaa, double &Paab, double &Pabb, double &Pbbb, int a, int b, int c){
  //Utilisation de la fonction decrite par Brian Mirtich 1996 (cf www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar)
  P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.;
  for(int i=0;i<size();i++){
    int j= (i+1)%(size());
    double a0 = CGAL::to_double(vertex[i].pos.operator[](a));
    double b0 = CGAL::to_double(vertex[i].pos.operator[](b));
    double a1 = CGAL::to_double(vertex[j].pos.operator[](a));
    double b1 = CGAL::to_double(vertex[j].pos.operator[](b));
    double Da = a1-a0;
    double Db = b1-b0;
    double a02 = a0*a0;
    double a03 = a0*a02;
    double a04 = a0*a03;
    double b02 = b0*b0;
    double b03 = b0*b02;
    double b04 = b0*b03;
    double a12 = a1*a1;
    double a13 = a1*a12;
    double b12 = b1*b1;
    double b13 = b1*b12;
    double C1 = a1+a0;
    double Ca = a1*C1+a02;
    double Caa = a1*Ca+a03;
    double Caaa = a1*Caa+a04;
    double Cb = b12+b1*b0+b02;
    double Cbb = b1*Cb+b03;
    double Cbbb = b1*Cbb+b04;
    double Cab = 3*a12+2*a1*a0+a02;
    double Kab = a12+2*a1*a0+3*a02;
    double Caab = a0*Cab+4*a13;
    double Kaab = a1*Kab+4*a03;
    double Cabb = 4*b13+3*b12*b0+2*b1*b02+b03;
    double Kabb = b13+2*b12*b0+3*b1*b02+4*b03;
    P1 += Db*C1;
    Pa += Db*Ca; Paa += Db*Caa; Paaa += Db*Caaa;
    Pb += Da*Cb; Pbb += Da*Cbb; Pbbb += Da*Cbbb;
    Pab += Db*(b1*Cab+b0*Kab);
    Paab += Db*(b1*Caab+b0*Kaab);
    Pabb += Da*(a1*Cabb+a0*Kabb);
  }
  P1 /= 2.;
  Pa /= 6.; Paa /= 12.; Paaa /= 20.;
  Pb /=-6.; Pbb /=-12.; Pbbb /=-20.;
  Pab /= 24.;
  Paab /= 60.;
  Pabb /= 60.;
}

void Face::compFaceIntegrals(double &Fa, double &Fb, double &Fc, double &Faa, double &Fbb, double &Fcc, double &Faaa, double &Fbbb, double &Fccc, double &Faab, double &Fbbc, double &Fcca, double na, double nb, double nc, int a, int b, int c){
  //Utilisation de la fonction decrite par Brian Mirtich 1996 (cf www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar)
  double P1,Pa,Pb,Paa,Pab,Pbb,Paaa,Paab,Pabb,Pbbb;
  Vector_3 p(Point_3(0,0,0),vertex[0].pos);
  double w = -CGAL::to_double(normale*p);
  double k1 = 1./nc;
  double k2 = k1*k1;
  double k3 = k1*k2;
  double k4 = k1*k3;
  compProjectionIntegrals(P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb, a, b, c);
  Fa = k1*Pa;
  Fb = k1*Pb;
  Fc = -k2*(na*Pa+nb*Pb+w*P1);
  Faa = k1*Paa;
  Fbb = k1*Pbb;
  Fcc = k3*(na*na*Paa+2*na*nb*Pab+nb*nb*Pbb+2*na*w*Pa+2*nb*w*Pb+w*w*P1);
  Faaa = k1*Paaa;
  Fbbb = k1*Pbbb;
  Fccc = -k4*(na*na*na*Paaa+3*na*na*nb*Paab+3*na*nb*nb*Pabb+nb*nb*nb*Pbbb+3*na*na*w*Paa+6*na*nb*w*Pab+3*nb*nb*w*Pbb+3*na*w*w*Pa+3*nb*w*w*Pb+w*w*w*P1);
  Faab = k1*Paab;
  Fbbc = -k2*(na*Pabb+nb*Pbbb+w*Pbb);
  Fcca = k3*(na*na*Paaa+2*na*nb*Paab+nb*nb*Pabb+2*na*w*Paa+2*nb*w*Pab+w*w*Pa);
}

void Particule::CompVolumeIntegrals(double &T1, double &Tx, double &Ty, double &Tz, double &Txx, double &Tyy, double &Tzz, double &Txy, double &Tyz, double &Tzx){
  //Utilisation de la fonction decrite par Brian Mirtich 1996 (cf www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar)
  T1 = Tx=Ty=Tz=Txx=Tyy=Tzz=Txy=Tyz=Tzx=0.;
  for(int i=0;i<size();i++){
    double Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxxx,Fyyy,Fzzz,Fxxy,Fyyz,Fzzx;
    double nx,ny,nz,na,nb,nc;
    nx=CGAL::to_double(faces[i].normale.operator[](0));
    ny=CGAL::to_double(faces[i].normale.operator[](1));
    nz=CGAL::to_double(faces[i].normale.operator[](2));
    //Choix d'une permutation orientee abc telle que nc soit maximale
    if(abs(nx)>abs(ny)){
      if(abs(nx)>abs(nz)){
	//Cas a=y,b=z,c=x
	na = ny; nb = nz; nc = nx;
	faces[i].compFaceIntegrals(Fy, Fz, Fx, Fyy, Fzz, Fxx, Fyyy, Fzzz, Fxxx, Fyyz, Fzzx, Fxxy,na,nb,nc,1,2,0);
      }
      else {
	//Cas a=x,b=y,c=z
	na = nx; nb = ny; nc = nz;
	faces[i].compFaceIntegrals(Fx, Fy, Fz, Fxx, Fyy, Fzz, Fxxx, Fyyy, Fzzz, Fxxy, Fyyz, Fzzx,na,nb,nc,0,1,2);
      }
    }
    else {
      if(abs(ny)>abs(nz)){
	//Cas a=z,b=x,c=y
	na = nz; nb = nx; nc = ny;
	faces[i].compFaceIntegrals(Fz, Fx, Fy, Fzz, Fxx, Fyy, Fzzz, Fxxx, Fyyy, Fzzx, Fxxy, Fyyz,na,nb,nc,2,0,1);
      }
      else{
	//Cas a=x,b=y,c=z
	na = nx; nb = ny; nc = nz;
	faces[i].compFaceIntegrals(Fx, Fy, Fz, Fxx, Fyy, Fzz, Fxxx, Fyyy, Fzzz, Fxxy, Fyyz, Fzzx,na,nb,nc,0,1,2);
      }
    }
    //Calcul des integrales
    T1 += nx*Fx;
    Tx += nx*Fxx;
    Ty += ny*Fyy;
    Tz += nz*Fzz;
    Txx += nx*Fxxx;
    Tyy += ny*Fyyy;
    Tzz += nz*Fzzz;
    Txy += nx*Fxxy;
    Tyz += ny*Fyyz;
    Tzx += nz*Fzzx;
  }
  Tx /= 2.;
  Ty /= 2.;
  Tz /= 2.;
  Txx /=3.;
  Tyy /=3.;
  Tzz /=3.;
  Txy /=2.;
  Tyz /=2.;
  Tzx /=2.;
}

//Resolution de l'equation de degre 3 ax3+bx2+cx+d=0, de solutions x1, x2 et x3
void solve_eq3(double a, double b, double c, double d, double &x1, double &x2, double &x3){
  //On trouve x1 par une methode de Newton
  x1 = 0.;
  double res = a*x1*x1*x1+b*x1*x1+c*x1+d;
  int i;
  for(i=0;i<10000 && res!=0.;i++){
    double p = 3.*a*x1*x1+2.*b*x1+c;
    if(abs(p)>1.e-5*abs(a)){
      x1 -= res/p;
    }
    else {
      if(abs(res)>1.e-4*abs(a)){
	x1 += pow(2.,-i);
      }
      else {
	double Delta = max(4.*b*b-12.*a*c,0.);
	double xm = -(2.*b+sqrt(Delta))/(6.*a);
	double xM = -(2.*b-sqrt(Delta))/(6.*a);
	if(abs(x1-xm)<abs(x1-xM)){
	  x1 = xm;
	}
	else {
	  x1 = xM;
	}
      }
    }
    res = a*x1*x1*x1+b*x1*x1+c*x1+d;
  }
  double a1 = a;
  double b1 = b+a*x1;
  double c1 = c+b*x1+a*x1*x1;
  //On resout a1*x2+b1*x+c1=0
  double Delta = max(b1*b1-4.*a1*c1,0.);
  if(Delta<eps*a*a){
    Delta = 0.;
  }
  x2 = (-b1+sqrt(Delta))/(2.*a1);
  x3 = (-b1-sqrt(Delta))/(2.*a1);
  if(abs(x1-x2)<2.e-3*x1){
    if(Delta < eps*a*a){
      //On vérifie si les trois solutions ne seraient pas identiques
      double b2 = a*3.*x2;
      double c2 = a*3.*x2*x2;
      double d2 = a*pow(x2,3);
      if(abs(b2-b)<eps*a && abs(c2-c)<eps*a && abs(d2-d)<eps*a){
	x1 = x3 = x2;
      }
    }
    else{
      if(abs(a*x1*x1*x1+b*x1*x1+c*x1+d)>abs(a*x2*x2*x2+b*x2*x2+c*x2+d)){
	x1 = x3;
	x3 = x2;
      }
    }
  }
  else if(abs(x1-x3)<1.e-3*x1){
    if(abs(a*x1*x1*x1+b*x1*x1+c*x1+d)>abs(a*x3*x3*x3+b*x3*x3+c*x3+d)){
      x1 = x2;
      x2 = x3;
    }
  }
}

void Particule::Inertie(){
  double T1,Tx,Ty,Tz,Txx,Tyy,Tzz,Txy,Tyz,Tzx;
  CompVolumeIntegrals(T1,Tx,Ty,Tz,Txx,Tyy,Tzz,Txy,Tyz,Tzx);
  double R[3][3];
  double xG = CGAL::to_double(x0.operator[](0));
  double yG = CGAL::to_double(x0.operator[](1));
  double zG = CGAL::to_double(x0.operator[](2));
  R[0][0] = rhos*(Tyy-2.*yG*Ty+yG*yG*T1+Tzz-2.*zG*Tz+zG*zG*T1);
  R[1][0] = R[0][1] = rhos*(Txy-yG*Tx-xG*Ty+xG*yG*T1);
  R[2][0] = R[0][2] = rhos*(Tzx-zG*Tx-xG*Tz+xG*zG*T1);
  R[1][1] = rhos*(Txx-2.*xG*Tx+xG*xG*T1+Tzz-2.*zG*Tz+zG*zG*T1);
  R[1][2] = R[2][1] = rhos*(Tyz-zG*Ty-yG*Tz+yG*zG*T1);
  R[2][2] = rhos*(Tyy-2.*yG*Ty+yG*yG*T1+Txx-2.*xG*Tx+xG*xG*T1);
  double A = R[0][0];
  double B = R[1][1];
  double C = R[2][2];
  double D = -R[1][2];
  double E = -R[0][2];
  double F = -R[0][1];
// 	cout << A << " " << -F << " " << -E << endl;
// 	cout << -F << " " << B << " " << -D << endl;
// 	cout << -E << " " << -D << " " << C << endl;
  //Masse et volume
  V = T1;
  m = rhos*T1;
  if(m<eps){
    cout<< "masse nulle " << m << endl;
    getchar();
  }
  //Calcul des moments d'inertie
  double a = -1.;
  double b = R[0][0] + R[1][1] + R[2][2];
  double c = R[0][1]*R[0][1] + R[0][2]*R[0][2] + R[1][2]*R[1][2] - R[0][0]*R[1][1] - R[0][0]*R[2][2] - R[1][1]*R[2][2];
  double d = R[0][0]*R[1][1]*R[2][2]-R[0][0]*R[1][2]*R[1][2]-R[1][1]*R[0][2]*R[0][2]-R[2][2]*R[0][1]*R[0][1]-2.*R[0][1]*R[0][2]*R[1][2];
	//cout << a << " " << b << " " << c << " " << d << endl;
  solve_eq3(a,b,c,d,I[0],I[1],I[2]);
  //Calcul des vecteurs propres associes
  if(abs(I[1]-I[2])>1.e-5*I[1]){
    for(int i=0;i<3;i++){
      double ux,uy,uz;
      if(abs(A-I[i])>eps){
	if(abs(B-I[i])>eps){
	  uz = (A-I[i])*(B-I[i])-F*F;
	  uy = D*(A-I[i])+E*F;
	  ux = (F*uy+E*uz)/(A-I[i]);
	} else {
	  uy = (A-I[i])*(C-I[i])-E*E;
	  uz = D*(A-I[i])+E*F;
	  ux = (F*uy+E*uz)/(A-I[i]);
	}
      }
      else if(abs(B-I[i])>eps){
	ux = (B-I[i])*(C-I[i])-D*D;
	uz = E*(B-I[i])+D*F;
	uy = (D*uz+F*ux)/(B-I[i]);
      }
      else if(abs(C-I[i])>eps){
	uy = (C-I[i])*(A-I[i])-E*E;
	ux = F*(C-I[i])+D*E;
	uz = (E*ux+D*uy)/(C-I[i]);
      }
      else{
	if(abs(D)<eps){
	  ux = 0.;
	  uy = -E;
	  uz = F;
	}
	else if(abs(E)<eps){
	  uy = 0.;
	  uz = -F;
	  ux = D;
	}
	else {
	  uz = 0.;
	  ux = -D;
	  uy = E;
	}
			}
      double norm = std::sqrt(ux*ux+uy*uy+uz*uz);
      ux /= norm;
      uy /= norm;
      uz /= norm;
      rotref[0][i] = ux;
      rotref[1][i] = uy;
      rotref[2][i] = uz;
      for(int j=0;j<3;j++){
	if(rotref[j][i]!=rotref[j][i]){
	  cout << "rotref "<< rotref[j][i] << " " << j << " " << i << " " << norm << endl;
	  cout << A-I[i] << " " << -F << " " << -E << endl;
	  cout << -F << " " << B-I[i] << " " << -D << endl;
	  cout << -E << " " << -D << " " << C-I[i] << endl;
	  getchar();
	}
      }
    }
  }
  else{
		//cout << "I " << I[0] << " " << I[1] << " " << I[2] << endl;
    if(abs(I[0]-I[1])>1.e-5*I[1]){
      for(int i=0;i<2;i++){
	double ux,uy,uz;
	if(abs(A-I[i])>eps){
	  if(abs(B-I[i])>eps){
	    uz = (A-I[i])*(B-I[i])-F*F;
	    uy = D*(A-I[i])+E*F;
	    ux = (F*uy+E*uz)/(A-I[i]);
	  } else {
	    uz = 1.;
	    uy = 0.;
	    ux = (F*uy+E*uz)/(A-I[i]);
	  }
	}
	else if(abs(B-I[i])>eps){
	  ux = (B-I[i])*(C-I[i])-D*D;
	  uz = E*(B-I[i])+D*F;
	  uy = (D*uz+F*ux)/(B-I[i]);
	}
	else if(abs(C-I[i])>eps){
	  uy = (C-I[i])*(A-I[i])-E*E;
	  ux = F*(C-I[i])+D*E;
	  uz = (E*ux+D*uy)/(C-I[i]);
	}
	else{
	  if(abs(D)<eps){
	    ux = 0.;
	    uy = -E;
	    uz = F;
	  }
	  else if(abs(E)<eps){
	    uy = 0.;
	    uz = -F;
	    ux = D;
	  }
	  else {
	    uz = 0.;
	    ux = -D;
	    uy = E;
	  }
	}
	double norm = sqrt(ux*ux+uy*uy+uz*uz);
	ux /= norm;
	uy /= norm;
	uz /= norm;
	rotref[0][i] = ux;
	rotref[1][i] = uy;
	rotref[2][i] = uz;
	for(int j=0;j<3;j++){
	  if(rotref[j][i]!=rotref[j][i]){
	    cout << "rotref "<< rotref[j][i] << " " << j << " " << i << " " << norm << endl;
	    getchar();
	  }
	}
      }
      //Test : produit scalaire des deux premieres colonnes
      double scal = rotref[0][0]*rotref[0][1]+rotref[1][0]*rotref[1][1]+rotref[2][0]*rotref[2][1];
// 			cout << "produit scalaire " << scal << endl;
			if(abs(scal)>eps){
					//on prend un vecteur normal a rotref[][0]
					rotref[0][1] = rotref[1][0];
					rotref[1][1] = -rotref[0][0];
					rotref[2][1] = 0.;
					double norm = sqrt(rotref[0][1]*rotref[0][1]+rotref[1][1]*rotref[1][1]+rotref[2][1]*rotref[2][1]);
					if(norm>eps){
						rotref[0][1] /= norm;
						rotref[1][1] /= norm;
						rotref[2][1] /= norm;
					}
					else {
						rotref[0][1] = -rotref[2][0];
						rotref[1][1] = 0.;
						rotref[2][1] = rotref[0][0];
						norm = sqrt(rotref[0][1]*rotref[0][1]+rotref[1][1]*rotref[1][1]+rotref[2][1]*rotref[2][1]);
						rotref[0][1] /= norm;
						rotref[1][1] /= norm;
						rotref[2][1] /= norm;
					}
			}
      rotref[0][2] = rotref[1][0]*rotref[2][1]-rotref[2][0]*rotref[1][1];
      rotref[1][2] = rotref[2][0]*rotref[0][1]-rotref[0][0]*rotref[2][1];
      rotref[2][2] = rotref[0][0]*rotref[1][1]-rotref[1][0]*rotref[0][1];
    }
    else{
      rotref[0][0] = rotref[1][1] = rotref[2][2] = 1.;
      rotref[0][1] = rotref[1][0] = rotref[0][2] = rotref[2][0] = rotref[1][2] = rotref[2][1] = 0.;
    }
  }
  //Test 26/11/12
  rotref[0][0] = rotref[1][1] = rotref[2][2] = 1.;
	rotref[0][1] = rotref[1][0] = rotref[0][2] = rotref[2][0] = rotref[1][2] = rotref[2][1] = 0.;
	// Fin Test 26/11/12
	
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      if(rotref[i][j]!=rotref[i][j]){
				cout << "rotref "<< rotref[i][j] << " " << i << " " << j << endl;
				getchar();
      }
    }
  }
  
  //Test sur le determinant de la matrice de rotation (1 ou -1)
  double det = rotref[0][2]*(rotref[1][0]*rotref[2][1]-rotref[2][0]*rotref[1][1]);
  det += rotref[1][2]*(rotref[2][0]*rotref[0][1]-rotref[0][0]*rotref[2][1]);
  det += rotref[2][2]*(rotref[0][0]*rotref[1][1]-rotref[1][0]*rotref[0][1]);
  if(det<0.){
    for(int i=0;i<3;i++){
      rotref[i][2] *= -1.;
    }
  }
  for(int i=0;i<3;i++){
    int j = (i+1)%3;
    if(abs(rotref[0][i]*rotref[0][j]+rotref[1][i]*rotref[1][j]+rotref[2][i]*rotref[2][j])>eps){
      cout << "erreur dans le calcul des moments d'inertie" << endl;
    }
  }

  //Calcul des moments d'inertie des faces (pour le calcul des torsions)
  for(int i=0;i<size();i++){
    faces[i].Inertie();
  }
  /*Test
  double R0[3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      R0[i][j] = 0.;
      for(int k=0;k<3;k++){
	for(int l=0;l<3;l++){
	  R0[i][j] += rotref[k][i]*R[k][l]*rotref[l][j];
	}
      }
      cout << R0[i][j] << " " ;
    }
    cout << endl;
  }
  cout << I[0] << " " << I[1] << " " << I[2] << endl;
  getchar();
  //Fin du test */
}

void Particule::Volume_libre(){
  Vl = 0.;
  for(int i=0;i<size();i++){
    if(faces[i].voisin == -1){
      Vector_3 v1(faces[i].vertex[0].pos,faces[i].vertex[1].pos);
      Vector_3 v2(faces[i].vertex[0].pos,faces[i].vertex[2].pos);
      Vector_3 v3(x0,faces[i].vertex[0].pos);
      Vl += 1./6.*CGAL::to_double(cross_product(v1,v2)*v3);
    }
  }
}

void Solide::impression(int n){ //Sortie au format vtk
  int nb_part = solide.size();

  int nb_triangles = 0.;
  for(int it=0; it<nb_part; it++){
    nb_triangles += solide[it].triangles.size();
  }

  //const char* solidevtk;
  //{
  std::ostringstream oss;
  oss << "resultats/solide" << n << ".vtk";
  string s = oss.str();
  //cout << s << endl;
  const char* const solidevtk = s.c_str();
  //}

  //Ouverture des flux en donne en ecriture
  std::ofstream vtk;
  vtk.open(solidevtk,ios::out);
  if(vtk.is_open())
  {
    // cout <<"ouverture de xt.vtk reussie" << endl;
  } else {
    cout <<"ouverture de solide" << n << ".vtk rate" << endl;
  }
  vtk << setprecision(15);
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
  vtk << "VECTORS deplacement double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(int it=0; it<nb_part; it++){
    for(int l= 0; l<solide[it].triangles.size(); l++)
    {
      vtk << solide[it].Dx.operator[](0) << " " << solide[it].Dx.operator[](1) << " " << solide[it].Dx.operator[](2) << endl;
    }
  }
  vtk << "\n";
  //Vitesse
  vtk << "VECTORS vitesse double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(int it=0; it<nb_part; it++){
    for(int l= 0; l<solide[it].triangles.size(); l++)
    {
      vtk << solide[it].u.operator[](0) << " " << solide[it].u.operator[](1) << " " << solide[it].u.operator[](2) << endl;
    }
  }
  vtk << "\n";
  //Rotation en x
  vtk << "VECTORS e double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(int it=0; it<nb_part; it++){
    for(int l= 0; l<solide[it].triangles.size(); l++)
    {
      vtk << solide[it].e.operator[](0) << " " << solide[it].e.operator[](1) << " " << solide[it].e.operator[](2) << endl;
    }
  }
  vtk << "\n";
  //Vitesse de rotation
  vtk << "VECTORS omega double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(int it=0; it<nb_part; it++){
    for(int l= 0; l<solide[it].triangles.size(); l++)
    {
      vtk << solide[it].omega.operator[](0) << " " << solide[it].omega.operator[](1) << " " << solide[it].omega.operator[](2) << endl;
    }
  }
  vtk << "\n";
  vtk.close();
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
      Faces[j] = Face(Vertex, voisin);
    }
    P[i] = Particule(centre, xmin, ymin, zmin, xmax, ymax, zmax, Faces);
    P[i].fixe = fixe;
    P[i].u = Vector_3(u,v,w);
    P[i].omega = Vector_3(theta,phi,psi);
		P[i].u_half = Vector_3(u,v,w);
		P[i].omega_half = Vector_3(theta,phi,psi);
  }
  //Boucle de mise a jour des particules sur les sommets du maillage
  //Mise a jour des distances a l'equilibre entre particules en meme temps
  for(int i=0;i<P.size();i++){
    for(int j=0;j<P[i].size();j++){
      for(int k=0;k<P[i].faces[j].size();k++){
	for(int l=0;l<P.size();l++){
	  if(points_particules[P[i].faces[j].vertex[k].num][l]){
	    P[i].faces[j].vertex[k].particules.push_back(l);
	    //cout << i << " " << j << " " << k << " " <<  P[i].faces[j].vertex[k].num << " " << l << endl;
	    //getchar();
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
  
  //Initialisation de la position et de l'inertie du solide
  for(int i=0; i<solide.size(); i++){
    solide[i].Dx = Vector_3(0.,0.,0.);
    solide[i].Dxprev = Vector_3(0.,0.,0.);
    solide[i].Fi = Vector_3(0.,0.,0.);
    solide[i].Ff = Vector_3(0.,0.,0.);
    solide[i].Ffprev = Vector_3(0.,0.,0.);
    solide[i].Mi = Vector_3(0.,0.,0.);
    solide[i].Mf = Vector_3(0.,0.,0.);
    solide[i].Mfprev = Vector_3(0.,0.,0.);
    solide[i].e = Vector_3(0.,0.,0.);
    solide[i].eprev = Vector_3(0.,0.,0.);
    solide[i].Inertie();
    solide[i].mvt_t = Aff_transformation_3(1,0,0,0,1,0,0,0,1);
    solide[i].mvt_tprev = Aff_transformation_3(1,0,0,0,1,0,0,0,1);
  }

  //En cas de reprise
  if(rep){
    std::ostringstream oss;
    oss << "resultats/solide" << numrep << ".vtk";
    string s = oss.str();
    const char* nom = s.c_str();
    std::ifstream init(nom,std::ios::in);
    string dump;
    int nb_part = solide.size();
    int nb_triangles = 0.;
    for(int it=0; it<nb_part; it++){
      nb_triangles += solide[it].triangles.size();
    }
    for(int it=0;it<5*nb_triangles+13;it++){
      getline(init,dump);
    }
    cout << dump << endl;
    //Recuperation du deplacement
    for(int it=0; it<nb_part; it++){
      double Dx,Dy,Dz;
      init >> Dx >> Dy >> Dz;
      solide[it].Dx = Vector_3(Dx,Dy,Dz);
      for(int l= 0; l<solide[it].triangles.size(); l++){
	getline(init,dump);
      }
    }
    getline(init,dump);
    getline(init,dump);
    cout << dump << endl;
    //Recuperation de la vitesse
    for(int it=0; it<nb_part; it++){
      double u,v,w;
      init >> u >> v >> w;
      solide[it].u = Vector_3(u,v,w);
      for(int l= 0; l<solide[it].triangles.size(); l++){
	getline(init,dump);
      }
    }
    getline(init,dump);
    getline(init,dump);
    cout << dump << endl;
    //Recuperation du vecteur de rotation
    for(int it=0; it<nb_part; it++){
      double ex,ey,ez;
      init >> ex >> ey >> ez;
      solide[it].e = Vector_3(ex,ey,ez);
      for(int l= 0; l<solide[it].triangles.size(); l++){
	getline(init,dump);
      }
    }
    getline(init,dump);
    getline(init,dump);
    cout << dump << endl;
    //Recuperation de la vitesse de rotation
    for(int it=0; it<nb_part; it++){
      double omegax,omegay,omegaz;
      init >> omegax >> omegay >> omegaz;
      solide[it].omega = Vector_3(omegax,omegay,omegaz);
      for(int l= 0; l<solide[it].triangles.size(); l++){
	getline(init,dump);
      }
    }
    init.close();
    //Mise a jour de differents parametres
    for(int i=0; i<solide.size(); i++){
      solide[i].Dxprev = solide[i].Dx;
      solide[i].eprev = solide[i].e;
      double rot[3][3];
      //Recuperation de la matrice de rotation
      double e0 = sqrt(1.-CGAL::to_double(solide[i].e.squared_length()));
      rot[0][0] = 1.-2.*CGAL::to_double(solide[i].e.operator[](1)*solide[i].e.operator[](1)+solide[i].e.operator[](2)*solide[i].e.operator[](2));
      rot[0][1] = 2.*CGAL::to_double(-e0*solide[i].e.operator[](2)+solide[i].e.operator[](0)*solide[i].e.operator[](1));
      rot[0][2] = 2.*CGAL::to_double(e0*solide[i].e.operator[](1)+solide[i].e.operator[](0)*solide[i].e.operator[](2));
      rot[1][0] = 2.*CGAL::to_double(e0*solide[i].e.operator[](2)+solide[i].e.operator[](1)*solide[i].e.operator[](0));
      rot[1][1] = 1.-2.*CGAL::to_double(solide[i].e.operator[](0)*solide[i].e.operator[](0)+solide[i].e.operator[](2)*solide[i].e.operator[](2));
      rot[1][2] = 2.*CGAL::to_double(-e0*solide[i].e.operator[](0)+solide[i].e.operator[](1)*solide[i].e.operator[](2));
      rot[2][0] = 2.*CGAL::to_double(-e0*solide[i].e.operator[](1)+solide[i].e.operator[](2)*solide[i].e.operator[](0));
      rot[2][1] = 2.*CGAL::to_double(e0*solide[i].e.operator[](0)+solide[i].e.operator[](2)*solide[i].e.operator[](1));
      rot[2][2] = 1.-2.*CGAL::to_double(solide[i].e.operator[](0)*solide[i].e.operator[](0)+solide[i].e.operator[](1)*solide[i].e.operator[](1));
      Aff_transformation_3 rotation(rot[0][0],rot[0][1],rot[0][2],rot[1][0],rot[1][1],rot[1][2],rot[2][0],rot[2][1],rot[2][2]);
      Aff_transformation_3 translation(CGAL::TRANSLATION,Vector_3(Point_3(0.,0.,0.),solide[i].x0)+solide[i].Dx);
      Aff_transformation_3 translation_inv(CGAL::TRANSLATION,Vector_3(solide[i].x0,Point_3(0.,0.,0.))-solide[i].Dx);
	  solide[i].mvt_tprev = solide[i].mvt_t;
      solide[i].mvt_t = translation*(rotation*translation_inv);
     
    }
    update_triangles();
  }
  
}


void Solide::solve_position(double dt){
  for(int i=0;i<size();i++){
    solide[i].solve_position(dt);
  }
  update_triangles();
	for(int i=0;i<size();i++){
		double x_min = solide[i].max_x, y_min=solide[i].max_y, z_min=solide[i].max_z, 
		       x_max =solide[i].max_x, y_max=solide[i].max_y, z_max =solide[i].max_z;
					 
		for(int j=0;j<solide[i].triangles.size();j++){
			
			x_min = min( x_min ,min(CGAL::to_double(solide[i].triangles[j].operator[](0).x()), min(CGAL::to_double(solide[i].triangles[j].operator[]                       (1).x()), CGAL::to_double(solide[i].triangles[j].operator[](2).x()) )) );
			
			y_min = min( y_min ,min(CGAL::to_double(solide[i].triangles[j].operator[](0).y()), min(CGAL::to_double(solide[i].triangles[j].operator[] (1).y()), CGAL::to_double(solide[i].triangles[j].operator[](2).y()) )) );
			
			z_min = min( z_min ,min(CGAL::to_double(solide[i].triangles[j].operator[](0).z()), min(CGAL::to_double(solide[i].triangles[j].operator[] (1).z()), CGAL::to_double(solide[i].triangles[j].operator[](2).z()) )) );
			
			x_max = max( x_max ,max(CGAL::to_double(solide[i].triangles[j].operator[](0).x()), max(CGAL::to_double(solide[i].triangles[j].operator[] (1).x()), CGAL::to_double(solide[i].triangles[j].operator[](2).x()) )) );
			
			y_max = max( y_max ,max(CGAL::to_double(solide[i].triangles[j].operator[](0).y()), max(CGAL::to_double(solide[i].triangles[j].operator[] (1).y()), CGAL::to_double(solide[i].triangles[j].operator[](2).y()) )) );
			
			z_max = max( z_max ,max(CGAL::to_double(solide[i].triangles[j].operator[](0).z()), max(CGAL::to_double(solide[i].triangles[j].operator[] (1).z()), CGAL::to_double(solide[i].triangles[j].operator[](2).z()) )) );
		}
		solide[i].min_x = x_min; solide[i].min_y = y_min; solide[i].min_z = z_min;
		solide[i].max_x = x_max; solide[i].max_y = y_max; solide[i].max_z = z_max;
	}
	
}

void Solide::solve_vitesse(double dt){
  for(int i=0;i<size();i++){
    solide[i].solve_vitesse(dt);
  }
}

// void Solide::Forces_fluide(double dt){
// 	for(int i=0;i<size();i++){
// 		solide[i].forces_fluide(dt);
// 	}
// }

void Solide::update_triangles(){
  for(int i=0;i<solide.size();i++){
    solide[i].triangles_prev = solide[i].triangles;
    solide[i].normales_prev = solide[i].normales;
    solide[i].fluide_prev = solide[i].fluide;
		for(int it=0;it<solide[i].triangles.size();it++){
			solide[i].Points_interface_prev[it] = solide[i].Points_interface[it];
			solide[i].Triangles_interface_prev[it] = solide[i].Triangles_interface[it];
			solide[i].Position_Triangles_interface_prev[it] = solide[i].Position_Triangles_interface[it];
			solide[i].Points_interface[it].erase(solide[i].Points_interface[it].begin(),solide[i].Points_interface[it].end());
			solide[i].Triangles_interface[it].erase(solide[i].Triangles_interface[it].begin(),solide[i].Triangles_interface[it].end());	solide[i].Position_Triangles_interface[it].erase(solide[i].Position_Triangles_interface[it].begin(),
                                                       solide[i].Position_Triangles_interface[it].end());
    }
    solide[i].triangles.erase(solide[i].triangles.begin(),solide[i].triangles.end());
    solide[i].normales.erase(solide[i].normales.begin(),solide[i].normales.end());
    solide[i].fluide.erase(solide[i].fluide.begin(),solide[i].fluide.end());
    //Calcul de la nouvelle position des triangles
    for(int f=0;f<solide[i].faces.size();f++){
      Point_3 s,r,v;
      vector<Point_3> si;
      si.push_back(solide[i].mvt_t.transform(solide[i].faces[f].centre));
      int j = solide[i].faces[f].voisin;
      if(j>=0){
	     si.push_back(solide[j].mvt_t.transform(solide[i].faces[f].centre));
      }
      s = centroid(si.begin(),si.end());
      for(int k=0;k<solide[i].faces[f].size();k++){
	int kp = (k+1)%(solide[i].faces[f].size());
	vector<Point_3> ri,vi;
	for(int part=0;part<solide[i].faces[f].vertex[k].size();part++){
	  int p = solide[i].faces[f].vertex[k].particules[part];
	  ri.push_back(solide[p].mvt_t.transform(solide[i].faces[f].vertex[k].pos));
	}
	r = centroid(ri.begin(),ri.end());
	for(int part=0;part<solide[i].faces[f].vertex[kp].size();part++){
	  int p = solide[i].faces[f].vertex[kp].particules[part];
	  vi.push_back(solide[p].mvt_t.transform(solide[i].faces[f].vertex[kp].pos));
	}
	v = centroid(vi.begin(),vi.end());
	Vector_3 vect0(s,r);
	Vector_3 vect1(s,v);
	Triangle_3 Tri(s,r,v);
	solide[i].triangles.push_back(Tri);
	Vector_3 normale = CGAL::cross_product(vect0,vect1);
	normale = normale*(1./sqrt(CGAL::to_double(normale.squared_length())));
	solide[i].normales.push_back(normale);
	if(solide[i].faces[f].voisin == -1){
	  solide[i].fluide.push_back(true);
	} else {
	  solide[i].fluide.push_back(false);
	}
      }
    }
  }
}

double Solide::Energie(){
  return Energie_cinetique()+Energie_potentielle();
}

double Solide::Energie_cinetique(){
  double E = 0.;
  for(int i=0;i<size();i++){
    double u2 = CGAL::to_double(solide[i].u.squared_length());
    E += 1./2.*solide[i].m*u2;
    //Calcul de -1/2*tr(D j(Q^T omega)) = 1/2*(I1*Omega1^2+I2*Omega2^2+I3*Omega3^2)
    double Q[3][3];
    double rot[3][3];
    double e0 = sqrt(1.-CGAL::to_double(solide[i].e.squared_length()));
    //Recuperation de la matrice de rotation
    rot[0][0] = 1.-2.*CGAL::to_double(solide[i].e.operator[](1)*solide[i].e.operator[](1)+solide[i].e.operator[](2)*solide[i].e.operator[](2));
    rot[0][1] = 2.*CGAL::to_double(-e0*solide[i].e.operator[](2)+solide[i].e.operator[](0)*solide[i].e.operator[](1));
    rot[0][2] = 2.*CGAL::to_double(e0*solide[i].e.operator[](1)+solide[i].e.operator[](0)*solide[i].e.operator[](2));
    rot[1][0] = 2.*CGAL::to_double(e0*solide[i].e.operator[](2)+solide[i].e.operator[](1)*solide[i].e.operator[](0));
    rot[1][1] = 1.-2.*CGAL::to_double(solide[i].e.operator[](0)*solide[i].e.operator[](0)+solide[i].e.operator[](2)*solide[i].e.operator[](2));
    rot[1][2] = 2.*CGAL::to_double(-e0*solide[i].e.operator[](0)+solide[i].e.operator[](1)*solide[i].e.operator[](2));
    rot[2][0] = 2.*CGAL::to_double(-e0*solide[i].e.operator[](1)+solide[i].e.operator[](2)*solide[i].e.operator[](0));
    rot[2][1] = 2.*CGAL::to_double(e0*solide[i].e.operator[](0)+solide[i].e.operator[](2)*solide[i].e.operator[](1));
    rot[2][2] = 1.-2.*CGAL::to_double(solide[i].e.operator[](0)*solide[i].e.operator[](0)+solide[i].e.operator[](1)*solide[i].e.operator[](1));
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	Q[j][k] = rot[j][0]*solide[i].rotref[0][k];
	Q[j][k] += rot[j][1]*solide[i].rotref[1][k];
	Q[j][k] += rot[j][2]*solide[i].rotref[2][k];
      }
    }  
    double Omega[3];
    Omega[0] = Omega[1] = Omega[2] = 0.;
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	Omega[j] += CGAL::to_double(solide[i].omega.operator[](k)*Q[k][j]);
      }
    }
    E += 1./2.*(solide[i].I[0]*Omega[0]*Omega[0]+solide[i].I[1]*Omega[1]*Omega[1]+solide[i].I[2]*Omega[2]*Omega[2]);
  }
  return E;
}

void Solide::Forces_internes(){
  //Initialisation
  for(int i=0;i<size();i++){
    solide[i].Fi = Vector_3(0.,0.,0.);
    solide[i].Mi = Vector_3(0.,0.,0.);
  }
  //Calcul de la d�formation volumique epsilon de chaque particule
  for(int i=0;i<size();i++){
    solide[i].Volume_libre();
    solide[i].epsilon = 0.;
    for(int j=0;j<solide[i].size();j++){
      if(solide[i].faces[j].voisin>=0){
	int part = solide[i].faces[j].voisin;
	Vector_3 Sn = 1./2.*cross_product(Vector_3(solide[i].faces[j].vertex[0].pos,solide[i].faces[j].vertex[1].pos),Vector_3(solide[i].faces[j].vertex[0].pos,solide[i].faces[j].vertex[2].pos));
	Point_3 c1 = solide[i].mvt_t.transform(solide[i].faces[j].centre);
	Point_3 c2 = solide[part].mvt_t.transform(solide[i].faces[j].centre);
	Vector_3 Delta_u(c1,c2);
	solide[i].epsilon += 1./2./(solide[i].V+N_dim*nu/(1.-2.*nu)*solide[i].Vl)*CGAL::to_double(Sn*Delta_u);
      }
    }
  }
  //Calcul des forces pour chaque particule
  for(int i=0;i<size();i++){
    for(int j=0;j<solide[i].size();j++){
      if(solide[i].faces[j].voisin>=0){
	int part = solide[i].faces[j].voisin;
	double S = 1./2.*sqrt(CGAL::to_double(cross_product(Vector_3(solide[i].faces[j].vertex[0].pos,solide[i].faces[j].vertex[1].pos),Vector_3(solide[i].faces[j].vertex[0].pos,solide[i].faces[j].vertex[2].pos)).squared_length()));
	Vector_3 X1X2(solide[i].mvt_t.transform(solide[i].x0),solide[part].mvt_t.transform(solide[part].x0));
	double DIJ = sqrt(CGAL::to_double(X1X2.squared_length()));
	Vector_3 nIJ = X1X2/DIJ;
	Point_3 c1 = solide[i].mvt_t.transform(solide[i].faces[j].centre);
	Point_3 c2 = solide[part].mvt_t.transform(solide[i].faces[j].centre);
	Vector_3 Delta_u(c1,c2);
	Vector_3 XC1(solide[i].x0,solide[i].faces[j].centre);
	Vector_3 XC2(solide[part].x0,solide[i].faces[j].centre);
	double alpha = sqrt(CGAL::to_double(XC1.squared_length()))/(solide[i].faces[j].D0);
	double epsilonIJ = alpha*solide[i].epsilon+(1.-alpha)*solide[part].epsilon;
	//Force de rappel elastique
	solide[i].Fi = solide[i].Fi + S/solide[i].faces[j].D0*E/(1.+nu)*Delta_u;
	//Force de deformation volumique
	solide[i].Fi = solide[i].Fi + S*E*nu/(1.+nu)/(1.-2.*nu)*epsilonIJ*(nIJ+Delta_u/DIJ-(Delta_u*nIJ)/DIJ*nIJ);
	//Moment des forces appliquees
	solide[i].Mi = solide[i].Mi + cross_product(solide[i].mvt_t.transform(XC1),S/solide[i].faces[j].D0*E/(1.+nu)*Delta_u);
	solide[i].Mi = solide[i].Mi + cross_product(solide[i].mvt_t.transform(XC1),S*E*nu/(1.+nu)/(1.-2.*nu)*epsilonIJ*(nIJ+Delta_u/DIJ-(Delta_u*nIJ)/DIJ*nIJ));
	//Moments de flexion/torsion
	double alphan = (1.+2.*nu)*E/4./(1.+nu)/S*(solide[i].faces[j].Is+solide[i].faces[j].It);
	double alphas = E/4./(1.+nu)/S*((3.+2.*nu)*solide[i].faces[j].Is-(1.+2.*nu)*solide[i].faces[j].It);
	double alphat = E/4./(1.+nu)/S*((3.+2.*nu)*solide[i].faces[j].It-(1.+2.*nu)*solide[i].faces[j].Is);
	solide[i].Mi = solide[i].Mi + S/solide[i].faces[j].D0*(alphan*cross_product(solide[i].mvt_t.transform(solide[i].faces[j].normale),solide[part].mvt_t.transform(solide[i].faces[j].normale))+alphas*cross_product(solide[i].mvt_t.transform(solide[i].faces[j].s),solide[part].mvt_t.transform(solide[i].faces[j].s))+alphat*cross_product(solide[i].mvt_t.transform(solide[i].faces[j].t),solide[part].mvt_t.transform(solide[i].faces[j].t)));
      }
    }
  }
}



double Solide::Energie_potentielle(){
  double Ep = 0.;
  //Calcul de la d�formation volumique epsilon de chaque particule
  for(int i=0;i<size();i++){
    solide[i].Volume_libre();
    solide[i].epsilon = 0.;
    for(int j=0;j<solide[i].size();j++){
      if(solide[i].faces[j].voisin>=0){
	int part = solide[i].faces[j].voisin;
	Vector_3 Sn = 1./2.*cross_product(Vector_3(solide[i].faces[j].vertex[0].pos,solide[i].faces[j].vertex[1].pos),Vector_3(solide[i].faces[j].vertex[0].pos,solide[i].faces[j].vertex[2].pos));
	Point_3 c1 = solide[i].mvt_t.transform(solide[i].faces[j].centre);
	Point_3 c2 = solide[part].mvt_t.transform(solide[i].faces[j].centre);
	Vector_3 Delta_u(c1,c2);
	solide[i].epsilon += 1./2./(solide[i].V+N_dim*nu/(1.-2.*nu)*solide[i].Vl)*CGAL::to_double(Sn*Delta_u);
      }
    }
    //Energie de deformation volumique de la particule
    Ep += E*nu/2./(1.+nu)/(1.-2.*nu)*(solide[i].V+N_dim*nu/(1.-2.*nu)*solide[i].Vl)*pow(solide[i].epsilon,2);
  }
  //Calcul de l'energie pour chaque lien
  for(int i=0;i<size();i++){
    for(int j=0;j<solide[i].size();j++){
      if(solide[i].faces[j].voisin>=0){
	int part = solide[i].faces[j].voisin;
	double S = 1./2.*sqrt(CGAL::to_double(cross_product(Vector_3(solide[i].faces[j].vertex[0].pos,solide[i].faces[j].vertex[1].pos),Vector_3(solide[i].faces[j].vertex[0].pos,solide[i].faces[j].vertex[2].pos)).squared_length()));
	Vector_3 X1X2(solide[i].mvt_t.transform(solide[i].x0),solide[part].mvt_t.transform(solide[part].x0));
	double DIJ = sqrt(CGAL::to_double(X1X2.squared_length()));
	Vector_3 nIJ = X1X2/DIJ;
	Point_3 c1 = solide[i].mvt_t.transform(solide[i].faces[j].centre);
	Point_3 c2 = solide[part].mvt_t.transform(solide[i].faces[j].centre);
	Vector_3 Delta_u(c1,c2);
	Vector_3 XC1(solide[i].x0,solide[i].faces[j].centre);
	Vector_3 XC2(solide[part].x0,solide[i].faces[j].centre);
	double alpha = sqrt(CGAL::to_double(XC1.squared_length()))/(solide[i].faces[j].D0);
	double epsilonIJ = alpha*solide[i].epsilon+(1.-alpha)*solide[part].epsilon;
	//Energie de rappel elastique
	Ep += 1./4.*S/solide[i].faces[j].D0*E/(1.+nu)*CGAL::to_double(Delta_u*Delta_u);
	//Moments de flexion/torsion
	double alphan = (1.+2.*nu)*E/4./(1.+nu)/S*(solide[i].faces[j].Is+solide[i].faces[j].It);
	double alphas = E/4./(1.+nu)/S*((3.+2.*nu)*solide[i].faces[j].Is-(1.+2.*nu)*solide[i].faces[j].It);
	double alphat = E/4./(1.+nu)/S*((3.+2.*nu)*solide[i].faces[j].It-(1.+2.*nu)*solide[i].faces[j].Is);
	Ep += S/2./solide[i].faces[j].D0*(alphan*(1.-CGAL::to_double(solide[i].mvt_t.transform(solide[i].faces[j].normale)*solide[part].mvt_t.transform(solide[i].faces[j].normale)))+alphas*(1.-CGAL::to_double(solide[i].mvt_t.transform(solide[i].faces[j].s)*solide[part].mvt_t.transform(solide[i].faces[j].s)))+alphat*(1.-CGAL::to_double(solide[i].mvt_t.transform(solide[i].faces[j].t)*solide[part].mvt_t.transform(solide[i].faces[j].t))));
      }
    }
  }
  return Ep;
}

double Solide::pas_temps(double t, double T){
  double dt = 10000.;
  //Restriction CFL sur la vitesse de rotation
  for(int i=0;i<size();i++){
    double dt1 = cfls/(sqrt(CGAL::to_double(solide[i].omega.squared_length()))+eps);
    dt = min(dt,dt1); 
  }
  //Restriction CFL li�e aux forces internes
  double cs = sqrt(E*(1.-nu)/rhos/(1.+nu)/(1.-2.*nu));
  for(int i=0;i<size();i++){
    for(int j=0;j<solide[i].size();j++){
      if(solide[i].faces[j].voisin>=0){
	        dt = min(dt,cfls*solide[i].faces[j].D0/cs);
      }
    }
  }
  dt = min(dt,T-t); 
  return dt;
}



#endif
