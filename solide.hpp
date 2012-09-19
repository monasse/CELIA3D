#include "intersections.hpp"

#ifndef SOLIDE_HPP
#define SOLIDE_HPP

class Vertex 
{
public:
  Vertex();
  Vertex(const Point_3 p, std::vector<int> & parts);
  Point_3 pos;
  int num;// numero du point dans le maillage de construction
  int size(){
	return particules.size();
  }
  std::vector<int> particules;
};

class Face
{
public:
  Face();//:vertex(std::vector<Vertex>(1)){}
  Face(std::vector<Vertex> & v, int part);
  Face(std::vector<Vertex> & v, int part, double dist);
  int size(){
	return vertex.size();
  }
  Point_3 centre;
  Vector_3 normale;
  std::vector<Vertex> vertex;
  int voisin;
  double D0; //Distance a l'equilibre avec la particule voisine
};

  

class Particule
{

 public:
   
  Particule();//:faces(std::vector<Face>(1)){}
  
  Particule(const double x_min, const double y_min, const double z_min, 
			const double x_max, const double y_max,const double z_max);
  
  Particule(Point_3 c, const double x_min, const double y_min, const double z_min, 
			const double x_max, const double y_max,const double z_max, 
			std::vector<Face> & F);
  ~Particule();
  
  bool box_inside_convex_polygon(const Particule& S, const Bbox& cell);  
  bool inside_convex_polygon(const Particule& S, const Point_3& P);  
  bool inside_box(const Bbox& cell, const Point_3& P);
  void Affiche();
  double volume(); 
  
  double min_x;
  double min_y;
  double min_z;
  double max_x;
  double max_y;
  double max_z;
  bool cube;

  int size(){
	return faces.size();
  }
  std::vector<Face> faces;
  Triangles triangles;
  std::vector<Vector_3> normales;
  std::vector<bool> fluide;
  Point_3 x0; //Position du centre de la particule a t=0
  Point_3 x; //Position du centre de la particule a t
  Point_3 xprev; //Position du centre de la particule a t-dt
  
};  
  
class Solide
{
	
public:
  
  Solide();//:solide(std::vector<Particule>(1)){}
  Solide(std::vector<Particule> & Part);
  ~Solide();
  void Affiche();
  int size(){
	return solide.size();
  }
  void impression(int n);
  void init(const char* s);
  // private :
  std::vector<Particule> solide;
};

#endif
