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
  void compFaceIntegrals(double &Fa, double &Fb, double &Fc, double &Faa, double &Fbb, double &Fcc, double &Faaa, double &Fbbb, double &Fccc, double &Faab, double &Fbbc, double &Fcca, double na,double nb, double nc, int a, int b, int c);
  void compProjectionIntegrals(double &P1, double &Pa, double &Pb, double &Paa, double &Pab, double &Pbb, double &Paaa, double &Paab, double &Pabb, double &Pbbb, int a, int b, int c);
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
	void Affiche();  //fonction auxilaire utile pour les tests
  double volume(); 
  void CompVolumeIntegrals(double &T1, double &Tx, double &Ty, double &Tz, double &Txx, double &Tyy, double &Tzz, double &Txy, double &Tyz, double &Tzx);
  void Inertie();
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
	std::vector< std::vector<Point_3> > Points_interface;
	std::vector< std::vector<Triangle_3> > Triangles_interface;
  double m; //masse de la particule
  double I[3]; //Moments d'inertie de la particule
  double rotref[3][3]; //Matrice de rotation Q_0 telle que la matrice d'inertie R s'ecrit : R = Q_0*R_0*Q_0^-1, avec R_0=diag(I1,I2,I3)
  Point_3 x0; //Position du centre de la particule a t=0
  Point_3 Dx; //Deplacement du centre de la particule a t
  Point_3 Dxprev; //Deplacement du centre de la particule a t-dt
  Vector_3 Fi; //Forces interieures du solide
  Vector_3 Ff; //Forces fluides exercees sur le solide
  Vector_3 Mi; //Moments interieurs du solide
  Vector_3 Mf; //Moments fluides exerces sur le solide
  double rot[3][3]; //Matrice de rotation de la particule
  
};  
  
class Solide
{
	
public:
  
  Solide();//:solide(std::vector<Particule>(1)){}
  Solide(std::vector<Particule> & Part);
  ~Solide();
	void Affiche();  //fonction auxilaire utile pour les test
  int size(){
	return solide.size();
  }
  void impression(int n);
  void init(const char* s);
  // private :
  std::vector<Particule> solide;
};

#endif
