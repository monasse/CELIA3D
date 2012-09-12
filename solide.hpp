#include "intersections.hpp"

#ifndef SOLIDE_HPP
#define SOLIDE_HPP

class Particule
{

 public:
   
  Particule();
	
	Particule(const double x_min, const double y_min, const double z_min, 
				 const double x_max, const double y_max,const double z_max);

	Particule(const double x_min, const double y_min, const double z_min, 
				const double x_max, const double y_max,const double z_max, 
				const Point_3 s1, const Point_3 r1, const Point_3 t1, const Point_3 v1,
				const Point_3 s2, const Point_3 r2, const Point_3 t2, const Point_3 v2);
 ~Particule();
 
 bool box_inside_convex_polygon(const Particule& S, const Bbox& cell);  
 bool inside_convex_polygon(const Particule& S, const Point_3& P);  
 bool inside_box(const Bbox& cell, const Point_3& P);
  void Affiche();
  double volume(); 
  
	int nb_faces;
    
	double min_x;
	double min_y;
	double min_z;
	double max_x;
	double max_y;
	double max_z;
	bool cube;
	
	Point_3 centre[6];	
	Triangles triangles;
	Vector_3 normales[12];    
};  
  
class Solide
{
	
	public:
		
		Solide():solide(std::vector<Particule>(1)){}
		Solide(std::vector<Particule> & Part);
		~Solide();
		void Affiche();
		int size(){
			return solide.size();
		}
		void impression(int n);
		// private :
		std::vector<Particule> solide;
};

#endif
