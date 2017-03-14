//Copyright 2017 Laurent Monasse

/*
  This file is part of CELIA3D.
  
  CELIA3D is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  CELIA3D is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with CELIA3D.  If not, see <http://www.gnu.org/licenses/>.
*/

/*!
 *  \file 
 \authors Laurent Monasse and Maria Adela Puscas
 *  \brief Definition of solid classes.
 * Specific coupling members are preceded by a "warning" sign.
 */

#include "intersections.hpp"
#ifndef SOLIDE_HPP
#define SOLIDE_HPP

//! Vertex class
class Vertex 
{
public:
  Vertex();
  Vertex(const Point_3 p, std::vector<int> & parts);
  Vertex & operator=(const  Vertex &V);
  Point_3 pos; //!< Vertex coordinates
  int num;//!< Index of the vertex in the construction mesh
  
  int size(){
    return particules.size();
  }
  std::vector<int> particules; //!< List of particles sharing the vertex
};

//! Face class
class Face
{
public:
  Face();
  Face(std::vector<Vertex> & v, int part);
  Face(std::vector<Vertex> & v, int part, double dist);
  Face & operator=(const  Face &F); 
  int size(){
    return vertex.size();
  }
  void compFaceIntegrals(double &Fa, double &Fb, double &Fc, double &Faa, double &Fbb, double &Fcc, double &Faaa, double &Fbbb, double &Fccc, double &Faab, double &Fbbc, double &Fcca, double na,double nb, double nc, int a, int b, int c);
  void compProjectionIntegrals(double &P1, double &Pa, double &Pb, double &Paa, double &Pab, double &Pbb, double &Paaa, double &Paab, double &Pabb, double &Pbbb, int a, int b, int c);
  void Inertie();
  Point_3 centre; //!< Face center
  Vector_3 normale; //!< Face exterior normal
  double S; //Face area
  double Is; //!< First inertia moment of the face
  double It; //!< Second inertia moment of the face
  Vector_3 s; //!< Eigenvector associated with the first inertia moment of the face
  Vector_3 t; //!< Eigenvector associated with the second inertia moment of the face
  std::vector<Vertex> vertex; //!< List of the face vertices
  int voisin; //!< Index of the neighbouring particle. The index takes the value -1 if there is no solid neighbour (possibly fluid)
  double D0; //!< Equilibrium distance with the neighbouring particle
};

  
//! Particule class
class Particule
{

public:
   
  Particule();
  
  Particule(const double x_min, const double y_min, const double z_min, 
	    const double x_max, const double y_max,const double z_max);
  
  Particule(Point_3 c, const double x_min, const double y_min, const double z_min, 
	    const double x_max, const double y_max,const double z_max, 
	    std::vector<Face> & F);
  ~Particule();
  Particule & operator=(const Particule &P); 
  void Affiche();  
  double volume(); 
  void CompVolumeIntegrals(double &T1, double &Tx, double &Ty, double &Tz, double &Txx, double &Tyy, double &Tzz, double &Txy, double &Tyz, double &Tzx);
  void Inertie();
  void Volume_libre();
  void solve_position(double dt);
  void solve_vitesse(double dt);
  Vector_3 vitesse_parois(const Point_3& X_f);  
  Vector_3 vitesse_parois_prev(const Point_3& X_f);
  bool cube; //!< = true if the particle is a cube, false otherwise
  Bbox bbox; //!< Bounding box of the particle
  
  std::vector<Face> faces; //!< List of the particle's faces

  std::vector<Point_3> vertices;//!< List of the particle's vertices
  
  Triangles triangles; //!< Triangulation of the particle faces at time t
    
  Triangles triangles_prev; //!< Triangulation of the particle faces at time t-dt
    
  std::vector<Vector_3> normales; //!< Exterior normals to \a Particule.triangles
    
  std::vector<Vector_3> normales_prev; //!< Exterior normals to \a Particule.triangles_prev
    
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  std::vector<bool> vide; //!< =true if \a Particule.triangles is in contact with void
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  std::vector<bool> fluide; //!< =true if \a Particule.triangles is in contact with fluid
    
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  std::vector<bool> fluide_prev; //!< =true if \a Particule.triangles_prev is in contact with fluid
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  std::vector< std::vector<Point_3> > Points_interface; //!< Listof intersection points of \a Particule.triangles with the fluid grid at time t 
    
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  std::vector< std::vector<Point_3> > Points_interface_prev; //!< List of intersection points of \a Particule.triangles_prev with the fluid grid at time t-dt 
  /*! 
   * \warning  <b> Specific coupling parameter! </b>
   */
  std::vector< std::vector<Triangle_3> > Triangles_interface; //!< Triangulation of \a Particule.triangles at time t
    
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  std::vector< std::vector< std::vector<int> > > Position_Triangles_interface; //!< Index of the cell where \a Triangles_interface is located at time t
    
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  std::vector< std::vector<Triangle_3> > Triangles_interface_prev; //!< Triangulation of \a Particule.triangles_prev at time t-dt  
    
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  std::vector< std::vector<std::vector<int> > > Position_Triangles_interface_prev; //!< Index of the cell where \a Triangles_interface is located at time t-dt

  int fixe; //!< =1 if the particle is fixed, 0 otherwise
  double m; //!< Particle mass
  double V; //!< Particle volume
  double Vl; //!< Free Volume of the particle (for the computation of epsilon)
  double epsilon; //!< Volumetric deformation of the particle
  double I[3]; //!< Inertia matrix of the particle
  double rotref[3][3]; //!<Rotation matrix \f$ Q_0 \f$ such that the inertia matrix \f$ R \f$ in the reference frame can be written :\f$ R = Q_0 R_0 Q_0^{-1}\f$, with \f$R_0=diag(I_1,I_2,I_3)\f$.
  Point_3 x0; //!<Position of the particle center at t=0
  Vector_3 Dx; //!<Displacement of the particle center at time t
  Vector_3 Dxprev; //!<Displacement of the particle center at time t-dt
  Vector_3 Fi; //!<Solid internal forces
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  Vector_3 Ff; //!<Fluid forces applied on the solid between times t and t+dt/2
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  Vector_3 Ffprev; //!< Fluid forces applied on the solide between times t-dt/2 and t
  Vector_3 Mi; //!< Interior torques of the solid
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  Vector_3 Mf; //!< Fluid torques applied on the solid between times t and t+dt/2
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  Vector_3 Mfprev; //!< FLuid torques applied on the solid between times t-dt/2 and t
  Vector_3 u; //!< Particle velocity at time t
  Vector_3 u_half; //!< Particle velocity at time t-dt/2
  Vector_3 omega; //!< Angular velocity at time t
  Vector_3 omega_half;//!< Angular velocity at time t-dt/2
  Vector_3 e; //!<Rotation vector at time t
  Vector_3 eprev; //!<Rotation vector at time t-dt
  Aff_transformation_3 mvt_t; //!<Affine transformation associated with the rigid body movement of the particle at time t
  Aff_transformation_3 mvt_tprev; //!<Affine transformation associated with the rigid body movement of the particle at time t-dt
}; 

//! Solide class
class Solide
{
	
public:
  
  Solide();
  Solide(std::vector<Particule> & Part);
  ~Solide();
  Solide & operator=(const Solide &S); 
  void Affiche();  
  int size(){
    return solide.size();
  }
  void Impression(int n);
  void Init(const char* s);
  void Solve_position(double dt);
  void Solve_vitesse(double dt);
  void Forces_internes();
  void update_triangles();
  void breaking_criterion();
  double Energie();
  double Energie_potentielle();
  double Energie_cinetique();
  double pas_temps(double t, double T);
  // private :
  std::vector<Particule> solide; //!< Solid mesh
};

bool inside_box(const Bbox& cell, const Point_3& P);
bool box_inside_convex_polygon(const Particule& S, const Bbox& cell);  
bool inside_convex_polygon(const Particule& S, const Point_3& P);  
double Error(Solide& S1, Solide& S2);
void Copy_f_m(Solide& S1, Solide& S2);
bool box_inside_tetra(const Tetrahedron &tetra, const Bbox& cell);
bool inside_tetra(const Tetrahedron &tetra, const Point_3& P);
#endif
