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
   \file
   \authors Laurent Monasse and Maria Adela Puscas
   \brief Definition of classes Cellule and Grille used in the resolution of the fluid.
  Specific coupling members are outlined with a "warning" sign.
 */

#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <cassert>
#include "parametres.hpp"
#include "solide.hpp"


#ifndef FLUIDE_HPP
#define FLUIDE_HPP

using namespace std;


//!\brief Definition of class Cellule 
class Cellule {

public :

 // Constructeur 
  Cellule();
  Cellule(double x, double y, double z);
  Cellule(double x, double y, double z, double dx, double dy, double dz);
  
  Cellule & operator=(const Cellule &cell);

  // Destructor
  ~Cellule();

  bool is_in_cell(double x,double y, double z);

  void Affiche ();  
	
  //!\brief (x,y,z) Position of the center of the cell.   
  double x;       
  double y;
  double z; 

  //!\brief (i,j,k) Index of the cell
  int i;
  int j;
  int k;

  //!\brief (dx, dy,dz) Size of the cell.
  double dx;       
  double dy;
  double dz;

  double rho;       //!< Density in the cell at time t.
  double rho1;      //!< Density in the cell at time t-dt.

  double u;        //!< Velocity in the cell in the x direction.
  double v;        //!< Velocity in the cell in the y direction.
  double w;        //!< Velocity in the cell in the z direction.

  double p;         //!< Pressure in the cell at time t.
  double p1;        //!< Pressure in the cell at time t-dt.
    

  double impx;      //!< Momentum in the x direction at time t.
  double impy;      //!< Momentum in the y direction at time t.
  double impz;      //!< Momentum in the z direction at time t.

  double impx0;      //!< Momentum in the x direction before flux computation.
  double impy0;      //!< Momentum in the y direction before flux computation.
  double impz0;      //!< Momentum in the z direction before flux computation.

  double rhoE;      //!< Energy density in the cell at time t.
  double rhoE0;     //!< Energy density in the cell before flux computation.
  double rho0;      //!< Density in the cell before flux computation.



  double lambda[5];     //!< Eigenvalues of the Euler system.
  double rp[5];         //!< Variables for the TVD flux limitor on the left. 
  double rm[5];        //!< Variables for the TVD flux limitor on the right. 
   
    /*! 
     * \warning  <b>  Specific coupling parameter ! </b>
     */
    double pdtx;     //!< Effective pressure in the x direction during the time-step. 
    /*! 
     * \warning  <b>  Specific coupling parameter ! </b>
     */
    double pdty;     //!< Effective pressure in the y direction during the time-step. 
    /*! 
     * \warning  <b>  Specific coupling parameter ! </b>
     */
    double pdtz;     //!< Effective pressure in the z direction during the time-step. 
    
    /*! 
     * \warning  <b>  Specific coupling parameter ! </b>
     */
  
  double Mrho;      //!< Density mixing term for small cut-cells.  
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double Mimpx;     //!< x-momentum mixing term for small cut-cells.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double Mimpy;     //!< y-momentum mixing term for small cut-cells.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double Mimpz;     //!< Momentum mixing term for small cut-cells.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double MrhoE;     //!< Energy mixing term for small cut-cells.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double cells;      //!< Volume of mixing cells if the cell is a target of the mixing procedure.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double alpha;      //!< Solid occupancy ratio in the cell at time t.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double alpha0;    //!< Solid occupancy ratio in the cell at time t-dt.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double kappai;    //!< Solid occupancy ratio on the face in the x-direction.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double kappaj; //!< Solid occupancy ratio on the face in the y direction.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double kappak; //!< Solid occupancy ratio on the face in the z direction.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double kappai0;    //!< Solid occupancy ratio on the face in the x direction at time t-dt.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double kappaj0; //!< Solid occupancy ratio on the face in the y direction at time t-dt.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double kappak0; //!< Solid occupancy ratio on the face in the z direction at time t-dt.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  int proche;         //!< Flag which is equal to 0 far from the interface and 1 near the interface.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  int proche1;       //!< Flag which is equal to 0 far from the interface and 1 near the interface at time t-dt.
  
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  bool vide; //!< Flag which indicates whether a cell is filled with void (true) or fluid (false)
  
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  double flux_modif[5]; //!< Modified fluxes for cut-cells.
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  double delta_w[5];    //!< Swept quantity.
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  double phi_x;        //!< Flux at the interface in the x direction.
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  double phi_y; //!< Flux at the interface in the y direction.
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  double phi_z; //!< Flux at the interface in the z direction.
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double phi_v;  //!< Flux at the interface. 

  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double cible_alpha;      //!< Mixing term for solid occupancy in the small cut-cells.  
  /*! 
   * \warning  <b>  Specific coupling parameter ! </b>
   */
  double cible_rho;      //!< Mixing term for density in the small cut-cells.  
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  double cible_impx;     //!< Mixing term for x-momentum in small cut-cells.
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  double cible_impy;     //!< Mixing term for y-momentum in small cut-cells.
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  double cible_impz;     //!< Mixing term for z-momentum in small cut-cells.
  /*! 
   * \warning  <b> Specific coupling parameter ! </b>
   */
  double cible_rhoE;     //!< Mixing term for energy in small cut-cells.
  
  double cible_i;      //!< Index of the target mixing cell. 
  double cible_j;
  double cible_k;
  
  double xi;      //!< Position of the effective face center. 
  double yj;
  double zk;
  
  double fluxi[5];     //!< Flux in the x direction.
  double fluxj[5];     //!< Flux in the y direction.
  double fluxk[5];     //!< Flux in the z direction.
  
  double S;                  //!< Physical entropy.
  double ve[5];             //!< Vector of entropic variables.

  double fex;               //!< Entropy flux in the x direction.
  double fey;               //!< Entropy flux in the y direction.
  double fez;               //!< Entropy flux in the z direction.
  
  
  double Qci[5];            //!< Entropy correction in the x direction.
  double Qcj[5];            //!< Entropy correction in the y direction.
  double Qck[5];            //!< Entropy correction in the z direction.
  
  double dtfxi[5];          //!< Flux in x multiplied by dt.
  double dtfyj[5];          //!< Flux in y multiplied by dt.
  double dtfzk[5];          //!< Flux in z multiplied by dt.
  
  double delw[5];           //!< Delta V in i+1/2.
  double delwnu[5];         //!< Factor for the limitor in i+1/2.
  
  double cf2[5];            //!< Higher-order corrections. 
  double cf3[5]; 
  double cf4[5]; 
  double cf5[5]; 
  double cf6[5]; 
  double cf7[5]; 
  double cf8[5]; 
  double cf9[5]; 
  double cf10[5]; 
  double cf11[5]; 


  double psic0[5];           //!< Centered higher-order corrections.
  double psic1[5]; 
  double psic2[5]; 
  double psic3[5]; 
  double psic4[5]; 

  double psid0[5];            //!< Decentered higher-order corrections.
  double psid1[5]; 
  double psid2[5]; 
  double psid3[5]; 
  double psid4[5]; 

  double vpr[5][5];            //!< Matrix of eigenvectors of the system.


  double psic0r[5];           //!< Centered higher-order corrections in the basis of eigenvectors. 
  double psic1r[5]; 
  double psic2r[5]; 
  double psic3r[5]; 
  double psic4r[5]; 

  double psid0r[5];            //!< Decentered higher-order corrections in the basis of eigenvectors. 
  double psid1r[5]; 
  double psid2r[5]; 
  double psid3r[5]; 
  double psid4r[5]; 

  double psid[5]; 
  double am[5];                 //!< Measure of monotonicity.
  double am1[5];  
  
  int ordre;                    //!< Order of the approximate flux in the cell.
  double co[11];               //!< Storage of order coefficients. 
  
  friend class Grille;         //!< Class Cellule is visible only through class Grille.
  
};

//! Definition of class Grille
class Grille
{

 public:

  Grille();
 Grille(int Nx0,int Ny0, int Nz0, double dx0, double x0, double dy0, double y0, double dz0, double z0);

  ~Grille();

  void affiche();  

  void affiche(string r);
  
  Cellule cellule(int i, int j, int k);
  Cellule in_cell(const Point_3& p);
  void in_cell(const Point_3& p, int &i, int& j, int& k, bool& interieur);
  void Init();

  double pas_temps(double t, double T);
  
  void BC();

  double Masse();
  double Impulsionx();
  double Impulsiony();
  double Impulsionz();
  double Energie();
 
  void Impression(int n);
  
  void solve_fluidx(const double dt);
  void solve_fluidy(const double dt);
  void solve_fluidz(const double dt);
  void melange(const double dt);
  void fnumx( const double sigma, double t);
  void fnumy(const double sigma, double t);
  void fnumz( const double sigma, double t);


  void corentx(double sigma);
  void corenty(double sigma);
  void corentz(double sigma);
  void Solve(const double dt, double t, int n, Solide& S); 

  void Forces_fluide(Solide& S, const double dt);
  void Modif_fnum(const double dt);  
  void Mixage(); 
  void Fill_cel(Solide& S);
  void swap_face(const Triangles& T3d_prev, const Triangles& T3d_n, const double dt,  Particule & P, double & volume_test);
  void swap_face_inexact(const Triangle_3& Tr_prev, const Triangle_3& Tr, const Triangles& T3d_prev, const Triangles& T3d_n, const double dt,  Particule & P, double & volume_test);
  void cells_intersection_face(int& in,int& jn,int& kn,int& in1,int& jn1,int& kn1, std::vector<Bbox>& box_cells, std::vector<Cellule>& Cells);
  void Swap_2d(const double dt, Solide& S);
  void Swap_3d(const double dt, Solide& S); 
  Cellule voisin_fluide(const Cellule &c,  bool &target); 
  Cellule voisin_mixt(const Cellule &c,  bool &target); 
  Cellule voisin(const Cellule &c); 
  Cellule cible(const Cellule &c, std::vector< std::vector<int> > & tab_cible );
  void Mixage_cible();
  bool Mixage_cible2(); 
  void Parois_particles(Solide& S,double dt);
  std::vector<Point_3> intersection(Triangle_3 t1, Triangle_3 t2);
//private :

  double x;          //!< Position of the origin of the fluid grid.
  double y;
  double z;
  double dx;          //!< Spatial discretization step.
  double dy;
  double dz;
  vector< vector< vector<Cellule > > > grille; //!< Fluid mesh.
 

};

#endif
