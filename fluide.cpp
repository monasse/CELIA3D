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
 *  \brief Definition of the methods for classes Cellule and Grille. 
 * Specific coupling procedures are indicated with a "warning" sign.
 */
#include <iostream> 
#include <stdio.h> 
#include <fstream> 
#include <math.h> 
#include "fluide.hpp"



#ifndef FLUIDE_CPP
#define FLUIDE_CPP 

using namespace std;  

inline double sign(const double x)
{
  return (x < 0.) ? -1. : 1. ;
}

//Definition of the methods for class Cellule

/*!\brief Default constructor. 
*/

Cellule::Cellule()
{

  x = y = z = 1.;
		
  i=j=k=0;
		
  dx = dy = dz = 1.;
    
  rho = rho1 = 1.;
    
  u = v = w = 0.;
    
  p = p1 = 1.;
		
  vide = false;
    
  pdtx = pdty = pdtz = 0.;
    
  impx = rho*u; impy = rho*v; impz = rho*w;
    
  rhoE=rho*u*u/2. + rho*v*v/2.  + rho*w*w/2. + p/(gam-1.);
    
  Mrho = Mimpx = Mimpy = Mimpz = MrhoE = 0.;
  cible_alpha = cible_rho = cible_impx = cible_impy = cible_impz = cible_rhoE = 0.;
  cible_i=cible_j=cible_k=0;
		
  rho0 = impx0 = impy0 = impz0 = rhoE0 = 0.;
    
  cells = alpha = alpha0 = 0.;
    
  kappai = kappaj = kappak = kappai0 = kappaj0 = kappak0 = 0.;
    
  proche = proche1 = 0;
    
  lambda[0] = lambda[1] = lambda[2] = lambda[3] = lambda[4] = 1.; 
    
  xi = yj = zk = 1.;
    
  phi_x = phi_y = phi_z = phi_v=0.;
    
  S = log(p) - gam*log(rho);
    
  ve[0] = (1.-gam)/p*rhoE-(S-gam-1.);
  ve[1] = (gam-1.)/p*impx;
  ve[2] = (gam-1.)/p*impy;
  ve[3] = (gam-1.)/p*impz;
  ve[4] = (1-gam)*rho/p;
    
  fex = -impx*S; fey = -impy*S; fez = -impz*S;
    
    
  for(int l=0;l<5;l++){
        
    rp[l] = rm[l] = 1.; 
        
    fluxi[l] = fluxj[l] = fluxk[l] = flux_modif[l] = delta_w[l] = 0.; 
        
    Qci[l] = Qcj[l] = Qck[l] = 0.; 
        
    delw[l] = delwnu[l] = 1.; 
        
    cf2[l] = cf3[l] = cf4[l] = cf5[l] = cf6[l] = cf7[l] = cf8[l] = cf9[l] = cf10[l] = cf11[l] = 0.;
        
        
    psic0[l] = psic1[l] = psic2[l] = psic3[l] = psic4[l] = 0.; 
        
    psid0[l] = psid1[l] = psid2[l] = psid3[l] = psid4[l] = 0.; 
        
    psic0r[l] = psic1r[l] = psic2r[l] = psic3r[l] = psic4r[l] = 0.; 
        
    psid0r[l] = psid1r[l] = psid2r[l] = psid3r[l] = psid4r[l] = 0.; 
        
    psid[l] = 0.; 
    am[l] = am1[l] = 0.;
        
    for(int m=0;m<5;m++){ 
      vpr[l][m] = 0.; 
    }
  } 
    
  ordre = ordremax; 
    
  for(int l=0; l<ordre;l++){ 
    co[l] = 1.; 
  } 
    
}

/*!\brief Overload of the constructor.
   \param (x0,y0,z0) Coordinates of the cell center.
*/
Cellule::Cellule(double x0, double y0, double z0)
{

  x=x0; y=y0; z=z0;
     
  dx = dy = dz = 1.;
		
  i=j=k=0;
		
  rho = rho1 = 1.;
    
  u = v = w = 0.;
    
  p = p1 = 1.;
		
  vide = false;
    
  pdtx = pdty = pdtz = 0.;
    
  impx = rho*u; impy = rho*v; impz = rho*w;
    
  rhoE=rho*u*u/2. + rho*v*v/2.  + rho*w*w/2. + p/(gam-1.);
    
  Mrho = Mimpx = Mimpy = Mimpz = MrhoE = 0.;
  cible_alpha = cible_rho = cible_impx = cible_impy = cible_impz = cible_rhoE = 0.;
  cible_i=cible_j=cible_k=0;
    
  rho0 = impx0 = impy0 = impz0 = rhoE0 = 0.;
    
  cells = alpha = alpha0 = 0.;
    
  kappai = kappaj = kappak = kappai0 = kappaj0 = kappak0 = 0.;
    
  proche = proche1 = 0;
    
  lambda[0] = lambda[1] = lambda[2] = lambda[3] = lambda[4] = 1.; 
    
  xi = yj = zk = 1.;
    
  phi_x = phi_y = phi_z = phi_v = 0.;
    
  S = log(p) - gam*log(rho);
    
  ve[0] = (1.-gam)/p*rhoE-(S-gam-1.);
  ve[1] = (gam-1.)/p*impx;
  ve[2] = (gam-1.)/p*impy;
  ve[3] = (gam-1.)/p*impz;
  ve[4] = (1-gam)*rho/p;
    
  fex = -impx*S; fey = -impy*S; fez = -impz*S;
    
    
  for(int l=0;l<5;l++){
        
    rp[l] = rm[l] = 1.; 
        
    fluxi[l] = fluxj[l] = fluxk[l] = flux_modif[l] = delta_w[l] = 0.; 
        
    Qci[l] = Qcj[l] = Qck[l] = 0.; 
        
    delw[l] = delwnu[l] = 1.; 
        
    cf2[l] = cf3[l] = cf4[l] = cf5[l] = cf6[l] = cf7[l] = cf8[l] = cf9[l] = cf10[l] = cf11[l] = 0.;
        
        
    psic0[l] = psic1[l] = psic2[l] = psic3[l] = psic4[l] = 0.; 
        
    psid0[l] = psid1[l] = psid2[l] = psid3[l] = psid4[l] = 0.; 
        
    psic0r[l] = psic1r[l] = psic2r[l] = psic3r[l] = psic4r[l] = 0.; 
        
    psid0r[l] = psid1r[l] = psid2r[l] = psid3r[l] = psid4r[l] = 0.; 
        
    psid[l] = 0.; 
    am[l] = am1[l] = 0.;
        
    for(int m=0;m<5;m++){ 
      vpr[l][m] = 0.; 
    }
  } 
    
  ordre = ordremax; 
    
  for(int l=0; l<ordre;l++){ 
    co[l] = 1.; 
  } 
    
}

/*!\brief Overload of the constructor.
   \param (x0,y0,z0) Coordinates of the cell center.
   \param (dx0,dy0,dz0) Size of the cell.
*/

Cellule::Cellule(double x0, double y0, double z0, double dx0, double dy0, double dz0)
{

  x = x0; y = y0; z = z0;
    
  dx = dx0 ; dy = dy0 ; dz = dz0 ;
    
  i=j=k=0;
		
  rho = rho1 = 1.;
    
  u = v = w = 0.;
    
  p = p1 = 1.;
		
  vide = false;
    
  pdtx = pdty = pdtz = 0.;
    
  impx = rho*u; impy = rho*v; impz = rho*w;
    
  rhoE=rho*u*u/2. + rho*v*v/2.  + rho*w*w/2. + p/(gam-1.);
    
  Mrho = Mimpx = Mimpy = Mimpz = MrhoE = 0.;
  cible_alpha = cible_rho = cible_impx = cible_impy = cible_impz = cible_rhoE = 0.;
  cible_i=cible_j=cible_k=0;
    
  rho0 = impx0 = impy0 = impz0 = rhoE0 = 0.;
    
  cells = alpha = alpha0 = 0.;
    
  kappai = kappaj = kappak = kappai0 = kappaj0 = kappak0 = 0.;
    
  proche = proche1 = 0;
    
  lambda[0] = lambda[1] = lambda[2] = lambda[3] = lambda[4] = 1.; 
    
  xi = yj = zk = 1.;
    
  phi_x = phi_y = phi_z = phi_v = 0.;
    
  S = log(p) - gam*log(rho);
    
  ve[0] = (1.-gam)/p*rhoE-(S-gam-1.);
  ve[1] = (gam-1.)/p*impx;
  ve[2] = (gam-1.)/p*impy;
  ve[3] = (gam-1.)/p*impz;
  ve[4] = (1.-gam)*rho/p;
    
  fex = -impx*S; fey = -impy*S; fez = -impz*S;
    
    
  for(int l=0;l<5;l++){
        
    rp[l] = rm[l] = 1.; 
        
    fluxi[l] = fluxj[l] = fluxk[l] = flux_modif[l] = delta_w[l] = 0.; 
        
    Qci[l] = Qcj[l] = Qck[l] = 0.; 
        
    delw[l] = delwnu[l] = 1.; 
        
    cf2[l] = cf3[l] = cf4[l] = cf5[l] = cf6[l] = cf7[l] = cf8[l] = cf9[l] = cf10[l] = cf11[l] = 0.;
        
        
    psic0[l] = psic1[l] = psic2[l] = psic3[l] = psic4[l] = 0.; 
        
    psid0[l] = psid1[l] = psid2[l] = psid3[l] = psid4[l] = 0.; 
        
    psic0r[l] = psic1r[l] = psic2r[l] = psic3r[l] = psic4r[l] = 0.; 
        
    psid0r[l] = psid1r[l] = psid2r[l] = psid3r[l] = psid4r[l] = 0.; 
        
    psid[l] = 0.; 
    am[l] = am1[l] = 0.;
        
    for(int m=0;m<5;m++){ 
      vpr[l][m] = 0.; 
    }
  } 
    
  ordre = ordremax; 
    
  for(int l=0; l<ordre;l++){ 
    co[l] = 1.; 
  } 
    
}

/*!\brief Operator = overload.
   \param c Cellule
   \return Cellule
*/

Cellule & Cellule:: operator=(const Cellule &c){
    
  assert(this != &c);	

  x = c.x; y = c.y; z = c.z;
		
  i = c.i; j = c.j; k = c.k;
    
  dx = c.dx ; dy = c.dy ; dz = c.dz ;
    
  rho = c.rho; rho1 = c.rho1;
    
  u=c.u; v=c.v; w=c.w;
    
  p=c.p; p1 = c.p1;
		
  vide = c.vide;
    
  pdtx = c.pdtx; pdty = c.pdty; pdtz = c.pdtz;
    
  impx = c.impx; impy = c.impy; impz = c.impz;  
    
  phi_x = c.phi_x; phi_y = c.phi_y; phi_z = c.phi_z; phi_v= c.phi_v;
    
  rhoE=c.rhoE;
    
  Mrho = c.Mrho; Mimpx = c.Mimpx; Mimpy = c.Mimpy; Mimpz = c.Mimpz; MrhoE = c.MrhoE;
		
  cible_alpha = c.cible_alpha ; cible_rho = c.cible_rho ; cible_impx = c.cible_impx ; cible_impy = c.cible_impy ; cible_impz = c.cible_impz; cible_rhoE = c. cible_rhoE; cible_i = c.cible_i ; cible_j = c.cible_j; cible_k = c.cible_k;
		
  rho0 = c.rho0; impx0 = c.impx0; impy0 = c.impy0; impz0 = c.impz0; rhoE0 = c.rhoE0;
    
  cells = c.cells; alpha = c.alpha; alpha0 = c.alpha0;
    
  kappai = c.kappai; kappaj = c.kappaj; kappak = c.kappak; 
  kappai0 = c.kappai0; kappaj0 = c.kappaj0; kappak0 = c.kappak0;
    
  proche = c.proche; proche1 = c.proche1;
    
  xi = c.xi; yj = c.yj; zk = c.zk;
    
  S = c.S;
    
  fex = c.fex; fey = c.fey; fez = c.fez;
    
  for(int l=0;l<5;l++){ 
        
    lambda[l] = c. lambda[l]; 
    rp[l] = c.rp[l]; rm[l] = c.rm[l]; 
        
    fluxi[l] = c.fluxi[l]; fluxj[l] = c.fluxj[l]; fluxk[l] = c.fluxk[l];  
    delta_w[l] =c.delta_w[l];
				 
    dtfxi[l] = c.dtfxi[l]; dtfyj[l] = c.dtfyj[l]; dtfzk[l] = c.dtfzk[l];
				 
    Qci[l] = c.Qci[l]; Qcj[l] = c.Qcj[l]; Qck[l] = c.Qck[l];
        
    delw[l] = c.delw[l]; delwnu[l] = c.delwnu[l];
        
    ve[l] = c.ve[l];  
        
    cf2[l] = c.cf2[l]; cf3[l] = c.cf3[l]; cf4[l] = c.cf4[l]; cf5[l] = c.cf5[l]; cf6[l] = c.cf6[l]; 
    cf7[l] = c.cf7[l]; cf8[l] = c.cf8[l]; cf9[l] = c.cf9[l]; cf10[l] = c.cf10[l]; cf11[l] = c.cf11[l]; 
        
    psic0[l] = c.psic0[l]; psic1[l] = c.psic1[l]; psic2[l] = c.psic2[l]; 
    psic3[l] = c.psic3[l]; psic4[l] = c.psic4[l]; 
        
    psid0[l] = c.psid0[l]; psid1[l] = c.psid1[l]; psid2[l] = c.psid2[l]; 
    psid3[l] = c.psid3[l]; psid4[l] = c.psid4[l]; 
        
        
        
        
    psic0r[l] = c.psic0r[l]; psic1r[l] = c.psic1r[l]; psic2r[l] = c.psic2r[l]; 
    psic3r[l] = c.psic3r[l]; psic4r[l] = c.psic4r[l]; 
        
    psid0r[l] = c.psid0r[l]; psid1r[l] = c.psid1r[l]; psid2r[l] = c.psid2r[l]; 
    psid3r[l] = c.psid3r[l]; psid4r[l] = c.psid4r[l]; 
        
    psid[l] = c.psid[l]; 
    am[l] = c.am[l]; am1[l] = c.am1[l]; 
        
    for(int m=0;m<5;m++){ 
      vpr[l][m] = c.vpr[l][m]; 
    }
        
  } 
    
  ordre = c.ordre; 
    
  for(int l=0; l<ordre;l++){ 
    co[l] = c.co[l]; 
  } 
    
  return *this; 
}

/*!\brief Destructor.
*/ 
//Destructor
Cellule::~Cellule(){
}

/*!\brief The function tests whether point (x0,y0,z0) is in the cell.
   \param (x0,y0,z0) coodinates of the point
   \return bool (true if (x0,y0,z0) is in the cell, false otherwise)
*/
bool Cellule :: is_in_cell(double x0,double y0, double z0)
{
  bool k = false;
    
  if( (( dx*(1./2.) - abs(x0-x) )>-1.*eps)  && (( dy*(1./2.) - abs(y0-y) )>-1.*eps) && (( dz*(1./2.) - abs(z0-z) )>-1.*eps) )
  { k= true;}
  return k;
}

/*!\brief Test display.
*/
void Cellule::Affiche (){
    
  cout<<" center "<< " x = "<< x<< " y = "<<y<< " z = "<<z<<endl;
  cout<< " p = "<< p<< " rho = "<<rho<<endl;
  cout<<" alpha ="<<alpha<<endl;
  cout<<" cell state: vide "<<vide<<endl;  
}


//Definition of the methods of class Grille 

/*!\brief Default constructor.
   \details Variable grille represents the fluid mesh, it is a 3D array of \a Cellule.
*/
Grille::Grille(): grille(Nx+2*marge, vector< vector<Cellule> >(Ny+2*marge, vector<Cellule>(Nz+2*marge)) ){
    
  x = X0; y = Y0; z = Z0;
  dx = deltax; dy = deltay; dz = deltaz;
    
  for(int i=0;i<Nx+2*marge;i++){
    for(int j=0;j<Ny+2*marge;j++){ 
      for(int k=0;k<Nz+2*marge;k++){ 
	grille[i][j][k] = Cellule(x+dx/2.+(i-marge)*dx,y+dy/2.+(j-marge)*dy, z+dz/2.+(k-marge)*dz,dx, dy, dz);

      }
    } 
  }
   
} 
/*!\brief Overload of the constructor.
   \param (x0, y0, z0) Position of the origin of the fluid domain
   \param (dx0,dy0,dz0) Spatial discretization step
   \param (Nx0, Ny0, Nz0) Number of fluid cells in the x, y et z directions.
*/
Grille::Grille(int Nx0, int Ny0, int Nz0, double dx0, double x0, double dy0,double y0, double dz0, double z0):grille
													      (Nx0+2*marge, vector< vector<Cellule> > (Ny0+2*marge, vector<Cellule>(Nz0+2*marge))){ 
     
  x = x0; y = y0; z = z0;
    
  dx = dx0; dy = dy0; dz = dz0;
    
  for(int i=0;i<Nx0+2*marge;i++){
    for(int j=0;j<Ny0+2*marge;j++){ 
      for(int k=0;k<Nz0+2*marge;k++){ 
	grille[i][j][k] =	Cellule(x+dx/2.+(i-marge)*dx,y+dy/2.+(j-marge)*dy, z+dz/2.+(k-marge)*dz, dx, dy, dz); 
								
      }
    } 
  }
    
}
/*!\brief Destructor.
*/ 
//Destructeur
Grille::~Grille(){
}

/*!\brief Test display
*/
void Grille::affiche()
{
  double vol=0.;
  int count=0.;
  double variation_delta_w_rho=0., variation_volume=0.;
  for(int i=0;i<Nx+2*marge;i++){
    for(int j=0;j<Ny+2*marge;j++){ 
      for(int k=0;k<Nz+2*marge;k++){
	if(grille[i][j][k].vide ){
	  count++;
	}
	vol +=(dx*dy*dz)*grille[i][j][k].alpha;
	variation_delta_w_rho +=grille[i][j][k].delta_w[0];
	variation_volume += (grille[i][j][k].alpha - grille[i][j][k].alpha0)*grille[i][j][k].rho1;
      }
    }
  }
  cout<<"number of void cells:= "<< count <<endl;
}
/*!\brief Test display
*/
void Grille:: affiche(string r)
{
  int s=0;
  double vol=0.;
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      Cellule& cb = grille[i][j][marge];
      for(int k=marge;k<Nz+marge;k++){ 
	s++;
	Cellule& c = grille[i][j][k];
	if(abs(c.w)>eps){
	  cout << r << " " << c.x << " " << c.y << " " << c.z << " w=" << c.w << endl;
	}
      }
    }
  }
}

/*!\brief Access cell (i,j, k).
   \param (i,j,k) index of the cell
   \return Cellule
*/
Cellule Grille::cellule(int i,int j, int k){ 
  return grille[i][j][k]; 
}
/*!\brief Access the cell containing point \a p.
   \param p a point
   \return Cellule (cell containing \a p)
*/
Cellule Grille::in_cell(const Point_3& p){
  int i,j,k;
  i = (int) (floor(CGAL::to_double((p.operator[](0)-x)/dx))+marge);
  j = (int) (floor(CGAL::to_double((p.operator[](1)-y)/dy))+marge);
  k = (int) (floor(CGAL::to_double((p.operator[](2)-z)/dz))+marge);
  if(i<0){
    i = 0;
  }
  if(i>Nx+2*marge-1){
    i = Nx+2*marge-1;
  }
  if(j<0){
    j = 0;
  }
  if(j>Ny+2*marge-1){
    j = Ny+2*marge-1;
  }
  if(k<0){
    k = 0;
  }
  if(k>Nz+2*marge-1){
    k = Nz+2*marge-1;
  }
  return cellule(i,j,k);
}
/*!\brief Access cell containing point \a p.
   \param p a point
   \param (i,j,k) index of the cell containing point \a p
   \param interieur = true if point \a p is inside the fluid domain, false otherwise 
   \return void
*/
void Grille::in_cell(const Point_3& p, int &i, int& j, int& k, bool& interieur){
	
  i = (int) (floor(CGAL::to_double((p.operator[](0)-x)/dx))+marge);
  j = (int) (floor(CGAL::to_double((p.operator[](1)-y)/dy))+marge);
  k = (int) (floor(CGAL::to_double((p.operator[](2)-z)/dz))+marge);


  if (i<0 || i>Nx+2*marge || j<0 || j>Ny+2*marge || k<0 || k>Nz+2*marge)
  {interieur= false ;}
  else {interieur = true;}
}
/*!\brief Definition of initial conditions. 
   \return void
*/
void Grille::Init(){
    
  for(int i=0;i<Nx+2*marge;i++){
    for(int j=0;j<Ny+2*marge;j++){
      for(int k=0;k<Nz+2*marge;k++){
	Cellule& c = grille[i][j][k];
	c.dx = deltax; c.dy = deltay; c.dz = deltaz;
	c.x = x+c.dx/2.+(i-marge)*c.dx;
	c.y = y+c.dy/2.+(j-marge)*c.dy;
	c.z = z+c.dz/2.+(k-marge)*c.dz;
	c.i= i; c.j=j; c.k=k;
	c.rho = Rho(c.x, c.y, c.z); 
	c.u = U(c.x,c.y, c.z);
	c.v = V(c.x,c.y, c.z);
	c.w = W(c.x,c.y, c.z);
	c.p   = P(c.x,c.y, c.z);
	c.impx = c.rho*c.u; c.impy = c.rho*c.v; c.impz = c.rho*c.w;
	c.rhoE= c.rho*c.u*c.u/2. + c.rho*c.v*c.v/2. + c.rho*c.w*c.w/2. + c.p/(gam-1.); 
	c.kappai = c.kappaj = c.kappak = c.alpha = 0.;
	if(c.rho >eps_vide && c.p>eps_vide){c.vide=false;}
	else {c.vide=true;}
      }
    }
  }         
}

/*!\brief Computation of the fluid time-step.
   \param T total simulation time
   \param t current simulation time
   \return double
*/
double Grille::pas_temps(double t, double T){ 
    
  double dt = 10000.;
  //CFL condition on the fluid
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      for(int k=marge;k<Nz+marge;k++){
	Cellule& c = grille[i][j][k];
	if(!c.vide){
	  double c2 = gam*c.p/c.rho;
	  double dt1 = cfl*min(c.dx/(sqrt(c2)+abs(c.u)),min(c.dy/(sqrt(c2)+abs(c.v)), c.dz/(sqrt(c2)+abs(c.w))));
	  dt = min(dt,dt1); 
	}
      }
    } 
  }
  dt = min(dt,T-t); 
  return dt; 
}
/*!\brief Computation of the total fluid mass.
   \return  double
*/
double Grille::Masse(){ 
  double m = 0.; 
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      for(int k=marge;k<Nz+marge;k++){
	Cellule& c = grille[i][j][k]; 
	m += c.rho*c.dx*c.dy*c.dz*(1.-c.alpha);
      }
    }
  } 
  return m; 
} 
/*!\brief Computation of the total fluid x-momentum.
   \return double
*/
double Grille::Impulsionx(){ 
  double impx = 0.; 
  //Cellule c ;
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      for(int k=marge;k<Nz+marge;k++){
	Cellule& c = grille[i][j][k]; 
	impx += c.impx*c.dx*c.dy*c.dz*(1.-c.alpha);
      }
    }
  } 
  return impx; 
}
/*!\brief Computation of the total fluid y-momentum.
   \return double
*/
double Grille::Impulsiony(){ 
  double impy = 0.; 
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      for(int k=marge;k<Nz+marge;k++){
	Cellule& c = grille[i][j][k];  
	impy += c.impy*c.dx*c.dy*c.dz*(1.-c.alpha);
      }
    }
  } 
  return impy; 
}
/*!\brief Computation of the total fluid z-momentum.
   \return double
*/
double Grille::Impulsionz(){ 
  double impz = 0.; 
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      for(int k=marge;k<Nz+marge;k++){
	Cellule& c = grille[i][j][k];  
	impz += c.impz*c.dx*c.dy*c.dz*(1.-c.alpha);
      }
    }
  } 
  return impz; 
}
/*!\brief Computation of the total fluid energy. 
   \return double 
*/
double Grille::Energie(){ 
  double E = 0.; 
  for(int i=marge;i<Nx+marge;i++){ 
    for(int j=marge;j<Ny+marge;j++){ 
      for(int k=marge;k<Nz+marge;k++){
	Cellule&  c = grille[i][j][k]; 
	E += c.rhoE*c.dx*c.dy*c.dz*(1.-c.alpha); 
      }
    } 
  }
  return E; 
} 
/*!\brief Mixing the cells with negative density or pressure with neighbours. 
   \param dt time-step
   \return void
*/
void Grille::melange(const double dt){ 
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      for(int k=marge;k<Nz+marge;k++){
	Cellule& c = grille[i][j][k];
	if(c.rho<0. && !c.vide){
	  c.rho = c.rho0;
	  c.u = c.impx0/c.rho;
	  c.v = c.impy0/c.rho;
	  c.w = c.impz0/c.rho;
	  c.p = (gam-1)*(c.rhoE0-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2.-c.rho*c.w*c.w/2.);
	  c.impx = c.impx0;
	  c.impy = c.impy0;
	  c.impz = c.impz0;
	  c.rhoE = c.rhoE0;
	}
	if(c.p<0. && !c.vide){
	  c.rho = c.rho0;
	  c.u = c.impx0/c.rho;
	  c.v = c.impy0/c.rho;
	  c.w = c.impz0/c.rho;
	  c.p = (gam-1)*(c.rhoE0-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2.-c.rho*c.w*c.w/2.);
	  c.impx = c.impx0;
	  c.impy = c.impy0;
	  c.impz = c.impz0;
	  c.rhoE = c.rhoE0;
	}
      }
    }
  }
}

/*!\brief Solve the fluid equations in the x direction.
   \warning Storage of the pressure used during the time-step \a pdtx . Specific coupling parameter !
   \param dt time-step
   \return void
*/
void Grille::solve_fluidx(const double dt){ 
    
  const double sigma = dt/dx; 
  double dw1 =0., dw2=0., dw3=0., dw4=0., dw5 = 0.;
    
  //Computation of variables at time t+dt
  for(int i=1;i<Nx+2*marge-1;i++){ 
    for(int j=1;j<Ny+2*marge-1;j++){ 
      for(int k=1;k<Nz+2*marge-1;k++){
                
	Cellule& c = grille[i][j][k]; 
	Cellule& ci = grille[i-1][j][k];    
                
	//Storage of pressure used during the time-step
	c.pdtx = dt*c.p;
                
	dw1 = -sigma*(c.fluxi[0]-ci.fluxi[0]); 
	dw2 = -sigma*(c.fluxi[1]-ci.fluxi[1]); 
	dw3 = -sigma*(c.fluxi[2]-ci.fluxi[2]); 
	dw4 = -sigma*(c.fluxi[3]-ci.fluxi[3]);
	dw5 = -sigma*(c.fluxi[4]-ci.fluxi[4]);
                
	for(int l=0; l<5; l++){
	  c.dtfxi[l] = sigma*c.fluxi[l];
	}

	c.rho  += dw1; c.impx += dw2; c.impy += dw3; c.impz += dw4; c.rhoE += dw5;
	if (std::abs(c.rho) > eps_vide){
	  c.u = c.impx/c.rho; c.v = c.impy/c.rho; c.w = c.impz/c.rho;
	  c.p = (gam-1.)*(c.rhoE-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2.-c.rho*c.w*c.w/2.); 
	}
	if( (c.rho <= eps_vide && c.rho >= 0.) || (c.p <= eps_vide && c.p >= 0.) ){
	  c.vide = true;
	  c.u = 0.; c.v = 0.; c.w = 0.;
	  c.p = 0.; 
	}
	else {c.vide = false;}
      }
    }
  }
  melange(dt);   
} 



/*!\brief Resolution of fluid equations in the y direction.
   \warning Storage of the pressure used during the time-step \a pdty . Specific coupling parameter !
   \param dt time-step
   \return void
*/
void Grille::solve_fluidy(const double dt){ 
  const double sigma = dt/dy; 
  double dw1 =0., dw2=0., dw3=0., dw4=0., dw5 = 0.;
    
  //Computation of variables at time t+dt 
  for(int i=1;i<Nx+2*marge-1;i++){ 
    for(int j=1;j<Ny+2*marge-1;j++){ 
      for(int k=1;k<Nz+2*marge-1;k++){ 
                
	Cellule& c = grille[i][j][k]; 
	Cellule& cj = grille[i][j-1][k];    
                
                
	//Storage of the pressure used during the time-step
	c.pdty = dt*c.p;
                
	dw1 = -sigma*(c.fluxj[0]-cj.fluxj[0]); 
	dw2 = -sigma*(c.fluxj[1]-cj.fluxj[1]); 
	dw3 = -sigma*(c.fluxj[2]-cj.fluxj[2]); 
	dw4 = -sigma*(c.fluxj[3]-cj.fluxj[3]);
	dw5 = -sigma*(c.fluxj[4]-cj.fluxj[4]);
                
	for(int l=0;l<5;l++){
	  c.dtfyj[l] = sigma*c.fluxj[l];
	}
                
	c.rho  += dw1; c.impx += dw2; c.impy += dw3; c.impz += dw4; c.rhoE += dw5;
                
	if (std::abs(c.rho) > eps_vide){
	  c.u = c.impx/c.rho; c.v = c.impy/c.rho; c.w = c.impz/c.rho;
	  c.p = (gam-1.)*(c.rhoE-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2.-c.rho*c.w*c.w/2.); 
	}
								
	if( (c.rho <= eps_vide && c.rho >= 0.) || (c.p <= eps_vide && c.p >= 0.)){
	  c.vide = true;
	  c.u = 0.; c.v = 0.; c.w = 0.;
	  c.p = 0.; 
	}
	else {c.vide = false;}
      }
    }
  }
  melange(dt);    
} 

/*!\brief Resolution of the fluid equations in the z direction.
   \warning Storage of the pressure used during the time-step \a pdtz . Specific coupling parameter !
   \param dt time-step
   \return void
*/
void Grille::solve_fluidz(const double dt){
    
  const double sigma = dt/dz;
  double dw1 =0., dw2=0., dw3=0., dw4=0., dw5 = 0.;
    
  //Computation of variables at time t+dt
    
  for(int i=1;i<Nx+2*marge-1;i++){
    for(int j=1;j<Ny+2*marge-1;j++){
      for(int k=1;k<Nz+2*marge-1;k++){
                
	Cellule& c = grille[i][j][k]; 
	Cellule& ck = grille[i][j][k-1]; 
                
                
	//Storage of the pressure used during the time-step
	c.pdtz = dt*c.p;
                
	dw1 = -sigma*(c.fluxk[0]-ck.fluxk[0]); 
	dw2 = -sigma*(c.fluxk[1]-ck.fluxk[1]); 
	dw3 = -sigma*(c.fluxk[2]-ck.fluxk[2]); 
	dw4 = -sigma*(c.fluxk[3]-ck.fluxk[3]);
	dw5 = -sigma*(c.fluxk[4]-ck.fluxk[4]);
                
	for(int l=0;l<5;l++){
	  c.dtfzk[l] = sigma*c.fluxk[l];
	}
	c.rho  += dw1; c.impx += dw2; c.impy += dw3; c.impz += dw4; c.rhoE += dw5;

	if (std::abs(c.rho) > eps_vide){
	  c.u = c.impx/c.rho; c.v = c.impy/c.rho; c.w = c.impz/c.rho;
	  c.p = (gam-1.)*(c.rhoE-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2.-c.rho*c.w*c.w/2.); 
	}

	if( (c.rho <= eps_vide && c.rho >= 0.) || (c.p <= eps_vide && c.p >= 0.)){
	  c.u = 0.; c.v = 0.; c.w = 0.;
	  c.p = 0.; 
	  c.vide = true;
											
	}
	else {c.vide = false;}
								
      }
    }
  }
  melange(dt);  
}

/*!\brief Entropy correction in the x direction.
   \param sigma = \a dt/dx : time-step/ spatial discretization step for the fluid in the x direction
   \return void
*/
void Grille::corentx(double sigma){
    
  //Initialization of variables
  for(int i=0;i<Nx+2*marge;i++){
    for(int j=0;j<Ny+2*marge;j++){
      for(int k=0;k<Nz+2*marge;k++){
	Cellule& cel = grille[i][j][k];
	if(!cel.vide){ 
	  cel.S = log(cel.p) - gam*log(cel.rho);
	  cel.ve[0] = (1.-gam)/cel.p*cel.rhoE-(cel.S-gam-1.);
	  cel.ve[1] = (gam-1.)/cel.p*cel.impx;
	  cel.ve[2] = (gam-1.)/cel.p*cel.impy;
	  cel.ve[3] = (gam-1.)/cel.p*cel.impz;
	  cel.ve[4] = (1.-gam)*cel.rho/cel.p;
	  cel.fex = -cel.impx*cel.S;
	  cel.fey = -cel.impy*cel.S;
	  cel.fez = -cel.impz*cel.S;
	  cel.Qci[0] = cel.Qci[1] = cel.Qci[2] = cel.Qci[3] = cel.Qci[4] = 0.;
	}
      }
    }
  }
  //Computation of the entropy corrector in the x direction
    
  double df0 = 0., df1=0., df2=0., df3=0., df4=0.;
  double F0 = 0., F1=0., F2=0., F3=0., F4=0.;
    
  for(int i=0;i<Nx+2*marge-1;i++){
    for(int j=0;j<Ny+2*marge;j++){ 
      for(int k=0;k<Nz+2*marge;k++){
	Cellule& c = grille[i][j][k];
	Cellule& cd = grille[i+1][j][k];
	if(!c.vide && !cd.vide){ 
	  double alpha = 0.;
	  //Computation of pe
	  double pe = (cd.ve[0]-c.ve[0])*(cd.rho-c.rho);
	  pe += (cd.ve[1]-c.ve[1])*(cd.impx-c.impx);
	  pe += (cd.ve[2]-c.ve[2])*(cd.impy-c.impy);
	  pe += (cd.ve[3]-c.ve[3])*(cd.impz-c.impz);
	  pe += (cd.ve[4]-c.ve[4])*(cd.rhoE-c.rhoE);
									
	  //Computation of flux differences between neighbouring cells
	  df0 = (cd.impx-c.impx);
	  df1 = (cd.rho*cd.u*cd.u+cd.p)-(c.rho*c.u*c.u+c.p);
	  df2 = (cd.rho*cd.u*cd.v)-(c.rho*c.u*c.v);
	  df3 = (cd.rho*cd.u*cd.w)-(c.rho*c.u*c.w);
	  df4 = (cd.rhoE*cd.u+cd.p*cd.u)-(c.rhoE*c.u+c.p*c.u);
									
	  //Computation of the centered flux
	  F0 = 1./2.*(cd.impx+c.impx);
	  F1 = 1./2.*((cd.rho*cd.u*cd.u+cd.p)+(c.rho*c.u*c.u+c.p));
	  F2 = 1./2.*((cd.rho*cd.u*cd.v)+(c.rho*c.u*c.v));
	  F3 = 1./2.*((cd.rho*cd.u*cd.w)+(c.rho*c.u*c.w));
	  F4 = 1./2.*((cd.rhoE*cd.u+cd.p*cd.u)+(c.rhoE*c.u+c.p*c.u));
									
	  //Computation of qef
	  double qef = cd.fex - c.fex;
	  qef -= 0.5*(cd.ve[0]+c.ve[0])*df0;
	  qef -= 0.5*(cd.ve[1]+c.ve[1])*df1;
	  qef -= 0.5*(cd.ve[2]+c.ve[2])*df2;
	  qef -= 0.5*(cd.ve[3]+c.ve[3])*df3;
	  qef -= 0.5*(cd.ve[4]+c.ve[4])*df4;
									
	  //Computation of q-q*
	  double qmqet = qef;
	  qmqet += (cd.ve[0]-c.ve[0])*(c.fluxi[0]-F0);
	  qmqet += (cd.ve[1]-c.ve[1])*(c.fluxi[1]-F1);
	  qmqet += (cd.ve[2]-c.ve[2])*(c.fluxi[2]-F2);
	  qmqet += (cd.ve[3]-c.ve[3])*(c.fluxi[3]-F3);
	  qmqet += (cd.ve[4]-c.ve[4])*(c.fluxi[4]-F4);
	  qmqet *= -2.*sigma;
									
	  //Computation of alpha
	  if(pe>eps){
	    alpha = 2.*max(qef,0.)/pe;
	  }
									
	  //Computation of the right entropy corrector
									
	  c.Qci[0] = alpha*(cd.rho-c.rho);
	  c.Qci[1] = alpha*(cd.impx-c.impx);
	  c.Qci[2] = alpha*(cd.impy-c.impy);
	  c.Qci[3] = alpha*(cd.impz-c.impz);
	  c.Qci[4] = alpha*(cd.rhoE-c.rhoE);
									
	  for(int l=0;l<5;l++){
	    c.fluxi[l] -= c.Qci[l]; //modification flux
	  }
	}
      }
    }
  }
}

/*!\brief Entropy correction in the y direction. 
   \param sigma = \a dt/dy : time-step/ spatial discretization step for the fluid in the y direction
   \return void
*/
void Grille::corenty(double sigma){
    
  //Initialization of variables
  for(int i=0;i<Nx+2*marge;i++){
    for(int j=0;j<Ny+2*marge;j++){ 
      for(int k=0;k<Nz+2*marge;k++){
	Cellule& cel = grille[i][j][k];
	if(!cel.vide){
	  cel.S = log(cel.p) - gam*log(cel.rho);
	  cel.ve[0] = (1.-gam)/cel.p*cel.rhoE-(cel.S-gam-1.);
	  cel.ve[1] = (gam-1.)/cel.p*cel.impx;
	  cel.ve[2] = (gam-1.)/cel.p*cel.impy;
	  cel.ve[3] = (gam-1.)/cel.p*cel.impz;
	  cel.ve[4] = (1.-gam)*cel.rho/cel.p;
	  cel.fex = -cel.impx*cel.S;
	  cel.fey = -cel.impy*cel.S;
	  cel.fez = -cel.impz*cel.S;
	  cel.Qcj[0] = cel.Qcj[1] = cel.Qcj[2] = cel.Qcj[3] = cel.Qcj[4] = 0.;
	}
      }
    }
  }
  //Computation of the entropy corrector in the y direction
    
  double df0 = 0., df1=0., df2=0., df3=0., df4=0.;
  double F0 = 0., F1=0., F2=0., F3=0., F4=0.;
    
  for(int i=0;i<Nx+2*marge;i++){
    for(int j=0;j<Ny+2*marge-1;j++){ 
      for(int k=0;k<Nz+2*marge;k++){
	Cellule& c = grille[i][j][k];
	Cellule& ch = grille[i][j+1][k];
	if(!c.vide && !ch.vide){
	  double alpha = 0.;
	  //Computation of pe
	  double pe = (ch.ve[0]-c.ve[0])*(ch.rho-c.rho);
	  pe += (ch.ve[1]-c.ve[1])*(ch.impx-c.impx);
	  pe += (ch.ve[2]-c.ve[2])*(ch.impy-c.impy);
	  pe += (ch.ve[3]-c.ve[3])*(ch.impz-c.impz);
	  pe += (ch.ve[4]-c.ve[4])*(ch.rhoE-c.rhoE);
									
	  //Computation of flux differences in neighbouring cells
	  df0 = (ch.impy-c.impy);
	  df1 = (ch.rho*ch.u*ch.v)-(c.rho*c.u*c.v);
	  df2 = (ch.rho*ch.v*ch.v+ch.p)-(c.rho*c.v*c.v+c.p);
	  df3 = (ch.rho*ch.v*ch.w)-(c.rho*c.v*c.w);
	  df4 = (ch.rhoE*ch.v+ch.p*ch.v)-(c.rhoE*c.v+c.p*c.v);
									
	  //Computation of the centered flux
	  F0 = 1./2.*(ch.impy+c.impy);
	  F1 = 1./2.*((ch.rho*ch.u*ch.v)+(c.rho*c.u*c.v));
	  F2 = 1./2.*((ch.rho*ch.v*ch.v+ch.p)+(c.rho*c.v*c.v+c.p));
	  F3 = 1./2.*((ch.rho*ch.v*ch.w)+(c.rho*c.v*c.w));
	  F4 = 1./2.*((ch.rhoE*ch.v+ch.p*ch.v)+(c.rhoE*c.v+c.p*c.v));
									
	  //Computation of qef
	  double qef = ch.fey - c.fey;
	  qef -= 0.5*(ch.ve[0]+c.ve[0])*df0;
	  qef -= 0.5*(ch.ve[1]+c.ve[1])*df1;
	  qef -= 0.5*(ch.ve[2]+c.ve[2])*df2;
	  qef -= 0.5*(ch.ve[3]+c.ve[3])*df3;
	  qef -= 0.5*(ch.ve[4]+c.ve[4])*df4;
									
	  //Computation of q-q*
	  double qmqet = qef;
	  qmqet += (ch.ve[0]-c.ve[0])*(c.fluxj[0]-F0);
	  qmqet += (ch.ve[1]-c.ve[1])*(c.fluxj[1]-F1);
	  qmqet += (ch.ve[2]-c.ve[2])*(c.fluxj[2]-F2);
	  qmqet += (ch.ve[3]-c.ve[3])*(c.fluxj[3]-F3);
	  qmqet += (ch.ve[3]-c.ve[3])*(c.fluxj[4]-F4);
	  qmqet *= -2.*sigma;
	  //Computation of alpha
	  if(pe>eps){
	    alpha = 2.*max(qef,0.)/pe;
	  }
	  //Computation of the entropy corrector
	  c.Qcj[0] = alpha*(ch.rho-c.rho);
	  c.Qcj[1] = alpha*(ch.impx-c.impx);
	  c.Qcj[2] = alpha*(ch.impy-c.impy);
	  c.Qcj[3] = alpha*(ch.impz-c.impz);
	  c.Qcj[4] = alpha*(ch.rhoE-c.rhoE);
									
	  for(int l=0;l<5;l++){
	    c.fluxj[l] -= c.Qcj[l]; //modification flux
	  }
	}
      }
    }
  }
    
}


/*!\brief Entropy correction in the z direction. 
   \param sigma = \a dt/dz : time-step/ spatial discretization step for the fluid in the z direction
   \return void
*/
void Grille::corentz(double sigma){
  //Initialization of variables
  for(int i=0;i<Nx+2*marge;i++){
    for(int j=0;j<Ny+2*marge;j++){ 
      for(int k=0;k<Nz+2*marge;k++){
	Cellule& cel = grille[i][j][k];
	if(!cel.vide){
	  cel.S = log(cel.p) - gam*log(cel.rho);
	  cel.ve[0] = (1.-gam)/cel.p*cel.rhoE-(cel.S-gam-1.);
	  cel.ve[1] = (gam-1.)/cel.p*cel.impx;
	  cel.ve[2] = (gam-1.)/cel.p*cel.impy;
	  cel.ve[3] = (gam-1.)/cel.p*cel.impz;
	  cel.ve[4] = (1.-gam)*cel.rho/cel.p;
	  cel.fex = -cel.impx*cel.S;
	  cel.fey = -cel.impy*cel.S;
	  cel.fez = -cel.impz*cel.S;
	  cel.Qck[0] = cel.Qck[1] = cel.Qck[2] = cel.Qck[3] = cel.Qck[4] = 0.;
	}
      }
    }
  }
  //Computation of the entropy corrector in the z direction
    
  double df0 = 0., df1=0., df2=0., df3=0., df4=0.;
  double F0 = 0., F1=0., F2=0., F3=0., F4=0.;
    
  for(int i=0;i<Nx+2*marge;i++){
    for(int j=0;j<Ny+2*marge;j++){ 
      for(int k=0;k<Nz+2*marge-1;k++){
	Cellule& c = grille[i][j][k];
	Cellule& ch = grille[i][j][k+1];
	if(!c.vide && !ch.vide){
	  double alpha = 0.;
	  //Computation of pe
	  double pe = (ch.ve[0]-c.ve[0])*(ch.rho-c.rho);
	  pe += (ch.ve[1]-c.ve[1])*(ch.impx-c.impx);
	  pe += (ch.ve[2]-c.ve[2])*(ch.impy-c.impy);
	  pe += (ch.ve[3]-c.ve[3])*(ch.impz-c.impz);
	  pe += (ch.ve[4]-c.ve[4])*(ch.rhoE-c.rhoE);
									
	  //Computation of flux differences in two neighbouring cells
	  df0 = (ch.impz-c.impz);
	  df1 = (ch.rho*ch.u*ch.w)-(c.rho*c.u*c.w);
	  df2 = (ch.rho*ch.v*ch.w)-(c.rho*c.v*c.w);
	  df3 = (ch.rho*ch.w*ch.w+ch.p)-(c.rho*c.w*c.w+c.p);
	  df4 = (ch.rhoE*ch.w+ch.p*ch.w)-(c.rhoE*c.w+c.p*c.w);
									
	  //Computation of centered flux
	  F0 = 1./2.*(ch.impz+c.impz);
	  F1 = 1./2.*((ch.rho*ch.u*ch.w)+(c.rho*c.u*c.w));
	  F3 = 1./2.*((ch.rho*ch.w*ch.w+ch.p)+(c.rho*c.w*c.w+c.p));
	  F2 = 1./2.*((ch.rho*ch.v*ch.w)+(c.rho*c.v*c.w));
	  F4 = 1./2.*((ch.rhoE*ch.w+ch.p*ch.w)+(c.rhoE*c.w+c.p*c.w));
									
	  //Computation of qef
	  double qef = ch.fez - c.fez;
	  qef -= 0.5*(ch.ve[0]+c.ve[0])*df0;
	  qef -= 0.5*(ch.ve[1]+c.ve[1])*df1;
	  qef -= 0.5*(ch.ve[2]+c.ve[2])*df2;
	  qef -= 0.5*(ch.ve[3]+c.ve[3])*df3;
	  qef -= 0.5*(ch.ve[4]+c.ve[4])*df4;
									
	  //Computation of q-q*
	  double qmqet = qef;
	  qmqet += (ch.ve[0]-c.ve[0])*(c.fluxk[0]-F0);
	  qmqet += (ch.ve[1]-c.ve[1])*(c.fluxk[1]-F1);
	  qmqet += (ch.ve[2]-c.ve[2])*(c.fluxk[2]-F2);
	  qmqet += (ch.ve[3]-c.ve[3])*(c.fluxk[3]-F3);
	  qmqet += (ch.ve[3]-c.ve[3])*(c.fluxk[4]-F4);
	  qmqet *= -2.*sigma;
	  //Computation of alpha
	  if(pe>eps){
	    alpha = 2.*max(qef,0.)/pe;
	  }
	  //Computation of the entropy corrector
	  c.Qck[0] = alpha*(ch.rho-c.rho);
	  c.Qck[1] = alpha*(ch.impx-c.impx);
	  c.Qck[2] = alpha*(ch.impy-c.impy);
	  c.Qck[3] = alpha*(ch.impz-c.impz);
	  c.Qck[4] = alpha*(ch.rhoE-c.rhoE);
									
	  for(int l=0;l<5;l++){
	    c.fluxk[l] -= c.Qck[l]; //Modification of fluxes
	  }	
	}
      }
    }
  }
    
    
    
}  

/*!\brief Computation of the numerical x-flux.
   \param sigma = \a dt/dx: time-step/fluid spatial discretization step in the x direction
   \param t current simulation time
   \return void
*/ 
void Grille::fnumx(const double sigma, double t){
  //Initialization to the centered flux
  for(int i=0; i<Nx+2*marge-1; i++){ 
    for(int j=0; j<Ny+2*marge-1; j++){
      for(int k=0; k<Nz+2*marge-1; k++){

	Cellule& c = grille[i][j][k];  
	Cellule& ci = grille[i+1][j][k]; 
	if(!c.vide && !ci.vide){  
	  //Order indicators
	  for(int l=0;l< c.ordre;l++){ 
	    c.co[l]=1.; 
	  } 
	  for(int l=c.ordre;l<ordremax;l++){ 
	    c.co[l]=0.;
	  } 
                
	  //Centered flux part
	  c.fluxi[0] = (c.impx+ci.impx)/2.; 
	  c.fluxi[1] = (c.rho*c.u*c.u+c.p+ci.rho*ci.u*ci.u+ci.p)/2.;
	  c.fluxi[2] = (c.rho*c.u*c.v+ci.rho*ci.u*ci.v)/2.;
	  c.fluxi[3] = (c.rho*c.u*c.w+ci.rho*ci.u*ci.w)/2.;
	  c.fluxi[4] = ((c.rhoE+c.p)*c.u+(ci.rhoE+ci.p)*ci.u)/2.; 
	}
	else if(std::abs(c.alpha-1.)>eps){
	  //Lax-Friedrichs flux near void cells
	  c.fluxi[0] = (c.impx+ci.impx)/2. + (c.rho -ci.rho)/2./sigma; 
	  c.fluxi[1] = (c.rho*c.u*c.u+c.p+ci.rho*ci.u*ci.u+ci.p)/2. + (c.impx -ci.impx)/2./sigma; 
	  c.fluxi[2] = (c.rho*c.u*c.v+ci.rho*ci.u*ci.v)/2. + (c.impy -ci.impy)/2./sigma; 
	  c.fluxi[3] = (c.rho*c.u*c.w+ci.rho*ci.u*ci.w)/2. + (c.impz -ci.impz)/2./sigma; 
	  c.fluxi[4] = ((c.rhoE+c.p)*c.u+(ci.rhoE+ci.p)*ci.u)/2. + (c.rhoE -ci.rhoE)/2./sigma; 
	}
	else{
	  c.fluxi[0] = 0.; 
	  c.fluxi[1] = 0.; 
	  c.fluxi[2] = 0.; 
	  c.fluxi[3] = 0.; 
	  c.fluxi[4] = 0.; 
	}
      }
    }
  } 
    
  //Loop on the cells: computation of preliminary limiter variables
  for(int i=0; i<Nx+2*marge-1; i++){
    for(int j=0; j<Ny+2*marge-1; j++){
      for(int k=0; k<Nz+2*marge-1; k++){
	Cellule& ci = grille[i+1][j][k]; 
	Cellule& c = grille[i][j][k]; 
	if(!c.vide && !ci.vide){ 
	  //Computation of the Roe variables 
	  double roe = sqrt(ci.rho/c.rho); 
	  double rhor = roe*c.rho; 
	  double ur = (roe*ci.u+c.u)/(1.+roe); 
	  double vr = (roe*ci.v+c.v)/(1.+roe); 
	  double wr = (roe*ci.w+c.w)/(1.+roe); 
	  double Hr = (roe*(ci.rho*ci.u*ci.u/2.+ ci.rho*ci.v*ci.v/2.+ ci.rho*ci.w*ci.w/2. + ci.p*gam/(gam-1.))/ci.rho + (c.rho*c.u*c.u/2.+c.rho*c.v*c.v/2.+ c.rho*c.w*c.w/2. + c.p*gam/(gam-1.))/c.rho)/(1.+roe);
	  double ur2 = ur*ur;
	  double vr2 = vr*vr;
	  double wr2 = wr*wr;
	  double cr2 = (gam-1.)*(Hr-ur2/2.-vr2/2.-wr2/2.); 
										
	  //Test on the sound velocity 
	  if(cr2<=0. && abs(c.alpha-1.)>eps){
	    cout << "x-flux computation" << endl;
	    cout << "i=" << i << " j=" << j << " k=" << k<< " negative speed of sound : c2=" << cr2 << endl;
	    cout << "x=" << c.x << " y=" << c.y << " z=" << c.z<< " alpha=" << c.alpha << endl;
	    cout << "t=" << t << endl; 
	    cout << "c.p=" << c.p << endl; 
	    cout << "c.rho=" << c.rho << endl; 
	    cout << "c.u=" << c.u << endl; 
	    cout << "c.v=" << c.v << endl; 
	    cout << "c.w=" << c.w << endl; 
	    cout << "ci.p=" << ci.p << endl; 
	    cout << "ci.rho=" << ci.rho << endl; 
	    cout << "ci.u=" << ci.u << endl;
	    cout << "ci.v=" << ci.v << endl;
	    cout << "ci.w=" << ci.w << endl;
	    cout << "ur=" << ur << endl; 
	    cout << "ur2=" << ur2 << endl; 
	    cout << "vr=" << vr << endl; 
	    cout << "vr2=" << vr2 << endl; 
	    cout << "wr2=" << wr2 << endl; 
	    cout << "Hr=" << Hr << endl;
	    getchar();
	  } 
										
	  //Speed of sound
	  double cr = sqrt(cr2);
										
	  //Eigenvalues 
	  c.lambda[0] = ur-cr; 
	  c.lambda[1] = ur;
	  c.lambda[2] = ur;
	  c.lambda[3] = ur;
	  c.lambda[4] = ur+cr; 
										
	  //Computation of differences between Wd and Wg
	  double drho = ci.rho - c.rho; 
	  double du = ci.u - c.u;
	  double dv = ci.v - c.v;
	  double dw = ci.w - c.w;
	  double dp = ci.p - c.p; 
										
	  //Computation of the deltaV (differences between Wd and Wg in the eigenvectors basis) 
	  double ros2c = rhor/cr/2.; 
	  c.delw[0] = dp/cr2/2. - ros2c*du; 
	  c.delw[1] = drho - dp/cr2;
	  c.delw[2] = 2.*ros2c*dv;
	  c.delw[3] = 2.*ros2c*dw;   
	  c.delw[4] = dp/cr2/2. + ros2c*du; 
										
	  //Computation of the complete correction in the eigenvectors basis
	  double xnu[5]; 
	  for(int l=0;l<5;l++){ 
	    xnu[l]  = sigma*abs(c.lambda[l]); 
	    c.delwnu[l] = abs(c.lambda[l])*(1.-xnu[l])*c.delw[l]; 
	    //Computation of higher-order corrective terms
	    c.cf2[l]  = c.co[1]*abs(c.lambda[l])*(1.-xnu[l]); 
	    c.cf3[l]  = c.co[2]*c.cf2[l]*(1.+xnu[l])/3.; 
	    c.cf4[l]  = c.co[3]*c.cf3[l]*(xnu[l]-2.)/4.; 
	    c.cf5[l]  = c.co[4]*c.cf4[l]*(xnu[l]+2.)/5.; 
	    c.cf6[l]  = c.co[5]*c.cf5[l]*(xnu[l]-3.)/6.; 
	    c.cf7[l]  = c.co[6]*c.cf6[l]*(xnu[l]+3.)/7.; 
	    c.cf8[l]  = c.co[7]*c.cf7[l]*(xnu[l]-4.)/8.; 
	    c.cf9[l]  = c.co[8]*c.cf8[l]*(xnu[l]+4.)/9.; 
	    c.cf10[l] = c.co[9]*c.cf9[l]*(xnu[l]-5.)/10.; 
	    c.cf11[l] = c.co[10]*c.cf10[l]*(xnu[l]+5.)/11.; 
	  } 
										
										
	  for(int l=0;l<5;l++){ 
	    //Computation of centered corrections 
	    c.psic0[l] = (c.cf2[l]-2.*c.cf4[l]+6.*c.cf6[l]-20.*c.cf8[l]+70.*c.cf10[l])*c.delw[l]; 
	    c.psic1[l] = (c.cf4[l]-4.*c.cf6[l]+15.*c.cf8[l]-56.*c.cf10[l])*c.delw[l]; 
	    c.psic2[l] = (c.cf6[l]-6.*c.cf8[l]+28.*c.cf10[l])*c.delw[l]; 
	    c.psic3[l] = (c.cf8[l]-8.*c.cf10[l])*c.delw[l]; 
	    c.psic4[l] = (c.cf10[l])*c.delw[l]; 
	    //Computation of decentered corrections
	    c.psid0[l] = (126.*c.cf11[l]-35.*c.cf9[l]+10.*c.cf7[l]-3.*c.cf5[l]+c.cf3[l])*c.delw[l]; 
	    c.psid1[l] = (84.*c.cf11[l]-21.*c.cf9[l]+5.*c.cf7[l]-c.cf5[l])*c.delw[l]; 
	    c.psid2[l] = (36.*c.cf11[l]-7.*c.cf9[l]+c.cf7[l])*c.delw[l]; 
	    c.psid3[l] = (9.*c.cf11[l]-c.cf9[l])*c.delw[l]; 
	    c.psid4[l] = (c.cf11[l])*c.delw[l]; 
	  } 
										
	  //Computation of left eigenvectors 
	  c.vpr[0][0] = 1.; 
	  c.vpr[1][0] = ur-cr; 
	  c.vpr[2][0] = vr;
	  c.vpr[3][0] = wr;
	  c.vpr[4][0] = Hr-ur*cr; 
										
	  c.vpr[0][1] = 1.; 
	  c.vpr[1][1] = ur;
	  c.vpr[2][1] = vr;
	  c.vpr[3][1] = wr;
	  c.vpr[4][1] = ur2/2. + vr2/2. + wr2/2.;
										
	  c.vpr[0][2] = 0.; 
	  c.vpr[1][2] = 0.;
	  c.vpr[2][2] = cr;
	  c.vpr[3][2] = 0.;
	  c.vpr[4][2] = vr*cr;
										
	  c.vpr[0][3] = 0.; 
	  c.vpr[1][3] = 0.;
	  c.vpr[2][3] = 0.;
	  c.vpr[3][3] = cr;
	  c.vpr[4][3] = wr*cr;
										
	  c.vpr[0][4] = 1.; 
	  c.vpr[1][4] = ur+cr;
	  c.vpr[2][4] = vr;
	  c.vpr[3][4] = wr;
	  c.vpr[4][4] = Hr+ur*cr; 
										
	  //Computation of corrections in the eigenvectors basis
	  for(int l=0;l<5;l++){ 
	    c.psic0r[l] = 0.; 
	    c.psic1r[l] = 0.; 
	    c.psic2r[l] = 0.; 
	    c.psic3r[l] = 0.; 
	    c.psic4r[l] = 0.; 
	    c.psid0r[l] = 0.; 
	    c.psid1r[l] = 0.; 
	    c.psid2r[l] = 0.; 
	    c.psid3r[l] = 0.; 
	    c.psid4r[l] = 0.; 
	  } 
	  for(int m=0;m<5;m++){ 
	    for(int l=0;l<5;l++){ 
	      c.psic0r[m] += c.psic0[l]*c.vpr[m][l]; 
	      c.psic1r[m] += c.psic1[l]*c.vpr[m][l]; 
	      c.psic2r[m] += c.psic2[l]*c.vpr[m][l]; 
	      c.psic3r[m] += c.psic3[l]*c.vpr[m][l]; 
	      c.psic4r[m] += c.psic4[l]*c.vpr[m][l]; 
	      c.psid0r[m] += c.psid0[l]*c.vpr[m][l]; 
	      c.psid1r[m] += c.psid1[l]*c.vpr[m][l]; 
	      c.psid2r[m] += c.psid2[l]*c.vpr[m][l]; 
	      c.psid3r[m] += c.psid3[l]*c.vpr[m][l]; 
	      c.psid4r[m] += c.psid4[l]*c.vpr[m][l]; 
	    } 
	  } 
	}
      }
    }
  }    
    
  //Computation of the monotonicity indicators
    
  for(int l=0;l<5;l++){ 
    for(int i=1;i<Nx+2*marge-1;i++){ 
      for(int j=1;j<Ny+2*marge-1;j++){
	for(int k=1;k<Nz+2*marge-1;k++){
	  Cellule& c = grille[i][j][k]; 
	  Cellule& cg = grille[i-1][j][k]; 
	  if(!c.vide && !cg.vide){
	    c.am[l] = c.lambda[l]*c.delw[l]-cg.lambda[l]*cg.delw[l];
	  }
	}
      }
    } 
    //Computation of dj^m4 
    for(int i=1;i<Nx+2*marge-2;i++){ 
      for(int j=0;j<Ny+2*marge;j++){
	for(int k=0;k<Nz+2*marge;k++){   
	  Cellule& c = grille[i][j][k]; 
	  Cellule& cd = grille[i+1][j][k]; 
	  if(!c.vide && !cd.vide){ 
	    double z1 = 4.*c.am[l]-cd.am[l]; 
	    double z2 = 4.*cd.am[l]-c.am[l]; 
	    double z3 = c.am[l]; 
	    double z4 = cd.am[l]; 
	    c.am1[l] = (sign(z1)+sign(z2))/2.*abs((sign(z1)+sign(z3))/2.)*(sign(z1)
									   + sign(z4))/2.*min(abs(z1),min(abs(z2),min(abs(z3),abs(z4))));
	  }
	}
      }
    } 
  } 
    
  //Computation of r+ and r- 
  for(int l=0;l<5;l++){ 
    for(int i=marge;i<Nx+2*marge-4;i++){ 
      for(int j=marge;j<Ny+2*marge-4;j++){ 
	for(int k=marge;k<Nz+2*marge-4;k++){
	  Cellule& c = grille[i][j][k]; 
	  Cellule& cd = grille[i+1][j][k]; 
	  Cellule& cg = grille[i-1][j][k]; 
	  if(!c.vide && !cd.vide){ 
	    c.rp[l] = sign(c.delw[l])*sign(cg.delw[l])*(abs(cg.delw[l])+eps)/(abs(c.delw[l])+eps); 
	    c.rm[l] = sign(c.delw[l])*sign(cd.delw[l])*(abs(cd.delw[l])+eps)/(abs(c.delw[l])+eps); 
	    //Higher-order corrections 
	    Cellule&  cg2 = grille[i-2][j][k]; 
	    Cellule& cg3 = grille[i-3][j][k]; 
	    Cellule& cg4 = grille[i-4][j][k]; 
	    Cellule& cg5 = grille[i-5][j][k]; 
	    Cellule& cd2 = grille[i+2][j][k]; 
	    Cellule& cd3 = grille[i+3][j][k]; 
	    Cellule& cd4 = grille[i+4][j][k]; 
	    c.psid[l] = -c.psid0[l]+cg.psid0[l]+cd.psid1[l]-cg2.psid1[l]-cd2.psid2[l]+cg3.psid2[l]
	      + cd3.psid3[l]-cg4.psid3[l]-cd4.psid4[l]+cg5.psid4[l];
	  }
	}
      }
    } 
  } 
    
  //Flux computation 
  for(int i=marge-1;i<Nx+marge;i++){ 
    for(int j=marge-1;j<Ny+marge;j++){ 
      for(int k=marge-1;k<Nz+marge;k++){ 
	Cellule&  c = grille[i][j][k]; 
	//Neighbouring cells 
	Cellule& cg = grille[i-1][j][k]; 
	Cellule& cg2 = grille[i-2][j][k]; 
	Cellule& cg3 = grille[i-3][j][k]; 
	Cellule& cg4 = grille[i-4][j][k]; 
	Cellule& cd = grille[i+1][j][k]; 
	Cellule& cd2 = grille[i+2][j][k]; 
	Cellule& cd3 = grille[i+3][j][k]; 
	Cellule& cd4 = grille[i+4][j][k]; 
                
	//TVD flux 
	double tvd[5]; 
	double psict[5]; 
	if(!c.vide && !cd.vide){   
	  //Initialization 
	  for(int l=0; l<5; l++){ 
	    tvd[l] = 0.;
	    //Centered part
	    psict[l] = c.psic0r[l] + cg.psic1r[l] + cd.psic1r[l] + cg2.psic2r[l] + cd2.psic2r[l] 
	      + cg3.psic3r[l] + cd3.psic3r[l] + cg4.psic4r[l] + cd4.psic4r[l]; 
	  } 
                
	  //Limiter 
	  double psic; 
	  for(int l=0; l<5; l++){ 
	    psic = c.psic0[l] + cg.psic1[l] + cd.psic1[l] + cg2.psic2[l] + cd2.psic2[l] 
	      + cg3.psic3[l] + cd3.psic3[l] + cg4.psic4[l] + cd4.psic4[l]; 
	    //Decentered part
	    double r; 
	    double xnum; 
	    double xnume; 
	    int is; 
	    double psi; 
	    if(c.lambda[l]>0.){ 
	      r = c.rp[l]; 
	      xnum = sigma*abs(cg.lambda[l]); 
	      xnume = max(xnum,eps); 
	      is = 1; 
	      psi = psic+c.psid[l]; 
	    } else { 
	      r = c.rm[l]; 
	      xnum = sigma*abs(cd.lambda[l]); 
	      xnume = max(xnum,eps); 
	      is = -1; 
	      psi = psic-cd.psid[l]; 
	    } 
                    
	    double xnu = sigma*abs(c.lambda[l]); 
	    xnu = max(xnu,eps); 
                    
	    //TVD limiter psitvd 
	    psi = (double) sign(c.delwnu[l])*psi/(abs(c.delwnu[l]+eps)); 
	    double psimax1 = 2.*r*(1.-xnume)/(xnu*(1.-xnu)); 
	    double psimax2 = 2./(1.-xnu); 
	    double psitvd = max(0.,min(psi,min(psimax1,psimax2))); 
                    
	    //Monotonicity criterion
	    if((c.delwnu[l] != 0.) && (abs(psi-psitvd)>eps)){ 
	      double dfo = psi*c.delwnu[l]/2.; 
	      double dabsf = psimax2*c.delwnu[l]/2.; 
	      double dful = psimax1*c.delwnu[l]/2.; 
	      double dfmd = dabsf/2.-c.am1[l]/2.; 
	      Cellule& camont = grille[i-is][j][k];   //Upwind cell
	      double dflc = dful/2.+((1.-xnume)/xnu)*camont.am1[l]/2.; 
	      double dfmin = max(min(0.,min(dabsf,dfmd)),min(0.,min(dful,dflc))); 
	      double dfmax = min(max(0.,max(dabsf,dfmd)),max(0.,max(dful,dflc))); 
	      if((dfmin-dfo)*(dfmax-dfo)>0.){ 
		psi = psitvd; 
	      } 
	    } 
                    
	    //Uncomment to use only the TVD part and not the MP part
	    //psi = psitvd; 
                    
	    //Uncomment to disable both TVD and MP
	    //psi = 0.; 
                    
	    double ctvd = psi*c.delwnu[l]/2.-abs(c.lambda[l])*c.delw[l]/2.; 
	    for(int m=0;m<5;m++){ 
	      tvd[m] += ctvd*c.vpr[m][l];
	    } 
	  } 
                
	  // Final computation of the flux 
	  for(int l=0;l<5;l++){ 
	    c.fluxi[l] += tvd[l]; 
	  }
                
	}
      }
    }
  } 
    
  //Entropy correction
    
  corentx(sigma);
  
//Boundary conditions

  //Reflecting boundary conditions
  if(BC_x_in ==  1 || BC_x_out ==  1){
    for(int j=0;j<Ny+2*marge;j++){
      for(int k=0;k<Nz+2*marge;k++){
                
	if(BC_x_in ==  1){
	  Cellule& c = grille[marge-1][j][k];
	  Cellule& cp = grille[marge][j][k];
	  double p0 = (gam-1.)*(cp.rhoE0-1./2.*(cp.impx0*cp.impx0 + cp.impy0*cp.impy0 
						+ cp.impz0*cp.impz0)/cp.rho0);
	  c.fluxi[0] = 0.;
	  c.fluxi[1] = p0;
	  c.fluxi[2] = 0.;
	  c.fluxi[3] = 0.;
	  c.fluxi[4] = 0.;
                    
	}
                
	if(BC_x_out ==  1){
	  Cellule& c2 = grille[Nx+marge-1][j][k];	  
	  double p02 = (gam-1.)*(c2.rhoE0-1./2.*(c2.impx0*c2.impx0 + c2.impy0*c2.impy0 
						 + c2.impz0*c2.impz0)/c2.rho0);;
	  c2.fluxi[0] = 0.;
	  c2.fluxi[1] = p02;
	  c2.fluxi[2] = 0.;
	  c2.fluxi[3] = 0.;
	  c2.fluxi[4] = 0.;
                    
	}
                
      }
    }
  }

  //Periodic boundary conditions
  if(BC_x_in ==  2 || BC_x_out ==  2){
    for(int i=0;i<Nx+2*marge;i++){
      for(int j=0;j<Ny+2*marge;j++){
	for(int k=0;k<Nz+2*marge;k++){
	  if(i==Nx+marge-1){
	    grille[i][j][k].fluxi[0] = grille[marge-1][j][k].fluxi[0];
	    grille[i][j][k].fluxi[1] = grille[marge-1][j][k].fluxi[1];
	    grille[i][j][k].fluxi[2] = grille[marge-1][j][k].fluxi[2];
	    grille[i][j][k].fluxi[3] = grille[marge-1][j][k].fluxi[3];
	    grille[i][j][k].fluxi[4] = grille[marge-1][j][k].fluxi[4];
	  }
	}
      }
    }
  }

  //Outflow boundary conditions (Poinsot-Lele)
  if(BC_x_out==3){
    int i=Nx+marge-1;
    for(int j=0;j<Ny+2*marge;j++){
      for(int k=0;k<Nz+2*marge;k++){
	Cellule& c = grille[i][j][k];
	if(c.y>0.2 || c.y<0.1 || c.z<0.09 || c.x>0.1){
	  if(c.p<eps || c.rho<eps){
	    cout << "p or rho negative : p " << c.p << " rho " << c.rho;
	    getchar();
	  }
	  double cr = sqrt(gam*c.p/c.rho);
	  Cellule& cg = grille[i-1][j][k];
	  Cellule& cg2 = grille[i-2][j][k];
	  double drho = (c.rho-cg2.rho)/dx/2.;
	  double du = (c.u-cg2.u)/dx/2.;
	  double dv = (c.v-cg2.v)/dx/2.;
	  double dw = (c.w-cg2.w)/dx/2.;
	  double dp = (c.p-cg2.p)/dx/2.;
	  double L = min(min(domainex,domainey),domainez);
	  double k = 0.278;
	  double alpha = k*(cr*cr-c.u*c.u)/L/cr;
	  double pinf = P(c.x,c.y,c.z,dx,dy,dz);
	  double L0 = c.u*(cr*cr*drho-dp);
	  double L1 = (c.u+cr)*(c.rho*cr*du+dp);
	  double L2 = c.rho*c.u*cr*dv;
	  double L3 = c.rho*c.u*cr*dw;
	  double L4 = alpha*(c.p-pinf);
	  if(c.u>0. && c.u-cr<0.){
	    //Subsonic outflow
	    L0 = c.u*(cr*cr*drho-dp);
	    L1 = (c.u+cr)*(c.rho*cr*du+dp);
	    L2 = c.rho*c.u*cr*dv;
	    L3 = c.rho*c.u*cr*dw;
	    L4 = alpha*(c.p-pinf);
	  }
	  else if(c.u>0. && c.u-cr>0.){
	    //Supersonic outflow
	    L0 = c.u*(cr*cr*drho-dp);
	    L1 = (c.u+cr)*(c.rho*cr*du+dp);
	    L2 = c.rho*c.u*cr*dv;
	    L3 = c.rho*c.u*cr*dw;
	    L4 = (c.u-cr)*(-c.rho*cr*du+dp);
	  }
	  else if(c.u<0. && c.u+cr>0.){
	    //Subsonic inflow
	    L0 = 0.;
	    L1 = alpha*(c.p-pinf);
	    L2 = 0.;
	    L3 = 0.;
	    L4 = L1;
	  }
	  else if(c.u<0. && c.u+cr<0.){
	    //Supersonic inflow
	    L0 = 0.;
	    L1 = alpha*(c.p-pinf);
	    L2 = 0.;
	    L3 = 0.;
	    L4 = L1;
	  }
	  double d0 = (L0+L1/2.+L4/2.)/cr/cr;
	  double d1 = (L1-L4)/2./c.rho/cr;
	  double d2 = L2/c.rho/cr;
	  double d3 = L3/c.rho/cr;
	  double d4 = (L1+L4)/2.;
	  d0 *= dx;
	  d1 *= dx;
	  d2 *= dx;
	  d3 *= dx;
	  d4 *= dx;
	  double kappa = 0.05;
	  c.fluxi[0] = cg.fluxi[0]+d0;
	  c.fluxi[1] = cg.fluxi[1]+(c.u*d0+c.rho*d1);
	  c.fluxi[2] = cg.fluxi[2]+(c.v*d0+c.rho*d2);
	  c.fluxi[3] = cg.fluxi[3]+(c.w*d0+c.rho*d3);
	  c.fluxi[4] = cg.fluxi[4]+((c.u*c.u+c.v*c.v+c.w*c.w)/2.*d0+c.rho*c.u*d1+c.rho*c.v*d2+c.rho*c.w*d3+d4/(gam-1.));
	}
      }
    }
  }
  if(BC_x_in==3){
    int i=marge-1;
    for(int j=0;j<Ny+2*marge;j++){
      for(int k=0;k<Nz+2*marge;k++){
	Cellule& c = grille[i][j][k];
	if(c.y>0.2 || c.y<0.1 || c.z<0.09 || c.x>0.1){
	  if(c.p<eps || c.rho<eps){
	    cout << "p or rho negative : p " << c.p << " rho " << c.rho;
	    getchar();
	  }
	  Cellule& cd = grille[i+1][j][k];
	  Cellule& cd2 = grille[i+2][j][k];
	  double cr = sqrt(gam*cd.p/cd.rho);
	  double drho = (cd2.rho-c.rho)/dx/2.;
	  double du = (cd2.u-c.u)/dx/2.;
	  double dv = (cd2.v-c.v)/dx/2.;
	  double dw = (cd2.w-c.w)/dx/2.;
	  double dp = (cd2.p-c.p)/dx/2.;
	  double L = min(min(domainex,domainey),domainez);
	  double k = 0.278;
	  double alpha = k*(cr*cr-c.u*c.u)/L/cr;
	  double pinf = P(cd.x,cd.y,cd.z,dx,dy,dz);
	  double L0 = cd.u*(cr*cr*drho-dp);
	  double L4 = (cd.u-cr)*(-cd.rho*cr*du+dp);
	  double L2 = cd.rho*cd.u*cr*dv;
	  double L3 = cd.rho*cd.u*cr*dw;
	  double L1 = alpha*(cd.p-pinf);
	  if(cd.u<0. && c.u-cr>0.){
	    //Subsonic outflow
	    L0 = cd.u*(cr*cr*drho-dp);
	    L4 = (cd.u-cr)*(-cd.rho*cr*du+dp);
	    L2 = cd.rho*cd.u*cr*dv;
	    L3 = cd.rho*cd.u*cr*dw;
	    L1 = alpha*(cd.p-pinf);
	  }
	  else if(cd.u<0. && cd.u-cr<0.){
	    //Supersonic outflow
	    L0 = cd.u*(cr*cr*drho-dp);
	    L4 = (cd.u-cr)*(-cd.rho*cr*du+dp);
	    L2 = cd.rho*cd.u*cr*dv;
	    L3 = cd.rho*cd.u*cr*dw;
	    L4 = (cd.u+cr)*(cd.rho*cr*du+dp);
	  }
	  else if(cd.u>0. && cd.u-cr<0.){
	    //Subsonic inflow
	    L0 = 0.;
	    L4 = alpha*(cd.p-pinf);
	    L2 = 0.;
	    L3 = 0.;
	    L1 = L4;
	  }
	  else if(cd.u>0. && cd.u-cr>0.){
	    //Supersonic inflow
	    L0 = 0.;
	    L4 = alpha*(c.p-pinf);
	    L2 = 0.;
	    L3 = 0.;
	    L1 = L4;
	  }
	  double d0 = (L0+L1/2.+L4/2.)/cr/cr;
	  double d1 = (L1-L4)/2./cd.rho/cr;
	  double d2 = L2/cd.rho/cr;
	  double d3 = L3/cd.rho/cr;
	  double d4 = (L1+L4)/2.;
	  d0 *= dx;
	  d1 *= dx;
	  d2 *= dx;
	  d3 *= dx;
	  d4 *= dx;
	  double kappa = 0.05;
	  c.fluxi[0] = cd.fluxi[0]-d0;
	  c.fluxi[1] = cd.fluxi[1]-(cd.u*d0+cd.rho*d1);
	  c.fluxi[2] = cd.fluxi[2]-(cd.v*d0+cd.rho*d2);
	  c.fluxi[3] = cd.fluxi[3]-(cd.w*d0+cd.rho*d3);
	  c.fluxi[4] = cd.fluxi[4]-((cd.u*cd.u+cd.v*cd.v+cd.w*cd.w)/2.*d0+cd.rho*cd.u*d1+cd.rho*cd.v*d2+cd.rho*cd.w*d3+d4/(gam-1.));
	}
      }
    }
  }
}


/*!\brief Computation of the numerical y-flux.
   \param sigma = \a dt/dy: time-step/fluid spatial discretization step in the y direction
   \param t current simulation time
   \return void
*/
void Grille::fnumy(const double sigma, double t){ 
  //Initialization to the centered flux
  for(int i=0; i<Nx+2*marge-1; i++){
    for(int j=0; j<Ny+2*marge-1; j++){ 
      for(int k=0; k<Nz+2*marge-1; k++){
	Cellule& c = grille[i][j][k];    
	Cellule& cj = grille[i][j+1][k];
	//Order indicators
	if(!c.vide && !cj.vide){
	  for(int l=0;l< c.ordre;l++){ 
	    c.co[l]=1.; 
	  } 
	  for(int l=c.ordre;l<ordremax;l++){ 
	    c.co[l]=0.;
	  } 
									
	  //Centered part of the flux  
	  c.fluxj[0] = (c.impy+cj.impy)/2.; 
	  c.fluxj[1] = (c.rho*c.u*c.v+cj.rho*cj.u*cj.v)/2.;
	  c.fluxj[2] = (c.rho*c.v*c.v +c.p +cj.rho*cj.v*cj.v + cj.p)/2.;
	  c.fluxj[3] = (c.rho*c.v*c.w+cj.rho*cj.v*cj.w)/2.;
	  c.fluxj[4] = ((c.rhoE+c.p)*c.v+(cj.rhoE+cj.p)*cj.v)/2.; 
	}
	else if(std::abs(c.alpha-1.)>eps){
									
	  //Lax-Friedrichs flux near void
	  c.fluxj[0] = (c.impy+cj.impy)/2. + (c.rho -cj.rho)/2./sigma; 
	  c.fluxj[1] = (c.rho*c.u*c.v+cj.rho*cj.u*cj.v)/2. + (c.impx -cj.impx)/2./sigma; 
	  c.fluxj[2] = (c.rho*c.v*c.v +c.p +cj.rho*cj.v*cj.v+cj.p)/2. + (c.impy -cj.impy)/2./sigma; 
	  c.fluxj[3] = (c.rho*c.v*c.w+cj.rho*cj.v*cj.w)/2. + (c.impz -cj.impz)/2./sigma; 
	  c.fluxj[4] = ((c.rhoE+c.p)*c.u+(cj.rhoE+cj.p)*cj.u)/2. + (c.rhoE -cj.rhoE)/2./sigma; 
	}
	else{
	  c.fluxj[0] = 0.;
	  c.fluxj[1] = 0.;
	  c.fluxj[2] = 0.;
	  c.fluxj[3] = 0.;
	  c.fluxj[4] = 0.;
	}
      }
    }
  } 
    
  //Preliminary computation of limiter variables 
  for(int i=0; i<Nx+2*marge-1; i++){
    for(int j=0; j<Ny+2*marge-1; j++){
      for(int k=0; k<Nz+2*marge-1; k++){
	Cellule& cj = grille[i][j+1][k]; 
	Cellule& c = grille[i][j][k]; 
	if(!c.vide && !cj.vide){ 
	  //Computation of the Roe variables 
	  double roe = sqrt(cj.rho/c.rho); 
	  double rhor = roe*c.rho; 
	  double ur = (roe*cj.u+c.u)/(1.+roe); 
	  double vr = (roe*cj.v+c.v)/(1.+roe); 
	  double wr = (roe*cj.w+c.w)/(1.+roe); 
	  double Hr = (roe*(cj.rho*cj.u*cj.u/2.+ cj.rho*cj.v*cj.v/2.+ cj.rho*cj.w*cj.w/2. + cj.p*gam/(gam-1.))/cj.rho + (c.rho*c.u*c.u/2.+c.rho*c.v*c.v/2.+ c.rho*c.w*c.w/2. + c.p*gam/(gam-1.))/c.rho)/(1.+roe);
	  double ur2 = ur*ur; 
	  double vr2 = vr*vr;
	  double wr2 = wr*wr;
	  double cr2 = (gam-1.)*(Hr-ur2/2.-vr2/2.-wr2/2.); 
									
	  //Test on the speed of sound
	  if(cr2<=0. && abs(c.alpha-1.)>eps){
	    cout << "computation of the y-flux" << endl;
	    cout << "i=" << i << " j=" << j <<" k = "<<k<< " negative speed of sound: c2=" << cr2 << endl;
	    cout << "x=" << c.x << " y=" << c.y << " z=" << c.z << " alpha=" << c.alpha << endl;
	    cout << "t=" << t << endl; 
	    cout << "c.p=" << c.p << endl; 
	    cout << "c.rho=" << c.rho << endl; 
	    cout << "c.u=" << c.u << endl; 
	    cout << "c.v=" << c.v << endl; 
	    cout << "c.w=" << c.w << endl; 
	    cout << "cj.p=" << cj.p << endl; 
	    cout << "cj.rho=" << cj.rho << endl; 
	    cout << "cj.u=" << cj.u << endl;
	    cout << "cj.v=" << cj.v << endl;
	    cout << "cj.w=" << cj.w << endl;
	    cout << "ur=" << ur << endl; 
	    cout << "ur2=" << ur2 << endl; 
	    cout << "vr=" << vr << endl; 
	    cout << "vr2=" << vr2 << endl; 
	    cout << "wr2=" << wr2 << endl; 
	    cout << "Hr=" << Hr << endl;
	    getchar(); 
	  } 
									
	  //Speed of sound 
	  double cr = sqrt(cr2); 
									
	  //Eigenvalues 
	  c.lambda[0] = vr-cr; 
	  c.lambda[1] = vr;
	  c.lambda[2] = vr;
	  c.lambda[3] = vr;
	  c.lambda[4] = vr+cr; 
									
	  //Computation of differences between Wd and Wg 
	  double drho = cj.rho - c.rho; 
	  double du = cj.u - c.u;
	  double dv = cj.v - c.v;
	  double dw = cj.w - c.w;
	  double dp = cj.p - c.p; 
									
	  //Computation of the deltaV (differences between Wd and Wg in the eigenvectors basis) 
	  double ros2c = rhor/cr/2.; 
	  c.delw[0] = dp/cr2/2. - ros2c*dv; 
	  c.delw[1] = drho - dp/cr2;
	  c.delw[2] = 2.*ros2c*du;
	  c.delw[3] = 2.*ros2c*dw;   
	  c.delw[4] = dp/cr2/2. + ros2c*dv; 
									
	  //Computation of the complete correction in the eigenvectors basis
	  double xnu[5]; 
	  for(int l=0;l<5;l++){ 
	    xnu[l]  = sigma*abs(c.lambda[l]); 
	    c.delwnu[l] = abs(c.lambda[l])*(1.-xnu[l])*c.delw[l]; 
	    //Higher-order correction coefficients
	    c.cf2[l]  = c.co[1]*abs(c.lambda[l])*(1.-xnu[l]); 
	    c.cf3[l]  = c.co[2]*c.cf2[l]*(1.+xnu[l])/3.; 
	    c.cf4[l]  = c.co[3]*c.cf3[l]*(xnu[l]-2.)/4.; 
	    c.cf5[l]  = c.co[4]*c.cf4[l]*(xnu[l]+2.)/5.; 
	    c.cf6[l]  = c.co[5]*c.cf5[l]*(xnu[l]-3.)/6.; 
	    c.cf7[l]  = c.co[6]*c.cf6[l]*(xnu[l]+3.)/7.; 
	    c.cf8[l]  = c.co[7]*c.cf7[l]*(xnu[l]-4.)/8.; 
	    c.cf9[l]  = c.co[8]*c.cf8[l]*(xnu[l]+4.)/9.; 
	    c.cf10[l] = c.co[9]*c.cf9[l]*(xnu[l]-5.)/10.; 
	    c.cf11[l] = c.co[10]*c.cf10[l]*(xnu[l]+5.)/11.; 
	  } 
									
									
	  for(int l=0;l<5;l++){ 
	    //Centered corrections
	    c.psic0[l] = (c.cf2[l]-2.*c.cf4[l]+6.*c.cf6[l]-20.*c.cf8[l]+70.*c.cf10[l])*c.delw[l]; 
	    c.psic1[l] = (c.cf4[l]-4.*c.cf6[l]+15.*c.cf8[l]-56.*c.cf10[l])*c.delw[l]; 
	    c.psic2[l] = (c.cf6[l]-6.*c.cf8[l]+28.*c.cf10[l])*c.delw[l]; 
	    c.psic3[l] = (c.cf8[l]-8.*c.cf10[l])*c.delw[l]; 
	    c.psic4[l] = (c.cf10[l])*c.delw[l]; 
	    //Decentered corrections 
	    c.psid0[l] = (126.*c.cf11[l]-35.*c.cf9[l]+10.*c.cf7[l]-3.*c.cf5[l]+c.cf3[l])*c.delw[l]; 
	    c.psid1[l] = (84.*c.cf11[l]-21.*c.cf9[l]+5.*c.cf7[l]-c.cf5[l])*c.delw[l]; 
	    c.psid2[l] = (36.*c.cf11[l]-7.*c.cf9[l]+c.cf7[l])*c.delw[l]; 
	    c.psid3[l] = (9.*c.cf11[l]-c.cf9[l])*c.delw[l]; 
	    c.psid4[l] = (c.cf11[l])*c.delw[l]; 
	  } 
									
	  //Left eigenvalues 
	  c.vpr[0][0] = 1.; 
	  c.vpr[1][0] = ur; 
	  c.vpr[2][0] = vr-cr;
	  c.vpr[3][0] = wr;
	  c.vpr[4][0] = Hr-vr*cr; 
									
	  c.vpr[0][1] = 1.; 
	  c.vpr[1][1] = ur;
	  c.vpr[2][1] = vr;
	  c.vpr[3][1] = wr;
	  c.vpr[4][1] = ur2/2. + vr2/2. + wr2/2.;
									
	  c.vpr[0][2] = 0.; 
	  c.vpr[1][2] = 0.;
	  c.vpr[2][2] = cr;
	  c.vpr[3][2] = 0.;
	  c.vpr[4][2] = ur*cr;
									
	  c.vpr[0][3] = 0.; 
	  c.vpr[1][3] = 0.;
	  c.vpr[2][3] = 0.;
	  c.vpr[3][3] = cr;
	  c.vpr[4][3] = wr*cr;
									
	  c.vpr[0][4] = 1.; 
	  c.vpr[1][4] = ur;
	  c.vpr[2][4] = vr+cr;
	  c.vpr[3][4] = wr;
	  c.vpr[4][4] = Hr+vr*cr; 
									
	  //Corrections in the eigenvectors basis
	  for(int l=0;l<5;l++){ 
	    c.psic0r[l] = 0.; 
	    c.psic1r[l] = 0.; 
	    c.psic2r[l] = 0.; 
	    c.psic3r[l] = 0.; 
	    c.psic4r[l] = 0.; 
	    c.psid0r[l] = 0.; 
	    c.psid1r[l] = 0.; 
	    c.psid2r[l] = 0.; 
	    c.psid3r[l] = 0.; 
	    c.psid4r[l] = 0.; 
	  } 
	  for(int m=0;m<5;m++){ 
	    for(int l=0;l<5;l++){ 
	      c.psic0r[m] += c.psic0[l]*c.vpr[m][l]; 
	      c.psic1r[m] += c.psic1[l]*c.vpr[m][l]; 
	      c.psic2r[m] += c.psic2[l]*c.vpr[m][l]; 
	      c.psic3r[m] += c.psic3[l]*c.vpr[m][l]; 
	      c.psic4r[m] += c.psic4[l]*c.vpr[m][l]; 
	      c.psid0r[m] += c.psid0[l]*c.vpr[m][l]; 
	      c.psid1r[m] += c.psid1[l]*c.vpr[m][l]; 
	      c.psid2r[m] += c.psid2[l]*c.vpr[m][l]; 
	      c.psid3r[m] += c.psid3[l]*c.vpr[m][l]; 
	      c.psid4r[m] += c.psid4[l]*c.vpr[m][l]; 
	    } 
	  } 
	}
      }
    }
  }    
    
  //Computation of the monotonicity indicators
  for(int l=0;l<5;l++){ 
    for(int i=1;i<Nx+2*marge-1;i++){
      for(int j=1;j<Ny+2*marge-1;j++){ 
	for(int k=1;k<Nz+2*marge-1;k++){
	  Cellule& c = grille[i][j][k]; 
	  Cellule& cg = grille[i][j-1][k]; 
	  if(!c.vide && !cg.vide){
	    c.am[l] = c.lambda[l]*c.delw[l]-cg.lambda[l]*cg.delw[l]; 
	  }
	}
      }
    } 
    //Computation of dj^m4 
    for(int i=0;i<Nx+2*marge;i++){
      for(int j=1;j<Ny+2*marge-2;j++){ 
	for(int k=0;k<Nz+2*marge;k++){   
	  Cellule& c = grille[i][j][k]; 
	  Cellule&  cd = grille[i][j+1][k]; 
	  if(!c.vide && !cd.vide){
	    double z1 = 4.*c.am[l]-cd.am[l]; 
	    double z2 = 4.*cd.am[l]-c.am[l]; 
	    double z3 = c.am[l]; 
	    double z4 = cd.am[l]; 
	    c.am1[l] = (sign(z1)+sign(z2))/2.*abs((sign(z1)+sign(z3))/2.)*(sign(z1)
									   + sign(z4))/2.*min(abs(z1),min(abs(z2),min(abs(z3),abs(z4))));
	  }
	}
      }
    } 
  } 
  //Computation of r+ and r- 
  for(int l=0;l<5;l++){ 
    for(int i=marge;i<Nx+2*marge-4;i++){ 
      for(int j=marge;j<Ny+2*marge-4;j++){ 
	for(int k=marge;k<Nz+2*marge-4;k++){
	  Cellule& c = grille[i][j][k]; 
	  Cellule& cd = grille[i][j+1][k]; 
	  Cellule& cg = grille[i][j-1][k]; 
	  if(!c.vide && !cd.vide){ 
	    c.rp[l] = sign(c.delw[l])*sign(cg.delw[l])*(abs(cg.delw[l])+eps)/(abs(c.delw[l])+eps); 
	    c.rm[l] = sign(c.delw[l])*sign(cd.delw[l])*(abs(cd.delw[l])+eps)/(abs(c.delw[l])+eps); 
	    //Higher-order corrections 
	    Cellule&  cg2 = grille[i][j-2][k]; 
	    Cellule& cg3 = grille[i][j-3][k]; 
	    Cellule& cg4 = grille[i][j-4][k]; 
	    Cellule& cg5 = grille[i][j-5][k]; 
	    Cellule& cd2 = grille[i][j+2][k]; 
	    Cellule& cd3 = grille[i][j+3][k]; 
	    Cellule& cd4 = grille[i][j+4][k]; 
	    c.psid[l] = - c.psid0[l]+cg.psid0[l]+cd.psid1[l]-cg2.psid1[l]-cd2.psid2[l]
	      + cg3.psid2[l]+cd3.psid3[l]-cg4.psid3[l]-cd4.psid4[l]+cg5.psid4[l];
	  }
	}
      }
    } 
  } 
		
  //Flux computation 
  for(int i=marge-1;i<Nx+marge;i++){ 
    for(int j=marge-1;j<Ny+marge;j++){ 
      for(int k=marge-1;k<Nz+marge;k++){ 
	Cellule&   c = grille[i][j][k]; 
	//Neighbouring cells 
	Cellule&  cg = grille[i][j-1][k]; 
	Cellule& cg2 = grille[i][j-2][k]; 
	Cellule& cg3 = grille[i][j-3][k]; 
	Cellule& cg4 = grille[i][j-4][k]; 
	Cellule& cd = grille[i][j+1][k]; 
	Cellule& cd2 = grille[i][j+2][k]; 
	Cellule& cd3 = grille[i][j+3][k]; 
	Cellule& cd4 = grille[i][j+4][k]; 
                
	//TVD flux 
	double tvd[5]; 
	double psict[5]; 
	if(!c.vide && !cd.vide){  
	  //Initialization 
	  for(int l=0;l<5;l++){ 
	    tvd[l] = 0.; 
	    //Centered part
	    psict[l] = c.psic0r[l] + cg.psic1r[l] + cd.psic1r[l] + cg2.psic2r[l] + cd2.psic2r[l]
	      + cg3.psic3r[l] + cd3.psic3r[l] + cg4.psic4r[l] + cd4.psic4r[l]; 
	  } 
									
	  //Limiter 
	  double psic; 
	  for(int l=0;l<5;l++){ 
	    psic = c.psic0[l] + cg.psic1[l] + cd.psic1[l] + cg2.psic2[l] + cd2.psic2[l]  + cg3.psic3[l] + cd3.psic3[l] + cg4.psic4[l] + cd4.psic4[l]; 
											
	    //Decentered part
	    double r; 
	    double xnum; 
	    double xnume; 
	    int is; 
	    double psi; 
	    if(c.lambda[l]>0.){ 
	      r = c.rp[l]; 
	      xnum = sigma*abs(cg.lambda[l]); 
	      xnume = max(xnum,eps); 
	      is = 1; 
	      psi = psic+c.psid[l]; 
	    } else { 
	      r = c.rm[l]; 
	      xnum = sigma*abs(cd.lambda[l]); 
	      xnume = max(xnum,eps); 
	      is = -1; 
	      psi = psic-cd.psid[l]; 
	    } 
	    
	    double xnu = sigma*abs(c.lambda[l]); 
	    xnu = max(xnu,eps); 
	    
	    //Computation of TVD limiter psitvd 
	    psi = (double) sign(c.delwnu[l])*psi/(abs(c.delwnu[l]+eps)); 
	    double psimax1 = 2.*r*(1.-xnume)/(xnu*(1.-xnu)); 
	    double psimax2 = 2./(1.-xnu); 
	    double psitvd = max(0.,min(psi,min(psimax1,psimax2))); 
											
	    //Monotonicity criterion
	    if((c.delwnu[l] != 0.) && (abs(psi-psitvd)>eps)){ 
	      double dfo = psi*c.delwnu[l]/2.; 
	      double dabsf = psimax2*c.delwnu[l]/2.; 
	      double dful = psimax1*c.delwnu[l]/2.; 
	      double dfmd = dabsf/2.-c.am1[l]/2.; 
	      Cellule& camont = grille[i][j-is][k];   //Upwind cell
	      double dflc = dful/2.+((1.-xnume)/xnu)*camont.am1[l]/2.; 
	      double dfmin = max(min(0.,min(dabsf,dfmd)),min(0.,min(dful,dflc))); 
	      double dfmax = min(max(0.,max(dabsf,dfmd)),max(0.,max(dful,dflc))); 
	      if((dfmin-dfo)*(dfmax-dfo)>0.){ 
		psi = psitvd; 
	      } 
	    } 
	    
	    //Uncomment to use only TVD without MP 
	    //psi = psitvd; 
	    
	    //Uncomment to use the scheme without TVD nor MP
	    //psi = 0.; 
	    
	    double ctvd = psi*c.delwnu[l]/2.-abs(c.lambda[l])*c.delw[l]/2.; 
	    for(int m=0;m<5;m++){ 
	      tvd[m] += ctvd*c.vpr[m][l];
	    } 
	  }
	  
	  //Final computation of the flux
	  for(int l=0;l<5;l++){ 
	    c.fluxj[l] += tvd[l]; 
	  }
	}
      }
    }
  }
  
  //Entropy correction
  corenty(sigma);
    

  //Boundary conditions
  //Reflecting boundary conditions
  if(BC_y_in ==  1 || BC_y_out ==  1){
    for(int i=0;i<Nx+2*marge;i++){
      for(int k=0;k<Nz+2*marge;k++){
                
	if(BC_y_in ==  1){
	  Cellule& c = grille[i][marge-1][k];
	  Cellule& cp = grille[i][marge][k];
	  double p0 = (gam-1.)*(cp.rhoE0 - 1./2.*(cp.impx0*cp.impx0 + cp.impy0*cp.impy0 + cp.impz0*cp.impz0)/cp.rho0);
	  c.fluxj[0] = 0.;
	  c.fluxj[1] = 0.;
	  c.fluxj[2] = p0;
	  c.fluxj[3] = 0.;
	  c.fluxj[4] = 0.;
        }
	if(BC_y_out ==  1){
	  Cellule& c2 = grille[i][Ny+marge-1][k];
	  double p02 = (gam-1.)*(c2.rhoE0-1./2.*(c2.impx0*c2.impx0 + c2.impy0*c2.impy0 + c2.impz0*c2.impz0)/c2.rho0);
	  c2.fluxj[0] = 0.;
	  c2.fluxj[1] = 0.;
	  c2.fluxj[2] = p02;
	  c2.fluxj[3] = 0.;
	  c2.fluxj[4] = 0.;
       	}
      }
    }
  }
  //Periodic boundary conditions
  if(BC_y_in ==  2 || BC_y_out ==  2){
    for(int i=0;i<Nx+2*marge;i++){
      for(int j=0;j<Ny+2*marge;j++){
	for(int k=0;k<Nz+2*marge;k++){
	  if(j==Ny+marge-1){
	    grille[i][j][k].fluxj[0] = grille[i][marge-1][k].fluxj[0];
	    grille[i][j][k].fluxj[1] = grille[i][marge-1][k].fluxj[1];
	    grille[i][j][k].fluxj[2] = grille[i][marge-1][k].fluxj[2];
	    grille[i][j][k].fluxj[3] = grille[i][marge-1][k].fluxj[3];
	    grille[i][j][k].fluxj[4] = grille[i][marge-1][k].fluxj[4];
	  }
	}
      }
    }
  }
  //Outflow conditions (Poinsot-Lele)
  if(BC_y_out==3){
    int j=Ny+marge-1;
    for(int k=0;k<Nz+2*marge;k++){
      for(int i=0;i<Nx+2*marge;i++){
	Cellule& c = grille[i][j][k];
	if(c.p<eps || c.rho<eps){
	  cout << "p or rho negative: p " << c.p << " rho " << c.rho;
	  getchar();
	}
	double cr = sqrt(gam*c.p/c.rho);
	Cellule& cg = grille[i][j-1][k];
	Cellule& cg2 = grille[i][j-2][k];
	double drho = (c.rho-cg2.rho)/dy/2.;
	double du = (c.u-cg2.u)/dy/2.;
	double dv = (c.v-cg2.v)/dy/2.;
	double dw = (c.w-cg2.w)/dy/2.;
	double dp = (c.p-cg2.p)/dy/2.;
	double L = min(min(domainex,domainey),domainez);
	double k = 0.278;
	double alpha = k*(cr*cr-c.v*c.v)/L/cr;
	double pinf = P(c.x,c.y,c.z,dx,dy,dz);
	double L0 = c.v*(cr*cr*drho-dp);
	double L1 = (c.v+cr)*(c.rho*cr*dv+dp);
	double L2 = c.rho*c.v*cr*dw;
	double L3 = c.rho*c.v*cr*du;
	double L4 = alpha*(c.p-pinf);
	if(c.v>0. && c.v-cr<0.){
	  //Subsonic outflow
	  L0 = c.v*(cr*cr*drho-dp);
	  L1 = (c.v+cr)*(c.rho*cr*dv+dp);
	  L2 = c.rho*c.v*cr*dw;
	  L3 = c.rho*c.v*cr*du;
	  L4 = alpha*(c.p-pinf);
	}
	else if(c.v>0. && c.v-cr>0.){
	  //Supersonic outflow
	  L0 = c.v*(cr*cr*drho-dp);
	  L1 = (c.v+cr)*(c.rho*cr*dv+dp);
	  L2 = c.rho*c.v*cr*dw;
	  L3 = c.rho*c.v*cr*du;
	  L4 = (c.v-cr)*(-c.rho*cr*dv+dp);
	}
	else if(c.v<0. && c.v+cr>0.){
	  //Subsonic inflow
	  L0 = 0.;
	  L1 = alpha*(c.p-pinf);
	  L2 = 0.;
	  L3 = 0.;
	  L4 = L1;
	}
	else if(c.v<0. && c.v+cr<0.){
	  //Supersonic inflow
	  L0 = 0.;
	  L1 = alpha*(c.p-pinf);
	  L2 = 0.;
	  L3 = 0.;
	  L4 = L1;
	}
	double d0 = (L0+L1/2.+L4/2.)/cr/cr;
	double d1 = (L1-L4)/2./c.rho/cr;
	double d2 = L2/c.rho/cr;
	double d3 = L3/c.rho/cr;
	double d4 = (L1+L4)/2.;
	d0 *= dy;
	d1 *= dy;
	d2 *= dy;
	d3 *= dy;
	d4 *= dy;
	double kappa = 0.05;
	c.fluxj[0] = cg.fluxj[0]+d0;
	c.fluxj[2] = cg.fluxj[2]+(c.v*d0+c.rho*d1);
	c.fluxj[3] = cg.fluxj[3]+(c.w*d0+c.rho*d2);
	c.fluxj[1] = cg.fluxj[1]+(c.u*d0+c.rho*d3);
	c.fluxj[4] = cg.fluxj[4]+((c.u*c.u+c.v*c.v+c.w*c.w)/2.*d0+c.rho*c.v*d1+c.rho*c.w*d2+c.rho*c.u*d3+d4/(gam-1.));
      }
    }
  }
  if(BC_y_in==3){
    int j=marge-1;
    for(int k=0;k<Nz+2*marge;k++){
      for(int i=0;i<Nx+2*marge;i++){
	Cellule& c = grille[i][j][k];
	if(c.p<eps || c.rho<eps){
	  cout << "p or rho negative: p " << c.p << " rho " << c.rho;
	  getchar();
	}
	Cellule& cd = grille[i][j+1][k];
	Cellule& cd2 = grille[i][j+2][k];
	double cr = sqrt(gam*cd.p/cd.rho);
	double drho = (cd2.rho-c.rho)/dy/2.;
	double du = (cd2.u-c.u)/dy/2.;
	double dv = (cd2.v-c.v)/dy/2.;
	double dw = (cd2.w-c.w)/dy/2.;
	double dp = (cd2.p-c.p)/dy/2.;
	double L = min(min(domainex,domainey),domainez);
	double k = 0.278;
	double alpha = k*(cr*cr-c.v*c.v)/L/cr;
	double pinf = P(cd.x,cd.y,cd.z,dx,dy,dz);
	double L0 = cd.v*(cr*cr*drho-dp);
	double L4 = (cd.v-cr)*(-cd.rho*cr*dv+dp);
	double L2 = cd.rho*cd.v*cr*dw;
	double L3 = cd.rho*cd.v*cr*du;
	double L1 = alpha*(cd.p-pinf);
	if(cd.v<0. && c.v-cr>0.){
	  //Subsonic outflow
	  L0 = cd.v*(cr*cr*drho-dp);
	  L4 = (cd.v-cr)*(-cd.rho*cr*dv+dp);
	  L2 = cd.rho*cd.v*cr*dw;
	  L3 = cd.rho*cd.v*cr*du;
	  L1 = alpha*(cd.p-pinf);
	}
	else if(cd.v<0. && cd.v-cr<0.){
	  //Supersonic outflow
	  L0 = cd.v*(cr*cr*drho-dp);
	  L4 = (cd.v-cr)*(-cd.rho*cr*dv+dp);
	  L2 = cd.rho*cd.v*cr*dw;
	  L3 = cd.rho*cd.v*cr*du;
	  L4 = (cd.v+cr)*(cd.rho*cr*dv+dp);
	}
	else if(cd.v>0. && cd.v-cr<0.){
	  //Subsonic inflow
	  L0 = 0.;
	  L4 = alpha*(cd.p-pinf);
	  L2 = 0.;
	  L3 = 0.;
	  L1 = L4;
	}
	else if(cd.v>0. && cd.v-cr>0.){
	  //Supersonic inflow
	  L0 = 0.;
	  L4 = alpha*(c.p-pinf);
	  L2 = 0.;
	  L3 = 0.;
	  L1 = L4;
	}
	double d0 = (L0+L1/2.+L4/2.)/cr/cr;
	double d1 = (L1-L4)/2./cd.rho/cr;
	double d2 = L2/cd.rho/cr;
	double d3 = L3/cd.rho/cr;
	double d4 = (L1+L4)/2.;
	d0 *= dy;
	d1 *= dy;
	d2 *= dy;
	d3 *= dy;
	d4 *= dy;
	double kappa = 0.05;
	c.fluxj[0] = cd.fluxj[0]-d0;
	c.fluxj[2] = cd.fluxj[2]-(cd.v*d0+cd.rho*d1);
	c.fluxj[3] = cd.fluxj[3]-(cd.w*d0+cd.rho*d2);
	c.fluxj[1] = cd.fluxj[1]-(cd.u*d0+cd.rho*d3);
	c.fluxj[4] = cd.fluxj[4]-((cd.u*cd.u+cd.v*cd.v+cd.w*cd.w)/2.*d0+cd.rho*cd.v*d1+cd.rho*cd.w*d2+cd.rho*cd.u*d3+d4/(gam-1.));
      }
    }
  }
}



/*!\brief Computation of the numerical z-flux.
   \param sigma = \a dt/dz: time-step/fluid spatial discretization step in the x direction
   \param t current simulation time
   \return void
*/
void Grille::fnumz(const double sigma, double t){ 
  //Initialization to the centered flux
  for(int i=0; i<Nx+2*marge-1; i++){
    for(int j=0; j<Ny+2*marge-1; j++){
      for(int k=0; k<Nz+2*marge-1; k++){ 
	Cellule&  c = grille[i][j][k]; 
	Cellule&  ck = grille[i][j][k+1];
	if(!c.vide && !ck.vide){
	  //Computation of the order indicators 
	  for(int l=0;l< c.ordre;l++){ 
	    c.co[l]=1.; 
	  } 
	  for(int l=c.ordre;l<ordremax;l++){ 
	    c.co[l]=0.;
	  }  
                
	  //Centered flux part  
	  c.fluxk[0] = (c.impz+ck.impz)/2.; 
	  c.fluxk[1] = (c.rho*c.u*c.w+ck.rho*ck.u*ck.w)/2.;
	  c.fluxk[2] = (c.rho*c.v*c.w+ck.rho*ck.v*ck.w)/2.;
	  c.fluxk[3] = (c.rho*c.w*c.w +c.p +ck.rho*ck.w*ck.w + ck.p)/2.;
	  c.fluxk[4] = ((c.rhoE+c.p)*c.w+(ck.rhoE+ck.p)*ck.w)/2.; 
	}
	else if(std::abs(c.alpha-1.)>eps){
	  //Lax-Friedrichs flux near void
	  c.fluxk[0] = (c.impz+ck.impz)/2. + (c.rho -ck.rho)/2./sigma; 
	  c.fluxk[1] = (c.rho*c.u*c.w+ck.rho*ck.u*ck.w)/2. + (c.impx -ck.impx)/2./sigma; 
	  c.fluxk[2] = (c.rho*c.v*c.w  +ck.rho*ck.v*ck.w)/2. + (c.impy -ck.impy)/2./sigma; 
	  c.fluxk[3] = (c.rho*c.w*c.w +c.p+ ck.rho*ck.w*ck.w +ck.p)/2. + (c.impz -ck.impz)/2./sigma; 
	  c.fluxk[4] = ((c.rhoE+c.p)*c.u+(ck.rhoE+ck.p)*ck.u)/2. + (c.rhoE -ck.rhoE)/2./sigma; 
	}
	else {
	  c.fluxk[0] = 0.; 
	  c.fluxk[1] = 0.; 
	  c.fluxk[2] = 0.; 
	  c.fluxk[3] = 0.; 
	  c.fluxk[4] = 0.; 
	}
      }
    }
  } 
    
  //Computation of the preliminary limiter variables 
    
  for(int i=0; i<Nx+2*marge-1; i++){
    for(int j=0; j<Ny+2*marge-1; j++){
      for(int k=0; k<Nz+2*marge-1; k++){ 
	Cellule& ck = grille[i][j][k+1]; 
	Cellule& c = grille[i][j][k]; 
	if(!c.vide && !ck.vide){ 
	  //Computation of the Roe variables 
	  double roe = sqrt(ck.rho/c.rho); 
	  double rhor = roe*c.rho; 
	  double ur = (roe*ck.u+c.u)/(1.+roe); 
	  double vr = (roe*ck.v+c.v)/(1.+roe); 
	  double wr = (roe*ck.w+c.w)/(1.+roe); 
	  double Hr = (roe*(ck.rho*ck.u*ck.u/2.+ ck.rho*ck.v*ck.v/2.+ ck.rho*ck.w*ck.w/2. 
			    + ck.p*gam/(gam-1.))/ck.rho+(c.rho*c.u*c.u/2.+c.rho*c.v*c.v/2.+ c.rho*c.w*c.w/2.
							 + c.p*gam/(gam-1.))/c.rho)/(1.+roe);
	  double ur2 = ur*ur; 
	  double vr2 = vr*vr;
	  double wr2 = wr*wr;
	  double cr2 = (gam-1.)*(Hr-ur2/2.-vr2/2.-wr2/2.); 
									
	  //Test on the speed of sound
	  if(cr2<=0. && abs(c.alpha-1.)>eps){
	    cout << "COmputation of the z-flux" << endl;
	    cout << "i=" << i << " j=" << j << " negative speed of sound: c2=" << cr2 << endl;
	    cout << "x=" << c.x << " y=" << c.y << " z=" << c.z << " alpha=" << c.alpha << endl;
	    cout << "t=" << t << endl; 
	    cout << "c.p=" << c.p << endl; 
	    cout << "c.rho=" << c.rho << endl; 
	    cout << "c.u=" << c.u << endl; 
	    cout << "c.v=" << c.v << endl; 
	    cout << "c.w=" << c.w << endl; 
	    cout << "ck.p=" << ck.p << endl; 
	    cout << "ck.rho=" << ck.rho << endl; 
	    cout << "ck.u=" << ck.u << endl;
	    cout << "ck.v=" << ck.v << endl;
	    cout << "ck.w=" << ck.w << endl;
	    cout << "ur=" << ur << endl; 
	    cout << "ur2=" << ur2 << endl; 
	    cout << "vr=" << vr << endl; 
	    cout << "vr2=" << vr2 << endl; 
	    cout << "wr2=" << wr2 << endl; 
	    cout << "Hr=" << Hr << endl;
	    getchar(); 
	  } 
									
	  //Speed of sound 
	  double cr = sqrt(cr2); 
									
	  //Eigenvalues 
	  c.lambda[0] = wr-cr; 
	  c.lambda[1] = wr;
	  c.lambda[2] = wr;
	  c.lambda[3] = wr;
	  c.lambda[4] = wr+cr; 
									
	  //Computation of the differences between Wd and Wg 
	  double drho = ck.rho - c.rho; 
	  double du = ck.u - c.u;
	  double dv = ck.v - c.v;
	  double dw = ck.w - c.w;
	  double dp = ck.p - c.p; 
									
	  //Computation of the deltaV (differences between Wd and Wg in the eigenvectors basis) 
	  double ros2c = rhor/cr/2.; 
	  c.delw[0] = dp/cr2/2. - ros2c*dw; 
	  c.delw[1] = drho - dp/cr2;
	  c.delw[2] = 2.*ros2c*du;
	  c.delw[3] = 2.*ros2c*dv;   // a verifie!!!!!
	  c.delw[4] = dp/cr2/2. + ros2c*dw; 
									
	  //Computation of the complete correction in the eigenvectors basis
	  double xnu[5]; 
	  for(int l=0;l<5;l++){ 
	    xnu[l]  = sigma*abs(c.lambda[l]); 
	    c.delwnu[l] = abs(c.lambda[l])*(1.-xnu[l])*c.delw[l]; 
	    //Higher-order corrective coefficients
	    c.cf2[l]  = c.co[1]*abs(c.lambda[l])*(1.-xnu[l]); 
	    c.cf3[l]  = c.co[2]*c.cf2[l]*(1.+xnu[l])/3.; 
	    c.cf4[l]  = c.co[3]*c.cf3[l]*(xnu[l]-2.)/4.; 
	    c.cf5[l]  = c.co[4]*c.cf4[l]*(xnu[l]+2.)/5.; 
	    c.cf6[l]  = c.co[5]*c.cf5[l]*(xnu[l]-3.)/6.; 
	    c.cf7[l]  = c.co[6]*c.cf6[l]*(xnu[l]+3.)/7.; 
	    c.cf8[l]  = c.co[7]*c.cf7[l]*(xnu[l]-4.)/8.; 
	    c.cf9[l]  = c.co[8]*c.cf8[l]*(xnu[l]+4.)/9.; 
	    c.cf10[l] = c.co[9]*c.cf9[l]*(xnu[l]-5.)/10.; 
	    c.cf11[l] = c.co[10]*c.cf10[l]*(xnu[l]+5.)/11.; 
	  } 
									
									
	  for(int l=0;l<5;l++){ 
	    //Centered corrections
	    c.psic0[l] = (c.cf2[l]-2.*c.cf4[l]+6.*c.cf6[l]-20.*c.cf8[l]+70.*c.cf10[l])*c.delw[l]; 
	    c.psic1[l] = (c.cf4[l]-4.*c.cf6[l]+15.*c.cf8[l]-56.*c.cf10[l])*c.delw[l]; 
	    c.psic2[l] = (c.cf6[l]-6.*c.cf8[l]+28.*c.cf10[l])*c.delw[l]; 
	    c.psic3[l] = (c.cf8[l]-8.*c.cf10[l])*c.delw[l]; 
	    c.psic4[l] = (c.cf10[l])*c.delw[l]; 
	    //Decentered corrections 
	    c.psid0[l] = (126.*c.cf11[l]-35.*c.cf9[l]+10.*c.cf7[l]-3.*c.cf5[l]+c.cf3[l])*c.delw[l]; 
	    c.psid1[l] = (84.*c.cf11[l]-21.*c.cf9[l]+5.*c.cf7[l]-c.cf5[l])*c.delw[l]; 
	    c.psid2[l] = (36.*c.cf11[l]-7.*c.cf9[l]+c.cf7[l])*c.delw[l]; 
	    c.psid3[l] = (9.*c.cf11[l]-c.cf9[l])*c.delw[l]; 
	    c.psid4[l] = (c.cf11[l])*c.delw[l]; 
	  } 
									
	  //Left eigenvectors
	  c.vpr[0][0] = 1.; 
	  c.vpr[1][0] = ur; 
	  c.vpr[2][0] = vr;
	  c.vpr[3][0] = wr-cr;
	  c.vpr[4][0] = Hr-wr*cr; 
									
	  c.vpr[0][1] = 1.; 
	  c.vpr[1][1] = ur;
	  c.vpr[2][1] = vr;
	  c.vpr[3][1] = wr;
	  c.vpr[4][1] = ur2/2. + vr2/2. + wr2/2.;
									
	  c.vpr[0][2] = 0.; 
	  c.vpr[1][2] = 0.;
	  c.vpr[2][2] = cr;
	  c.vpr[3][2] = 0.;
	  c.vpr[4][2] = ur*cr;
									
	  c.vpr[0][3] = 0.; 
	  c.vpr[1][3] = 0.;
	  c.vpr[2][3] = 0;
	  c.vpr[3][3] = cr;
	  c.vpr[4][3] = vr*cr;
									
	  c.vpr[0][4] = 1.; 
	  c.vpr[1][4] = ur;
	  c.vpr[2][4] = vr;
	  c.vpr[3][4] = wr+cr;
	  c.vpr[4][4] = Hr+wr*cr; 
									
	  //Corrections in the eigenvectors basis
	  for(int l=0;l<5;l++){ 
	    c.psic0r[l] = 0.; 
	    c.psic1r[l] = 0.; 
	    c.psic2r[l] = 0.; 
	    c.psic3r[l] = 0.; 
	    c.psic4r[l] = 0.; 
	    c.psid0r[l] = 0.; 
	    c.psid1r[l] = 0.; 
	    c.psid2r[l] = 0.; 
	    c.psid3r[l] = 0.; 
	    c.psid4r[l] = 0.; 
	  } 
	  for(int m=0;m<5;m++){ 
	    for(int l=0;l<5;l++){ 
	      c.psic0r[m] += c.psic0[l]*c.vpr[m][l]; 
	      c.psic1r[m] += c.psic1[l]*c.vpr[m][l]; 
	      c.psic2r[m] += c.psic2[l]*c.vpr[m][l]; 
	      c.psic3r[m] += c.psic3[l]*c.vpr[m][l]; 
	      c.psic4r[m] += c.psic4[l]*c.vpr[m][l]; 
	      c.psid0r[m] += c.psid0[l]*c.vpr[m][l]; 
	      c.psid1r[m] += c.psid1[l]*c.vpr[m][l]; 
	      c.psid2r[m] += c.psid2[l]*c.vpr[m][l]; 
	      c.psid3r[m] += c.psid3[l]*c.vpr[m][l]; 
	      c.psid4r[m] += c.psid4[l]*c.vpr[m][l]; 
	    } 
	  } 
	}
      }
    }
  }  
    
  //Monotonicity indicators
  for(int l=0;l<5;l++){ 
    for(int i=1;i<Nx+2*marge-1;i++){
      for(int j=1;j<Ny+2*marge-1;j++){
	for(int k=1;k<Nz+2*marge-1;k++){ 
	  Cellule& c = grille[i][j][k]; 
	  Cellule& cg = grille[i][j][k-1]; 
	  if(!c.vide && !cg.vide){
	    c.am[l] = c.lambda[l]*c.delw[l]-cg.lambda[l]*cg.delw[l]; 
	  }
	}
      }
    } 
    //Computation of dj^m4 
    for(int i=0;i<Nx+2*marge;i++){
      for(int j=0;j<Ny+2*marge;j++){   
	for(int k=1;k<Nz+2*marge-2;k++){ 
	  Cellule& c = grille[i][j][k]; 
	  Cellule& cd = grille[i][j][k+1]; 
	  if(!c.vide && !cd.vide){
	    double z1 = 4.*c.am[l]-cd.am[l]; 
	    double z2 = 4.*cd.am[l]-c.am[l]; 
	    double z3 = c.am[l]; 
	    double z4 = cd.am[l]; 
	    c.am1[l] = (sign(z1)+sign(z2))/2.*abs((sign(z1)+sign(z3))/2.)*(sign(z1)
									   + sign(z4))/2.*min(abs(z1),min(abs(z2),min(abs(z3),abs(z4)))); 
	  }
	}
      }
    } 
  } 
    
  //Computation of r+ and r- 
  for(int l=0;l<5;l++){ 
    for(int i=marge;i<Nx+2*marge-4;i++){ 
      for(int j=marge;j<Ny+2*marge-4;j++){
	for(int k=marge;k<Nz+2*marge-4;k++){ 
	  Cellule& c = grille[i][j][k]; 
	  Cellule& cd = grille[i][j][k+1]; 
	  Cellule&  cg = grille[i][j][k-1];
	  if(!c.vide && !cd.vide){ 
	    c.rp[l] = sign(c.delw[l])*sign(cg.delw[l])*(abs(cg.delw[l])+eps)/(abs(c.delw[l])+eps); 
	    c.rm[l] = sign(c.delw[l])*sign(cd.delw[l])*(abs(cd.delw[l])+eps)/(abs(c.delw[l])+eps); 
	    //Higher-order corrections 
	    Cellule& cg2 = grille[i][j][k-2]; 
	    Cellule& cg3 = grille[i][j][k-3]; 
	    Cellule& cg4 = grille[i][j][k-4]; 
	    Cellule& cg5 = grille[i][j][k-5]; 
	    Cellule& cd2 = grille[i][j][k+2]; 
	    Cellule& cd3 = grille[i][j][k+3]; 
	    Cellule& cd4 = grille[i][j][k+4]; 
	    c.psid[l] = -c.psid0[l]+cg.psid0[l]+cd.psid1[l]-cg2.psid1[l]-cd2.psid2[l]
	      + cg3.psid2[l]+cd3.psid3[l]-cg4.psid3[l]-cd4.psid4[l]+cg5.psid4[l]; 
	  }
	}
      }
    } 
  } 
    
  //Flux computation 
  for(int i=marge-1;i<Nx+marge;i++){ 
    for(int j=marge-1;j<Ny+marge;j++){ 
      for(int k=marge-1;k<Nz+marge;k++){ 
	Cellule& c = grille[i][j][k]; 
	//Neighbouring cells
	Cellule& cg = grille[i][j][k-1]; 
	Cellule& cg2 = grille[i][j][k-2]; 
	Cellule& cg3 = grille[i][j][k-3]; 
	Cellule& cg4 = grille[i][j][k-4]; 
	Cellule& cd = grille[i][j][k+1]; 
	Cellule& cd2 = grille[i][j][k+2]; 
	Cellule& cd3 = grille[i][j][k+3]; 
	Cellule& cd4 = grille[i][j][k+4]; 
	if(!c.vide && !cd.vide){   
	  //TVD flux 
	  double tvd[5]; 
	  double psict[5]; 
									
	  //Initialization 
	  for(int l=0;l<5;l++){ 
	    tvd[l] = 0.; 
	    //Centered part
	    psict[l] = c.psic0r[l] + cg.psic1r[l] + cd.psic1r[l] + cg2.psic2r[l] + cd2.psic2r[l] 
	      + cg3.psic3r[l] + cd3.psic3r[l] + cg4.psic4r[l] + cd4.psic4r[l]; 
	  } 
									
	  //Limiter 
	  double psic; 
	  for(int l=0;l<5;l++){ 
	    psic = c.psic0[l] + cg.psic1[l] + cd.psic1[l] + cg2.psic2[l] 
	      + cd2.psic2[l] + cg3.psic3[l] + cd3.psic3[l] + cg4.psic4[l] + cd4.psic4[l]; 
											
	    //Decentered part
	    double r; 
	    double xnum; 
	    double xnume; 
	    int is; 
	    double psi; 
	    if(c.lambda[l]>0.){ 
	      r = c.rp[l]; 
	      xnum = sigma*abs(cg.lambda[l]); 
	      xnume = max(xnum,eps); 
	      is = 1; 
	      psi = psic+c.psid[l]; 
	    } else { 
	      r = c.rm[l]; 
	      xnum = sigma*abs(cd.lambda[l]); 
	      xnume = max(xnum,eps); 
	      is = -1; 
	      psi = psic-cd.psid[l]; 
	    } 
											
	    double xnu = sigma*abs(c.lambda[l]); 
	    xnu = max(xnu,eps); 
											
	    //TVD limiter psitvd 
	    psi = (double) sign(c.delwnu[l])*psi/(abs(c.delwnu[l]+eps)); 
	    double psimax1 = 2.*r*(1.-xnume)/(xnu*(1.-xnu)); 
	    double psimax2 = 2./(1.-xnu); 
	    double psitvd = max(0.,min(psi,min(psimax1,psimax2))); 
											
	    //Monotonicity criterion
	    if((c.delwnu[l] != 0.) && (abs(psi-psitvd)>eps)){ 
	      double dfo = psi*c.delwnu[l]/2.; 
	      double dabsf = psimax2*c.delwnu[l]/2.; 
	      double dful = psimax1*c.delwnu[l]/2.; 
	      double dfmd = dabsf/2.-c.am1[l]/2.; 
	      Cellule& camont = grille[i][j][k-is];   //Upwind cell
	      double dflc = dful/2.+((1.-xnume)/xnu)*camont.am1[l]/2.; 
	      double dfmin = max(min(0.,min(dabsf,dfmd)),min(0.,min(dful,dflc))); 
	      double dfmax = min(max(0.,max(dabsf,dfmd)),max(0.,max(dful,dflc))); 
	      if((dfmin-dfo)*(dfmax-dfo)>0.){ 
		psi = psitvd; 
	      } 
	    } 
											
	    //Uncomment to use only TVD without MP 
	    //psi = psitvd; 
											
	    //Uncomment to use the scheme without TVD nor MP 
	    //psi = 0.; 
											
	    double ctvd = psi*c.delwnu[l]/2.-abs(c.lambda[l])*c.delw[l]/2.; 
	    for(int m=0;m<5;m++){ 
	      tvd[m] += ctvd*c.vpr[m][l];
	    } 
	  } 
									
	  //Final computation of the flux
	  for(int l=0;l<5;l++){ 
	    c.fluxk[l] += tvd[l]; 
	  }
	}
      }
    }
  } 
    
  //Entropy correction
  corentz(sigma);

  //Boundary conditions
  if(!flag_2d){

    //Reflecting boundary conditions
    if(BC_z_in ==  1 || BC_z_out ==  1){
      for(int i=0;i<Nx+2*marge;i++){
	for(int j=0;j<Ny+2*marge;j++){
                
	  if(BC_z_in ==  1){
	    Cellule& c = grille[i][j][marge-1];
	    Cellule& cp = grille[i][j][marge];
	    double p0 = (gam-1.)*(cp.rhoE0 - 1./2.*(cp.impx0*cp.impx0 + cp.impy0*cp.impy0 
						    + cp.impz0*cp.impz0)/cp.rho0);
	    c.fluxk[0] = 0.;
	    c.fluxk[1] = 0.;
	    c.fluxk[2] = 0.;
	    c.fluxk[3] = p0;
	    c.fluxk[4] = 0.;
  	  }
                
	  if(BC_z_out ==  1){
	    Cellule& c2 = grille[i][j][Nz+marge-1]; 
	    double p02 = (gam-1.)*(c2.rhoE0-1./2.*(c2.impx0*c2.impx0 + c2.impy0*c2.impy0 
						   + c2.impz0*c2.impz0)/c2.rho0);
	    c2.fluxk[0] = 0.;
	    c2.fluxk[1] = 0.;
	    c2.fluxk[2] = 0.;
	    c2.fluxk[3] = p02;
	    c2.fluxk[4] = 0.;
	  }
	} 
      }
    }
    
    //Periodic boundary conditions
    if(BC_z_in ==  2 || BC_z_out ==  2){
      for(int i=0;i<Nx+2*marge;i++){
	for(int j=0;j<Ny+2*marge;j++){
	  for(int k=0;k<Nz+2*marge;k++){
	    if(k==Nz+marge-1){
	      grille[i][j][k].fluxk[0] = grille[i][j][marge-1].fluxk[0];
	      grille[i][j][k].fluxk[1] = grille[i][j][marge-1].fluxk[1];
	      grille[i][j][k].fluxk[2] = grille[i][j][marge-1].fluxk[2];
	      grille[i][j][k].fluxk[3] = grille[i][j][marge-1].fluxk[3];
	      grille[i][j][k].fluxk[4] = grille[i][j][marge-1].fluxk[4];
	    }
	  }
	}
      }
    }  
  }	

  //Outflow boundary conditions (Poinsot-Lele)
  if(BC_z_out==3){
    int k=Nz+marge-1;
    for(int i=0;i<Nx+2*marge;i++){
      for(int j=0;j<Ny+2*marge;j++){
	Cellule& c = grille[i][j][k];
	if(c.p<eps || c.rho<eps){
	  cout << "p or rho negative: p " << c.p << " rho " << c.rho;
	  getchar();
	}
	double cr = sqrt(gam*c.p/c.rho);
	Cellule& cg = grille[i][j][k-1];
	Cellule& cg2 = grille[i][j][k-2];
	double drho = (c.rho-cg2.rho)/dz/2.;
	double du = (c.u-cg2.u)/dz/2.;
	double dv = (c.v-cg2.v)/dz/2.;
	double dw = (c.w-cg2.w)/dz/2.;
	double dp = (c.p-cg2.p)/dz/2.;
	double L = min(min(domainex,domainey),domainez);
	double k = 0.278;
	double alpha = k*(cr*cr-c.w*c.w)/L/cr;
	double pinf = P(c.x,c.y,c.z,dx,dy,dz);
	double L0 = c.w*(cr*cr*drho-dp);
	double L1 = (c.w+cr)*(c.rho*cr*dw+dp);
	double L2 = c.rho*c.w*cr*du;
	double L3 = c.rho*c.w*cr*dv;
	double L4 = alpha*(c.p-pinf);
	if(c.w>0. && c.w-cr<0.){
	  //Subsonic outflow
	  L0 = c.w*(cr*cr*drho-dp);
	  L1 = (c.w+cr)*(c.rho*cr*dw+dp);
	  L2 = c.rho*c.w*cr*du;
	  L3 = c.rho*c.w*cr*dv;
	  L4 = alpha*(c.p-pinf);
	}
	else if(c.w>0. && c.w-cr>0.){
	  //Supersonic outflow
	  L0 = c.w*(cr*cr*drho-dp);
	  L1 = (c.w+cr)*(c.rho*cr*dw+dp);
	  L2 = c.rho*c.w*cr*du;
	  L3 = c.rho*c.w*cr*dv;
	  L4 = (c.w-cr)*(-c.rho*cr*dw+dp);
	}
	else if(c.w<0. && c.w+cr>0.){
	  //Subsonic inflow
	  L0 = 0.;
	  L1 = alpha*(c.p-pinf);
	  L2 = 0.;
	  L3 = 0.;
	  L4 = L1;
	}
	else if(c.w<0. && c.w+cr<0.){
	  //Supersonic inflow
	  L0 = 0.;
	  L1 = alpha*(c.p-pinf);
	  L2 = 0.;
	  L3 = 0.;
	  L4 = L1;
	}
	double d0 = (L0+L1/2.+L4/2.)/cr/cr;
	double d1 = (L1-L4)/2./c.rho/cr;
	double d2 = L2/c.rho/cr;
	double d3 = L3/c.rho/cr;
	double d4 = (L1+L4)/2.;
	d0 *= dz;
	d1 *= dz;
	d2 *= dz;
	d3 *= dz;
	d4 *= dz;
	double kappa = 0.05;
	c.fluxk[0] = cg.fluxk[0]+d0;
	c.fluxk[3] = cg.fluxk[3]+(c.w*d0+c.rho*d1);
	c.fluxk[1] = cg.fluxk[1]+(c.u*d0+c.rho*d2);
	c.fluxk[2] = cg.fluxk[2]+(c.v*d0+c.rho*d3);
	c.fluxk[4] = cg.fluxk[4]+((c.u*c.u+c.v*c.v+c.w*c.w)/2.*d0+c.rho*c.w*d1+c.rho*c.u*d2+c.rho*c.v*d3+d4/(gam-1.));
      }
    }
  }
  if(BC_z_in==3){
    int k=marge-1;
    for(int i=0;i<Nx+2*marge;i++){
      for(int j=0;j<Ny+2*marge;j++){
	Cellule& c = grille[i][j][k];
	if(c.p<eps || c.rho<eps){
	  cout << "p or rho negative: p " << c.p << " rho " << c.rho;
	  getchar();
	}
	Cellule& cd = grille[i][j][k+1];
	Cellule& cd2 = grille[i][j][k+2];
	double cr = sqrt(gam*cd.p/cd.rho);
	double drho = (cd2.rho-c.rho)/dz/2.;
	double du = (cd2.u-c.u)/dz/2.;
	double dv = (cd2.v-c.v)/dz/2.;
	double dw = (cd2.w-c.w)/dz/2.;
	double dp = (cd2.p-c.p)/dz/2.;
	double L = min(min(domainex,domainey),domainez);
	double k = 0.278;
	double alpha = k*(cr*cr-c.w*c.w)/L/cr;
	double pinf = P(cd.x,cd.y,cd.z,dx,dy,dz);
	double L0 = cd.w*(cr*cr*drho-dp);
	double L4 = (cd.w-cr)*(-cd.rho*cr*dw+dp);
	double L2 = cd.rho*cd.w*cr*du;
	double L3 = cd.rho*cd.w*cr*dv;
	double L1 = alpha*(cd.p-pinf);
	if(cd.w<0. && c.w-cr>0.){
	  //Subsonic outflow
	  L0 = cd.w*(cr*cr*drho-dp);
	  L4 = (cd.w-cr)*(-cd.rho*cr*dw+dp);
	  L2 = cd.rho*cd.w*cr*du;
	  L3 = cd.rho*cd.w*cr*dv;
	  L1 = alpha*(cd.p-pinf);
	}
	else if(cd.w<0. && cd.w-cr<0.){
	  //Supersonic outflow
	  L0 = cd.w*(cr*cr*drho-dp);
	  L4 = (cd.w-cr)*(-cd.rho*cr*dw+dp);
	  L2 = cd.rho*cd.w*cr*du;
	  L3 = cd.rho*cd.w*cr*dv;
	  L4 = (cd.w+cr)*(cd.rho*cr*dw+dp);
	}
	else if(cd.w>0. && cd.w-cr<0.){
	  //Subsonic inflow
	  L0 = 0.;
	  L4 = alpha*(cd.p-pinf);
	  L2 = 0.;
	  L3 = 0.;
	  L1 = L4;
	}
	else if(cd.w>0. && cd.w-cr>0.){
	  //Supersonic inflow
	  L0 = 0.;
	  L4 = alpha*(c.p-pinf);
	  L2 = 0.;
	  L3 = 0.;
	  L1 = L4;
	}
	double d0 = (L0+L1/2.+L4/2.)/cr/cr;
	double d1 = (L1-L4)/2./cd.rho/cr;
	double d2 = L2/cd.rho/cr;
	double d3 = L3/cd.rho/cr;
	double d4 = (L1+L4)/2.;
	d0 *= dz;
	d1 *= dz;
	d2 *= dz;
	d3 *= dz;
	d4 *= dz;
	double kappa = 0.05;
	c.fluxk[0] = cd.fluxk[0]-d0;
	c.fluxk[1] = cd.fluxk[1]-(cd.w*d0+cd.rho*d1);
	c.fluxk[2] = cd.fluxk[2]-(cd.u*d0+cd.rho*d2);
	c.fluxk[3] = cd.fluxk[3]-(cd.v*d0+cd.rho*d3);
	c.fluxk[4] = cd.fluxk[4]-((cd.u*cd.u+cd.v*cd.v+cd.w*cd.w)/2.*d0+cd.rho*cd.w*d1+cd.rho*cd.u*d2+cd.rho*cd.v*d3+d4/(gam-1.));
      }
    }
  }
}




/*!\brief Boundary conditions. 
   \details Types of BC:  1 = reflecting; 2 = periodic; 3= outflow.
   \return void
*/
void Grille::BC(){ 
  // Inner Boundary Condition for x
  for(int i=0;i<marge;i++){
    for(int j=0;j<Ny+2*marge;j++){
      for(int k=0;k<Nz+2*marge;k++){
	Cellule& c = grille[i][j][k];
	Cellule&  cm = grille[2*marge-i-1][j][k];   //Mirror cell
	Cellule& cp = grille[Nx+i][j][k];          //Periodic cell
	
	if(BC_x_in ==  1){
	  c.rho = cm.rho;
	  c.u   = -cm.u;
	  c.v = cm.v;
	  c.w = cm.w;
	  c.p = cm.p;
	  c.impx = -cm.impx;
	  c.impy = cm.impy;
	  c.impz = cm.impz;
	  c.rhoE = cm.rhoE;
	  c.vide = cm.vide;
	}
	if(BC_x_in ==  3){
	  c.rho = cm.rho;
	  c.u   = cm.u;
	  c.v = cm.v;
	  c.w = cm.w;
	  c.p = cm.p;
	  c.impx = cm.impx;
	  c.impy = cm.impy;
	  c.impz = cm.impz;
	  c.rhoE = cm.rhoE;
	  c.vide = cm.vide;
	}
	else if(BC_x_in ==  2){
	  c.rho = cp.rho;
	  c.u   = cp.u;
	  c.v = cp.v;
	  c.w = cp.w;
	  c.p   = cp.p;
	  c.impx = cp.impx;
	  c.impy = cp.impy;
	  c.impz = cp.impz;
	  c.rhoE= cp.rhoE;
	  c.vide = cp.vide;
	}
      }
    }
  }
    
  // Outer Boundary Condition for x
  for(int i=Nx+marge;i<Nx+2*marge;i++){
    for(int j=0;j<Ny+2*marge;j++){
      for(int k=0;k<Nz+2*marge;k++){
	Cellule&  c = grille[i][j][k];
	Cellule&  cm = grille[2*Nx+2*marge-i-1][j][k];  //Mirror cell
	Cellule&  cp = grille[i-Nx][j][k];              //Periodic cell
                
	if(BC_x_out == 1){
	  c.rho = cm.rho;
	  c.u = -cm.u;
	  c.v = cm.v;
	  c.w = cm.w;
	  c.p   = cm.p;
	  c.impx = -cm.impx;
	  c.impy = cm.impy;
	  c.impz = cm.impz;
	  c.rhoE= cm.rhoE;
	  c.vide = cm.vide;
	}
	if(BC_x_out == 3){
	  c.rho = cm.rho;
	  c.u = cm.u;
	  c.v = cm.v;
	  c.w = cm.w;
	  c.p   = cm.p;
	  c.impx = cm.impx;
	  c.impy = cm.impy;
	  c.impz = cm.impz;
	  c.rhoE= cm.rhoE;
	  c.vide = cm.vide;
	}
	else if(BC_x_out == 2){
	  c.rho = cp.rho;
	  c.u = cp.u;
	  c.v = cp.v;
	  c.w = cp.w;
	  c.p   = cp.p;
	  c.impx = cp.impx;
	  c.impy = cp.impy;
	  c.impz = cp.impz;
	  c.rhoE= cp.rhoE;
	  c.vide = cp.vide;
	}
      }
    }
  }
    
  // Inner Boundary Condition for y
  for(int i=0;i<Nx+2*marge;i++){
    for(int j=0;j<marge;j++){
      for(int k=0;k<Nz+2*marge;k++){
	Cellule& c = grille[i][j][k];
	Cellule& cm = grille[i][2*marge-j-1][k];             //Mirror cell
	Cellule& cp = grille[i][Ny+j][k];                   //Periodic cell
                
	if(BC_y_in == 1){
	  c.rho = cm.rho;
	  c.u   = cm.u;
	  c.v = -cm.v;
	  c.w   = cm.w;
	  c.p   = cm.p;
	  c.impx = cm.impx;
	  c.impy = -cm.impy;
	  c.impz = cm.impz;
	  c.rhoE= cm.rhoE;
	  c.vide = cm.vide;
	}
	if(BC_y_in == 3){
	  c.rho = cm.rho;
	  c.u   = cm.u;
	  c.v = cm.v;
	  c.w   = cm.w;
	  c.p   = cm.p;
	  c.impx = cm.impx;
	  c.impy = cm.impy;
	  c.impz = cm.impz;
	  c.rhoE= cm.rhoE;
	  c.vide = cm.vide;
	}
	else if(BC_y_in == 2){
	  c.rho = cp.rho;
	  c.u   = cp.u;
	  c.v = cp.v;
	  c.w = cp.w;
	  c.p   = cp.p;
	  c.impx = cp.impx;
	  c.impy = cp.impy;
	  c.impz = cp.impz;
	  c.rhoE= cp.rhoE;
	  c.vide = cp.vide;
	}
      }
    }
  }
    
  // Outer Boundary Condition for y
  for(int i=0;i<Nx+2*marge;i++){
    for(int j=Ny+marge;j<Ny+2*marge;j++){
      for(int k=0;k<Nz+2*marge;k++){
	Cellule&  c = grille[i][j][k];
	Cellule&  cm = grille[i][2*Ny+2*marge-j-1][k];      //Mirror cell
	Cellule&  cp = grille[i][j-Ny][k];                  //Priodic cell
                
	if(BC_y_out == 1){
	  c.rho = cm.rho;
	  c.u   = cm.u;
	  c.v = -cm.v;
	  c.w = cm.w;
	  c.p   = cm.p;
	  c.impx = cm.impx;
	  c.impy = -cm.impy;
	  c.impz = cm.impz;
	  c.rhoE= cm.rhoE;
	  c.vide = cm.vide;
	}
	if(BC_y_out == 3){
	  c.rho = cm.rho;
	  c.u   = cm.u;
	  c.v = cm.v;
	  c.w = cm.w;
	  c.p   = cm.p;
	  c.impx = cm.impx;
	  c.impy = cm.impy;
	  c.impz = cm.impz;
	  c.rhoE= cm.rhoE;
	  c.vide = cm.vide;
	}
	else if(BC_y_out == 2){
	  c.rho = cp.rho;
	  c.u = cp.u;
	  c.v = cp.v;
	  c.w = cp.w;
	  c.p   = cp.p;
	  c.impx = cp.impx;
	  c.impy = cp.impy;
	  c.impz = cp.impz;
	  c.rhoE= cp.rhoE;
	  c.vide = cp.vide;
	}
      }
    }
  }
    
  // Inner Boundary Condition for z
  if(!flag_2d){  
    for(int i=0;i<Nx+2*marge;i++){
      for(int j=0;j<Ny+2*marge;j++){
	for(int k=0;k<marge;k++){
	  Cellule&  c = grille[i][j][k];
	  Cellule&  cm = grille[i][j][2*marge-k-1];             //Mirror cell
	  Cellule&  cp = grille[i][j][Nz+k];                   //Periodic cell
                
	  if(BC_z_in == 1){
	    c.rho = cm.rho;
	    c.u = cm.u;
	    c.v = cm.v;
	    c.w = -cm.w;
	    c.p = cm.p;
	    c.impx = cm.impx;
	    c.impy = cm.impy;
	    c.impz = -cm.impz;
	    c.rhoE= cm.rhoE;
	    c.vide = cm.vide;
	  }
	  if(BC_z_in == 3){
	    c.rho = cm.rho;
	    c.u = cm.u;
	    c.v = cm.v;
	    c.w = cm.w;
	    c.p = cm.p;
	    c.impx = cm.impx;
	    c.impy = cm.impy;
	    c.impz = cm.impz;
	    c.rhoE= cm.rhoE;
	    c.vide = cm.vide;
	  }
	  else if(BC_z_in == 2){
	    c.rho = cp.rho;
	    c.u = cp.u;
	    c.v = cp.v;
	    c.w = cp.w;
	    c.p = cp.p;
	    c.impx = cp.impx;
	    c.impy = cp.impy;
	    c.impz = cp.impz;
	    c.rhoE= cp.rhoE;
	    c.vide = cp.vide;
	  }
	}
      }
    }
    
    // Outer Boundary Condition for z
    for(int i=0;i<Nx+2*marge;i++){
      for(int j=0;j<Ny+2*marge;j++){
	for(int k=Nz+marge;k<Nz+2*marge;k++){
	  Cellule&  c = grille[i][j][k];
	  Cellule&  cm = grille[i][j][2*Nz+2*marge-k-1];      //Mirror cell
	  Cellule& cp = grille[i][j][k-Nz];                  //Periodic cell
                
	  if(BC_z_out == 1){
	    c.rho = cm.rho;
	    c.u = cm.u;
	    c.v = cm.v;
	    c.w = -cm.w;
	    c.p = cm.p;
	    c.impx = cm.impx;
	    c.impy = cm.impy;
	    c.impz = -cm.impz;
	    c.rhoE= cm.rhoE;
	    c.vide = cm.vide;
	  }
	  if(BC_z_out == 3){
	    c.rho = cm.rho;
	    c.u = cm.u;
	    c.v = cm.v;
	    c.w = cm.w;
	    c.p = cm.p;
	    c.impx = cm.impx;
	    c.impy = cm.impy;
	    c.impz = cm.impz;
	    c.rhoE= cm.rhoE;
	    c.vide = cm.vide;
	  }
	  else if(BC_z_out == 2){
	    c.rho = cp.rho;
	    c.u = cp.u;
	    c.v = cp.v;
	    c.w = cp.w;
	    c.p   = cp.p;
	    c.impx = cp.impx;
	    c.impy = cp.impy;
	    c.impz = cp.impz;
	    c.rhoE= cp.rhoE;
	    c.vide = cp.vide;
	  }
	}
      }
    }
  }
	
  else{
    for(int i=0;i<Nx+2*marge;i++){
      for(int j=0;j<Ny+2*marge;j++){
	for(int k=0;k<marge;k++){
	  Cellule&  c = grille[i][j][k];
	  Cellule&  cm = grille[i][j][marge]; 
					
	  c.rho = cm.rho;
	  c.u = cm.u;
	  c.v = cm.v;
	  c.w = -cm.w;
	  c.p = cm.p;
	  c.impx = cm.impx;
	  c.impy = cm.impy;
	  c.impz = -cm.impz;
	  c.rhoE= cm.rhoE;
	  c.vide = cm.vide;

	}
      }
    }
		
    // Outer Boundary Condition for z
    for(int i=0;i<Nx+2*marge;i++){
      for(int j=0;j<Ny+2*marge;j++){
	for(int k=Nz+marge;k<Nz+2*marge;k++){
	  Cellule&  c = grille[i][j][k];
	  Cellule&  cm = grille[i][j][marge];      //Mirror cell

	  c.rho = cm.rho;
	  c.u = cm.u;
	  c.v = cm.v;
	  c.w = -cm.w;
	  c.p = cm.p;
	  c.impx = cm.impx;
	  c.impy = cm.impy;
	  c.impz = -cm.impz;
	  c.rhoE= cm.rhoE;
	  c.vide = cm.vide;
	}
      }
    }
  }
} 


/*!\brief Print out results. 
   \param n index of the output file
   \return void
*/
void Grille::Impression(int n){
  //Output of the vtk file
  std::ostringstream oss;
  oss << "resultats/fluide" << n << ".vtk";
  string s = oss.str();
  const char* const fluidevtk = s.c_str();
    
  //Open the output flux
  std::ofstream vtk(fluidevtk,ios::out);
  if(!vtk){
    cout <<"Opening of fluide" << n << ".vtk failed" << endl;
  }
  //Initialization of the vtk file
  vtk << "# vtk DataFile Version 3.0" << endl;
  vtk << "#Simulation Euler" << endl;
  vtk << "ASCII" << endl;
  vtk<<"\n";
  vtk << "DATASET UNSTRUCTURED_GRID" << endl;
  vtk << "POINTS " << (Nx+1)*(Ny+1) *(Nz+1)<< " DOUBLE" << endl;
    
  for(int i=0; i<Nx+1; i++){
    for(int j=0; j<Ny+1; j++){ 
      for(int k=0; k<Nz+1; k++){ 
	vtk << i*dx << " " << j*dy << " " <<k*dz << " "<< endl;
      }
    }
  }
  vtk<<"\n";

  //Computation of the number of true fluid cells
  int Nfluides = 0;
  for(int i=marge; i<Nx+marge; i++){
    for(int j=marge; j<Ny+marge; j++){ 
      for(int k=marge; k<Nz+marge; k++){
	Cellule& c = grille[i][j][k]; 
	if(abs(c.alpha-1.)>eps){
	  Nfluides++;
	}
      }
    }
  }
	
  vtk << "CELLS " << Nfluides << " " << 9*Nfluides<< endl;
    
  for(int i=marge; i<Nx+marge; i++){
    for(int j=marge; j<Ny+marge; j++){ 
      for(int k=marge; k<Nz+marge; k++){
	Cellule& c = grille[i][j][k]; 
	if(abs(c.alpha-1.)>eps){
	  vtk << 8 << " " << (k-marge)+(j-marge)*(Nz+1)+(i-marge)*(Nz+1)*(Ny+1) << " " << ((k-marge)+1)+(j-marge)*(Nz+1)+(i-marge)*(Nz+1)*(Ny+1) << " " << ((k-marge)+1)+((j-marge)+1)*(Nz+1) + (i-marge)*(Nz+1)*(Ny+1) << " "<< (k-marge)+((j-marge)+1)*(Nz+1) + (i-marge)*(Nz+1)*(Ny+1)<< " " <<  (k-marge)+(j-marge)*(Nz+1)+((i-marge)+1)*(Nz+1)*(Ny+1)<< " "<< ((k-marge)+1)+(j-marge)*(Nz+1) + ((i-marge)+1)*(Nz+1)*(Ny+1) << " " <<  ((k-marge)+1)+((j-marge)+1)*(Nz+1)+((i-marge)+1)*(Nz+1)*(Ny+1)<< " " << (k-marge)+((j-marge)+1)*(Nz+1)+((i-marge)+1)*(Nz+1)*(Ny+1)<< endl;
	}
      }
    }
  }
  vtk<<"\n";
  vtk << "CELL_TYPES " <<Nfluides<<endl;
  for(int k=0; k<Nfluides; k++){ 
    vtk<<12<<endl;
  }
    
  vtk<<"\n";
  vtk << "CELL_DATA " << Nfluides << endl;
  //Pressure
  vtk << "SCALARS pressure double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(int i=marge; i<Nx+marge; i++){
    for(int j=marge; j<Ny+marge; j++){ 
      for(int k=marge; k<Nz+marge; k++){ 
	Cellule& c = grille[i][j][k]; 
	if(abs(c.alpha-1.)>eps){ 
	  vtk << grille[i][j][k].p << endl;
	} 				
      }
    }
  }
    
  vtk<<"\n";
  //Density
  vtk << "SCALARS density double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(int i=marge; i<Nx+marge; i++){
    for(int j=marge; j<Ny+marge; j++){ 
      for(int k=marge; k<Nz+marge; k++){ 
	Cellule& c = grille[i][j][k]; 
	if(abs(c.alpha-1.)>eps){ 
	  vtk << grille[i][j][k].rho << endl;
	}
      }
    }
  }
    
  vtk<<"\n";
  //Velocity x-component u
  vtk << "SCALARS u double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(int i=marge; i<Nx+marge; i++){
    for(int j=marge; j<Ny+marge; j++){ 
      for(int k=marge; k<Nz+marge; k++){ 
	Cellule& c = grille[i][j][k]; 
	if(abs(c.alpha-1.)>eps){ 
	  vtk << grille[i][j][k].u << endl;
	} 
      }
    }
  }
  vtk<<"\n";
  //Velocity y-component v
  vtk << "SCALARS v double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(int i=marge; i<Nx+marge; i++){
    for(int j=marge; j<Ny+marge; j++){ 
      for(int k=marge; k<Nz+marge; k++){ 
	Cellule& c = grille[i][j][k]; 
	if(abs(c.alpha-1.)>eps){ 
	  vtk << grille[i][j][k].v << endl;
	}       
      }
    }
  }
  vtk<<"\n";
  //Velocity z-component w
  vtk << "SCALARS w double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(int i=marge; i<Nx+marge; i++){
    for(int j=marge; j<Ny+marge; j++){ 
      for(int k=marge; k<Nz+marge; k++){ 
	Cellule& c = grille[i][j][k]; 
	if(abs(c.alpha-1.)>eps){ 
	  vtk << grille[i][j][k].w << endl;
	}       
      }
    }
  }
}

Cellule Grille::voisin_fluide(const Cellule &c, bool &target){
  double dir = 0.; 
  int i= c.i; 
  int j= c.j;
  int k= c.k;
	
  dir = std::min(std::min(std::min(std::min(std::min(c.kappai, c.kappaj), c.kappak), grille[i-1][j][k].kappai), grille[i][j-1][k].kappaj), grille[i][j][k-1].kappak);
	
	
  if ( ((i+1)>=marge) && ((i+1)< (Nx+marge)) &&  (std::abs(dir - c.kappai)<eps) && (grille[i+1][j][k].alpha <eps) && (grille[i+1][j][k].p > 0.) && (grille[i+1][j][k].rho > 0.) && (!grille[i+1][j][k].vide))
  {
    return grille[i+1][j][k];
  }
  else if ( ((j+1)>=marge) && ( (j+1)< (Ny+marge)) &&  (std::abs(dir - c.kappaj)<eps) && (grille[i][j+1][k].alpha  <eps) && (grille[i][j+1][k].p > 0.) && (grille[i][j+1][k].rho > 0.) && (!grille[i][j+1][k].vide))
  {
    return grille[i][j+1][k];
  }
  else if( ((k+1)>=marge) && ((k+1)< (Nz+marge)) &&  (std::abs(dir - c.kappak)<eps) && (grille[i][j][k+1].alpha  <eps) && (grille[i][j][k+1].p > 0.) && (grille[i][j][k+1].rho > 0.) && (!grille[i][j][k+1].vide))
  {
    return grille[i][j][k+1];
  }
  else if ( ((i-1)>=marge) && ((i-1)< (Nx+marge)) &&  (std::abs(dir - grille[i-1][j][k].kappai)<eps) && (grille[i-1][j][k].alpha< eps) && (grille[i-1][j][k].p > 0.) && (grille[i-1][j][k].rho > 0.) && (!grille[i-1][j][k].vide) )
  {
    return grille[i-1][j][k];
  }
  else if ( ((j-1)>=marge) && ((j-1)< (Ny+marge)) &&  (std::abs(dir - grille[i][j-1][k].kappaj)<eps) && (grille[i][j-1][k].alpha <eps) && (grille[i][j-1][k].p > 0.) && (grille[i][j-1][k].rho > 0.) && (!grille[i][j-1][k].vide))
  {
    return grille[i][j-1][k];
  }
  else if ( ((k-1)>=marge) && ((k-1)< (Nz+marge)) && (std::abs(dir - grille[i][j][k-1].kappak)<eps) && (grille[i][j][k-1].alpha  <eps) && (grille[i][j][k-1].p > 0.) && (grille[i][j][k-1].rho > 0.) && (!grille[i][j][k-1].vide))
  {
    return grille[i][j][k-1];
  }
  else{
		
    target=false;
    return c;
  }
}

Cellule Grille::voisin_mixt(const Cellule &c, bool &target){
	
  double dir = 0.; 
  int i= c.i; 
  int j= c.j;
  int k= c.k;
  dir = std::min(std::min(std::min(std::min(std::min(c.kappai, c.kappaj), c.kappak), grille[i-1][j][k].kappai), grille[i][j-1][k].kappaj), grille[i][j][k-1].kappak);

	
  if ( ((i+1)>=marge) && ((i+1)< (Nx+marge)) && (std::abs(dir - c.kappai)<eps) && ( (c.alpha > grille[i+1][j][k].alpha) || (std::abs(c.alpha - grille[i+1][j][k].alpha)<eps) ) 
       && (grille[i+1][j][k].p > 0.) && (grille[i+1][j][k].rho > 0.) && (!grille[i+1][j][k].vide))
  {
    return grille[i+1][j][k];
  }
  else if ( ((j+1)>=marge) && ( (j+1)< (Ny+marge)) && (std::abs(dir - c.kappaj)<eps) && ((c.alpha > grille[i][j+1][k].alpha) || (std::abs(c.alpha - grille[i][j+1][k].alpha)<eps))                     && (grille[i][j+1][k].p > 0.) && (grille[i][j+1][k].rho > 0.)  && (!grille[i][j+1][k].vide))
  {
    return grille[i][j+1][k];
  }
  else if( ((k+1)>=marge) && ((k+1)< (Nz+marge)) && (std::abs(dir - c.kappak)<eps) && ((c.alpha > grille[i][j][k+1].alpha) || (std::abs(c.alpha - grille[i][j][k+1].alpha)<eps))
	   && (grille[i][j][k+1].p > 0.) && (grille[i][j][k+1].rho > 0.) && (!grille[i][j][k+1].vide))
  {
    return grille[i][j][k+1];
  }
  else if ( ((i-1)>=marge) && ((i-1)< (Nx+marge)) && (std::abs(dir - grille[i-1][j][k].kappai)<eps) && ((c.alpha > grille[i-1][j][k].alpha) || (std::abs(c.alpha - grille[i-1][j][k].alpha)<eps)) && (grille[i-1][j][k].p > 0.) && (grille[i-1][j][k].rho > 0.) && (!grille[i-1][j][k].vide))
  {
    return grille[i-1][j][k];
  }
  else if ( ((j-1)>=marge) && ((j-1)< (Ny+marge)) && (std::abs(dir - grille[i][j-1][k].kappaj)<eps) && ((c.alpha > grille[i][j-1][k].alpha) || (std::abs(c.alpha - grille[i][j-1][k].alpha)<eps)) && (grille[i][j-1][k].p > 0.) && (grille[i][j-1][k].rho > 0.) && (!grille[i][j-1][k].vide))
  {
    return grille[i][j-1][k];
  }
  else if ( ((k-1)>=marge) && ((k-1)< (Nz+marge)) && (std::abs(dir - grille[i][j][k-1].kappak)<eps) && ((c.alpha > grille[i][j][k-1].alpha) || (std::abs(c.alpha - grille[i][j][k-1].alpha)<eps)) && (grille[i][j][k-1].p > 0.) && (grille[i][j][k-1].rho > 0.) && (!grille[i][j][k-1].vide) )
  {
    return grille[i][j][k-1];
  }
  else if ( ((i+1)>=marge) && ((i+1)< (Nx+marge)) && (std::abs(dir - c.kappai)<eps) || ( (c.alpha > grille[i+1][j][k].alpha) || (std::abs(c.alpha - grille[i+1][j][k].alpha)<eps) ) 
	    && (grille[i+1][j][k].p > 0.) && (grille[i+1][j][k].rho > 0.) && (!grille[i+1][j][k].vide))
  {
    return grille[i+1][j][k];
  }
  else if ( ((j+1)>=marge) && ( (j+1)< (Ny+marge)) && (std::abs(dir - c.kappaj)<eps) || ((c.alpha > grille[i][j+1][k].alpha) || (std::abs(c.alpha - grille[i][j+1][k].alpha)<eps))                     && (grille[i][j+1][k].p > 0.) && (grille[i][j+1][k].rho > 0.) && (!grille[i][j+1][k].vide))
  {
    return grille[i][j+1][k];
  }
  else if( ((k+1)>=marge) && ((k+1)< (Nz+marge)) && (std::abs(dir - c.kappak)<eps) || ((c.alpha > grille[i][j][k+1].alpha) || (std::abs(c.alpha - grille[i][j][k+1].alpha)<eps))
	   && (grille[i][j][k+1].p > 0.) && (grille[i][j][k+1].rho > 0.) && (!grille[i][j][k+1].vide))
  {
    return grille[i][j][k+1];
  }
  else if ( ((i-1)>=marge) && ((i-1)< (Nx+marge)) && (std::abs(dir - grille[i-1][j][k].kappai)<eps) || ((c.alpha > grille[i-1][j][k].alpha) || (std::abs(c.alpha - grille[i-1][j][k].alpha)<eps)) && (grille[i-1][j][k].p > 0.) && (grille[i-1][j][k].rho > 0.) && (!grille[i-1][j][k].vide))
  {
    return grille[i-1][j][k];
  }
  else if ( ((j-1)>=marge) && ((j-1)< (Ny+marge)) && (std::abs(dir - grille[i][j-1][k].kappaj)<eps) || ((c.alpha >= grille[i][j-1][k].alpha) || (std::abs(c.alpha - grille[i][j-1][k].alpha)<eps)) && (grille[i][j-1][k].p > 0.) && (grille[i][j-1][k].rho > 0.) && (!grille[i][j-1][k].vide))
  {
    return grille[i][j-1][k];
  }
  else if ( ((k-1)>=marge) && ((k-1)< (Nz+marge)) && (std::abs(dir - grille[i][j][k-1].kappak)<eps) || ((c.alpha > grille[i][j][k-1].alpha) || (std::abs(c.alpha - grille[i][j][k-1].alpha)<eps)) && (grille[i][j][k-1].p > 0.) && (grille[i][j][k-1].rho > 0.) && (!grille[i][j][k-1].vide) )
  {
    return grille[i][j][k-1];
  }
  else{
		
    target=false;
    return  voisin(c);
  }
}

Cellule Grille::voisin(const Cellule &c) {
	
  int i= c.i; 
  int j= c.j;
  int k= c.k;
  double dir = i; 
	
  dir = std::min(std::min(std::min(std::min(std::min(c.kappai, c.kappaj), c.kappak), grille[i-1][j][k].kappai), grille[i][j-1][k].kappaj), grille[i][j][k-1].kappak);
	
  if (((i+1)>=marge) && ((i+1)< (Nx+marge)) && (std::abs(dir - c.kappai)<eps)  && (grille[i+1][j][k].p > 0.) && (grille[i+1][j][k].rho > 0.) && (!grille[i+1][j][k].vide))
  {
    return grille[i+1][j][k];
  }
  else if ( ((j+1)>=marge) && ( (j+1)< (Ny+marge)) && (std::abs(dir - c.kappaj)<eps) && (grille[i][j+1][k].p > 0.) && (grille[i][j+1][k].rho > 0.) && (!grille[i][j+1][k].vide))
  {
    return grille[i][j+1][k];
  }
  else if ( ((k+1)>=marge) && ((k+1)< (Nz+marge)) && (std::abs(dir - c.kappak)<eps) && (grille[i][j][k+1].p > 0.) && (grille[i][j][k+1].rho > 0.) && (!grille[i][j][k+1].vide))
  {
    return grille[i][j][k+1];
  }
  else if (((i-1)>=marge) && ((i-1)< (Nx+marge)) && (std::abs(dir - grille[i-1][j][k].kappai)) && (grille[i-1][j][k].p > 0.) && (grille[i-1][j][k].rho > 0.) && (!grille[i-1][j][k].vide))
  {
    return grille[i-1][j][k];
  }
  else if (((j-1)>=marge) && ( (j-1)< (Ny+marge)) && (std::abs(dir - grille[i][j-1][k].kappaj)) && (grille[i][j-1][k].p > 0.) && (grille[i][j-1][k].rho > 0.) && (!grille[i][j-1][k].vide))
  {
    return grille[i][j-1][k];
  }
  else if( ((k-1)>=marge) && ((k-1)< (Nz+marge)) && (std::abs(dir - grille[i][j][k-1].kappak)) && (grille[i][j][k-1].p > 0.) && (grille[i][j][k-1].rho > 0.) && (!grille[i][j][k-1].vide))
  {
    return grille[i][j][k-1];
  }
  //else if( ((i+1)>=marge) && ((i+1)< (Nx+marge)) && (std::abs(dir - c.kappai)<eps) && (!grille[i+1][j][k].vide))
  else if( ((i+1)>=marge) && ((i+1)< (Nx+marge)) && (grille[i+1][j][k].p>0.) && (grille[i+1][j][k].rho>0.) && (!grille[i+1][j][k].vide))
  {
    return grille[i+1][j][k];
  }
  else if ( ((j+1)>=marge) && ((j+1)< (Ny+marge)) && (grille[i][j+1][k].p>0.) && (grille[i][j+1][k].rho>0.) && (!grille[i][j+1][k].vide))
  {
    return grille[i][j+1][k];
  }
  else if (((k+1)>=marge) && ((k+1)< (Nz+marge)) && (grille[i][j][k+1].p>0.) && (grille[i][j][k+1].rho>0.) && (!grille[i][j][k+1].vide))
  {
    return grille[i][j][k+1];
  }
  else if ( ((i-1)>=marge) && ((i-1)< (Nx+marge)) && (grille[i-1][j][k].p>0.) && (grille[i-1][j][k].rho>0.) && (!grille[i-1][j][k].vide))
  {
    return grille[i-1][j][k];
  }
  else if(((j-1)>=marge) && ((j-1)< (Ny+marge)) && (grille[i][j-1][k].p>0.) && (grille[i][j-1][k].rho>0.) && (!grille[i][j-1][k].vide))
  {
    return grille[i][j-1][k];
  }
  else if (((k-1)>=marge) && ((k-1)< (Nz+marge)) && (grille[i][j][k-1].p>0.) && (grille[i][j][k-1].rho>0.) && (!grille[i][j][k-1].vide))
  {
    return grille[i][j][k-1];
  }
  else {
    return grille[i][j][k];
  }
}
Cellule Grille::cible(const Cellule &c, std::vector< std::vector<int> > & tab_cible ){
	
  bool target = true;
  Cellule cell_cible;
  cell_cible = voisin_fluide(c, target);
  if(target){ 
    return cell_cible;
  }
  else{
    target = true;
    cell_cible= voisin_mixt(c,target);
    int l=cell_cible.i, m=cell_cible.j, n=cell_cible.k;
    bool cycle=false;
    for(int iter=0; iter<tab_cible.size(); iter++ ){
      if(tab_cible[iter][0]==l && tab_cible[iter][1]==m && tab_cible[iter][2]==n){
	cycle =true;
      }
    }
    if(cycle){
      int iter= tab_cible.size()-1;
      return grille[tab_cible[iter][0]][tab_cible[iter][1]][tab_cible[iter][2]];
    }
    else{
      std::vector<int> poz(3); poz[0]= l; poz[1] = m; poz[2] = n; tab_cible.push_back(poz);
      if (target && l>=marge && l< Nx+marge && m>=marge && m< Ny+marge && n>=marge && n< Nz+marge) 
      {	
	return cible(cell_cible, tab_cible); 
      }
      else{
	return cell_cible;
      }  
    }
  }
}

#endif
