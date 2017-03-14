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
  \authors Maria Adela Puscas and Laurent Monasse
  \brief Definition of specific coupling functions.
  \details Computation of fluid forces and torques acting on the solid, modifications of the fluid fluxes, filling of ghost cells, computation of the swept quantity.
  \warning  <b> Specific coupling procedures ! </b>
*/


#include "fluide.hpp"
#include "intersections.cpp"

/*!\brief Computation of fluid forces (\a Particule.Ff) and torques (\a Particule.Mf) applied on the \a Solide.
  \details Let \a f an interface element, the pressure force exerted by the fluid on interface \a f is given by:
 \f{eqnarray*}{
 F_f  =  (- p^x \, A_f n^{x}_f, \,- p^y \,A_f n^{y}_f, \,- p^z \,A_f n^{z}_f )^t
 \f} \n
 where  \f$ A_f \f$ is the area of interface f,  \f$ n_f \f$  the exterior normal to f, and \f$ p^x, p^y, p^z \f$ are the effective pressures in the x, y and z directions during the time-step (\a Cellule.pdtx, \a  Cellule.pdty and \a Cellule.pdtz).
 The fluid torque exerted on f is given by:
 \f{eqnarray*}{
 M_f = F_f  \wedge (X_f - X_I),
 \f}
 where \f$ X_f \f$ is the center of interface f and \f$  X_I \f$ is the center of the particle containing f.
 These forces are transmitted to the solid as being constant during the time-step.
  \param S Solid
  \param dt Time-step
  \warning <b> Specific coupling procedure ! </b>
  \return void
 */
void Grille::Forces_fluide(Solide& S, const double dt){
	
  Vector_3 Ffluide(0.,0.,0.);
  //Update of fluid forces and torques on the solid
  for(int iter_s=0; iter_s<S.size(); iter_s++){ 
		
    S.solide[iter_s].Ffprev = S.solide[iter_s].Ff;
    S.solide[iter_s].Mfprev = S.solide[iter_s].Mf;
    Point_3 Xn = S.solide[iter_s].x0 + S.solide[iter_s].Dx;
    double fx=0.; double fy=0.; double fz=0.;
    Kernel::FT mx = 0.,my = 0. ,mz = 0.;
		
    for(int it=0; it<S.solide[iter_s].triangles.size(); it++){
      for(int iter=0; iter<S.solide[iter_s].Position_Triangles_interface[it].size(); iter++)
      {  
	double aire= std::sqrt(CGAL::to_double(S.solide[iter_s].Triangles_interface[it][iter].squared_area()));
	if(dt>eps){	
	  int i= S.solide[iter_s].Position_Triangles_interface[it][iter][0]; 
	  int j= S.solide[iter_s].Position_Triangles_interface[it][iter][1]; 
	  int k= S.solide[iter_s].Position_Triangles_interface[it][iter][2]; 
	  if(i>=marge && i<Nx+marge && j>=marge && j<Ny+marge && k>=marge && k<Nz+marge){
	    double tempx = (grille[i][j][k].pdtx/dt) * aire * (CGAL::to_double(S.solide[iter_s].normales[it].x()));
	    double tempy = (grille[i][j][k].pdty/dt) * aire * (CGAL::to_double(S.solide[iter_s].normales[it].y()));
	    double tempz = (grille[i][j][k].pdtz/dt) * aire * (CGAL::to_double(S.solide[iter_s].normales[it].z()));
	    Vector_3 temp_Mf = cross_product(Vector_3(Xn,Point_3(centroid(S.solide[iter_s].Triangles_interface[it][iter].operator[](0),
									  S.solide[iter_s].Triangles_interface[it][iter].operator[](1),
									  S.solide[iter_s].Triangles_interface[it][iter].operator[](2)))), 
					     Vector_3(-tempx,-tempy,-tempz));
	    fx-= tempx; fy-= tempy; fz-= tempz;
	    mx+= temp_Mf.x(); my+= temp_Mf.y(); mz+= temp_Mf.z();
	  }
	}
      }
    }
    S.solide[iter_s].Ff = Vector_3(fx,fy,fz);
    S.solide[iter_s].Mf = Vector_3(CGAL::to_double(mx),CGAL::to_double(my),CGAL::to_double(mz)); 
    Ffluide = Ffluide + S.solide[iter_s].Ff;
  }
  cout<<"Fluid forces "<<Ffluide<<endl;
  for(int it=0; it<S.solide.size(); it++){
    for(int i=0; i<S.solide[it].faces.size(); i++){
      if(S.solide[it].faces[i].voisin == -2){
	S.solide[it].faces[i].voisin = -1;
      }
    }
  }
}	

/*!\brief Modification of fluid fluxes and discrete balance on a cell (cut-cell method).
   \details The fluid sees the presence of the solid through this function. The final value of the state \f$ U^{n+1}_{i, j, k}\f$ in the cell is computed using:
 \f{eqnarray*}{
 \left( 1-  \Lambda_{i,j,k}^{n+1} \right) U^{n+1}_{i,j,k}   = \left( 1-  \Lambda_{i,j,k} ^{n+1}\right) U^n_{i,j,k}  + \Delta t \, \left(   \frac{(1-\lambda_{i-1/2,j,k}^{n+1} )}{\Delta x_{ i,j,k}} F_{i-1/2, j, k}^{n+1/2} -\frac{(1-\lambda_{i+1/2,j,k}^{n+1} )}{\Delta x_{ i,j,k}} F_{i+1/2, j, k}^{n+1/2}  + ...\right)	\f}
 \f{eqnarray*}{+  \frac{\Delta t}{V_{i,j,k}} \sum_{f \in C_{i,j,k}}{A_{f}} {\phi}_{f_{ i,j,k}}  +   \sum_{f \in C_{i,j,k}} \Delta U^{n, n+1}_{f_{ i,j,k}}  .
 \f}
 where \f$ \Lambda_{i,j,k}^{n+1} \f$: \a Cellule.alpha (solid occupancy ratio in cell (i,j,k)),\n
 \f$ \lambda_{i+1/2,j,k}^{n+1} \f$: \a Cellule.kappai (solid occupancy ratio on the faces of cell (i,j,k)),\n
 \f$ F_{i+1/2, j, k}^{n+1/2} \f$:   \a  Cellule.fluxi (flux to the right of cell (i,j,k)), \n
 \f$ V_{i,j,k} \f$: volume of cell (i,j,k), \n
 \f{eqnarray*}{ \phi_{f_{ i,j,k}} = (0, \Pi_x, \Pi_y,\Pi_z, \Pi_v)^t. 	\f}
 \f$ \Pi_x, \Pi_y,\Pi_z, \Pi_v\f$: \a Cellule.phi_x, \a Cellule.phi_y, Cellule.phi_z, \a Cellule.phi_v (boundary fluxes),\n
 \f$ \Delta U^{n, n+1}_{f_{ i,j,k}} 	\f$:  \a Cellule.delta_w (swept quantity).
 	\param dt Time-step
 	\warning <b> Specific coupling procedure ! </b>
 	\return void
 */
void Grille::Modif_fnum(const double dt){
	
  double phi_x=0., phi_y=0., phi_z=0.;
  double vol=deltax*deltay*deltaz;
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){ 
      for(int k=marge;k<Nz+marge;k++){
	Cellule& c = grille[i][j][k];
	if(std::abs(c.alpha-1.)>eps && !c.vide){
	  Cellule& ci = grille[i-1][j][k];   
	  Cellule& cj = grille[i][j-1][k];   
	  Cellule& ck = grille[i][j][k-1];   
           
	  c.flux_modif[0] = 0.;
	  c.flux_modif[1] = c.phi_x;
	  c.flux_modif[2] = c.phi_y;
	  c.flux_modif[3] = c.phi_z;
	  c.flux_modif[4] = c.phi_v;
	  for(int l=0.; l<5; l++){  
	    c.flux_modif[l] -= (1.-c.kappai)*c.dtfxi[l] - (1.-ci.kappai)*ci.dtfxi[l]
	      + (1.-c.kappaj)*c.dtfyj[l] - (1.-cj.kappaj)*cj.dtfyj[l]
	      + (1.-c.kappak)*c.dtfzk[l] - (1.-ck.kappak)*ck.dtfzk[l] - c.delta_w[l];
	    c.flux_modif[l] /= (1.-c.alpha);
	  }
	  //Update of the cell state
	  c.rho = c.rho0  +  c.flux_modif[0];
	  c.impx = c.impx0 + c.flux_modif[1];
	  c.impy = c.impy0 + c.flux_modif[2];
	  c.impz = c.impz0 + c.flux_modif[3];
	  c.rhoE = c.rhoE0 + c.flux_modif[4];
	  if(std::abs(c.rho) > eps_vide){
	    c.u = c.impx/c.rho;
	    c.v = c.impy/c.rho;
	    c.w = c.impz/c.rho;
	    c.p = (gam-1.)*(c.rhoE-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2. - c.rho*c.w*c.w/2.);
	    if(std::abs(c.p) > eps_vide){
	      c.vide = false;
	    }
	    phi_x+=c.phi_x*vol/dt; phi_y+=c.phi_y*vol/dt; phi_z+=c.phi_z*vol/dt;
	  }
	  if( (abs(c.rho) <= eps_vide) || (abs(c.p) <= eps_vide)){
	    c.u = 0.; c.v = 0.; c.w = 0.; c.p = 0.;
	    c.impx=0.; c.impy=0.; c.impz=0.; c.rhoE=0.;
	    c.vide = true;
	  }
	  else{c.vide = false;}
	}      
      }
    }
  }
  Vector_3 Phi(phi_x, phi_y, phi_z); cout<<" Boundary flux "<<Phi<<endl;
}

/*!\brief Conservative mixing of small cut-cells.
   \details We define a small cut-cell as a cell such that \f$ alpha > epsa \f$ (\a Cellule.alpha solid occupancy ratio, and \a epsa: limit ratio for small cut-cells defined in parametres.hpp). In order not to modify the time-step which ensuring the CFL condition, the small cut-cells are merged with their neighbours. Denote \a p a small cut-cell and \a g a neighbouring cell such that \a g is totally fluid~(\f$ alpha_g = 0 \f$ ). Define the following exchange terms:
 \f{eqnarray*}{ E_{pg} = \frac{1}{2 - alpha_p} (U_g - U_{p}), \quad  E_{gp} = \frac{1-  alpha_p}{2 - alpha_p} (U_p - U_{g}) \f}
 and set:
 \f{eqnarray*}{
 U_p = U_{p} + E_{pg}, \quad  \quad U_g = U_{g} + E_{gp} \f}
 	\warning <b> Specific coupling procedure ! </b>
 	\return void
 */
void Grille::Mixage(){
	
	
  bool test_fini = true;
	
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){ 
      for(int k=marge;k<Nz+marge;k++){
	Cellule& cp = grille[i][j][k];
	bool test=true;
	if( (cp.alpha>epsa ||cp.p <0. || cp.rho<0.) && abs(cp.alpha-1.)>eps && !cp.vide){
					
	  for(int ii=-1; ii<=1 && test; ii++){
	    for(int jj=-1; jj<=1 && test; jj++){
	      for(int kk=-1; kk<=1 && test; kk++){
		if (grille[i+ii][j+jj][k+kk].alpha <eps && grille[i+ii][j+jj][k+kk].p>0. && grille[i+ii][j+jj][k+kk].rho>0. && i+ii>=marge && i+ii<Nx+marge && j+jj>=marge && j+jj<Ny+marge && k+kk>=marge && k+kk<Nz+marge && !grille[i+ii][j+jj][k+kk].vide)
		{
		  test=false;
		  Cellule& cg = grille[i+ii][j+jj][k+kk];
		  double temp_rhop= cp.rho;
		  double temp_rhog=cg.rho;
		  cp.Mrho =  (cg.rho - cp.rho)/(2. - cp.alpha) ;
		  cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
		  cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
		  cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
		  cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
									
		  cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
		  cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
		  cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
		  cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
		  cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
									
									
		  cp.rho += cp.Mrho;
		  cp.impx += cp.Mimpx;
		  cp.impy += cp.Mimpy;
		  cp.impz += cp.Mimpz;
		  cp.rhoE += cp.MrhoE;
		  cp.u = cp.impx/cp.rho;
		  cp.v = cp.impy/cp.rho;
		  cp.w = cp.impz/cp.rho;
		  cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
									
		  cg.rho += cg.Mrho;
		  cg.impx += cg.Mimpx;
		  cg.impy += cg.Mimpy;
		  cg.impz += cg.Mimpz;
		  cg.rhoE += cg.MrhoE;
		  cg.u = cg.impx/cg.rho;
		  cg.v = cg.impy/cg.rho;
		  cg.w = cg.impz/cg.rho;
		  cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
		}						
	      }
	    }
	  }
					
	  if(test){
						
	    if (grille[i-2][j][k].alpha == 0. && grille[i-2][j][k].p>0. && grille[i-2][j][k].rho>0. && i-2>=marge && !grille[i-2][j][k].vide)
	    {
	      Cellule& cg = grille[i-2][j][k];
							
	      cp.Mrho =  (cg.rho - cp.rho)/(2. - cp.alpha) ;
	      cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
	      cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
	      cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
	      cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
							
	      cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
	      cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
	      cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
	      cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
	      cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
							
							
	      cp.rho += cp.Mrho;
	      cp.impx += cp.Mimpx;
	      cp.impy += cp.Mimpy;
	      cp.impz += cp.Mimpz;
	      cp.rhoE += cp.MrhoE;
	      cp.u = cp.impx/cp.rho;
	      cp.v = cp.impy/cp.rho;
	      cp.w = cp.impz/cp.rho;
	      cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
							
	      cg.rho += cg.Mrho;
	      cg.impx += cg.Mimpx;
	      cg.impy += cg.Mimpy;
	      cg.impz += cg.Mimpz;
	      cg.rhoE += cg.MrhoE;
	      cg.u = cg.impx/cg.rho;
	      cg.v = cg.impy/cg.rho;
	      cg.w = cg.impz/cg.rho;
	      cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
		
	      test = false;
	    }
	    else if (grille[i+2][j][k].alpha == 0. && grille[i+2][j][k].p>0. && grille[i+2][j][k].rho>0. &&  i+2<Nx+marge && !grille[i+2][j][k].vide)
	    {
	      Cellule& cg = grille[i+2][j][k];
							
	      cp.Mrho =  (cg.rho - cp.rho)/(2. - cp.alpha) ;
	      cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
	      cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
	      cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
	      cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
							
	      cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
	      cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
	      cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
	      cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
	      cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
							
	      cp.rho += cp.Mrho;
	      cp.impx += cp.Mimpx;
	      cp.impy += cp.Mimpy;
	      cp.impz += cp.Mimpz;
	      cp.rhoE += cp.MrhoE;
	      cp.u = cp.impx/cp.rho;
	      cp.v = cp.impy/cp.rho;
	      cp.w = cp.impz/cp.rho;
	      cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
							
	      cg.rho += cg.Mrho;
	      cg.impx += cg.Mimpx;
	      cg.impy += cg.Mimpy;
	      cg.impz += cg.Mimpz;
	      cg.rhoE += cg.MrhoE;
	      cg.u = cg.impx/cg.rho;
	      cg.v = cg.impy/cg.rho;
	      cg.w = cg.impz/cg.rho;
	      cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
		
	      test = false;
	    }
						
	    else if (grille[i][j-2][k].alpha == 0. && grille[i][j-2][k].p>0. && grille[i][j-2][k].rho>0. && j-2>=marge && !grille[i][j-2][k].vide)
	    {
	      Cellule& cg = grille[i][j-2][k];
							
	      cp.Mrho =  (cg.rho - cp.rho)/(2. - cp.alpha) ;
	      cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
	      cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
	      cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
	      cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
							
	      cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
	      cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
	      cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
	      cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
	      cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
							
							
	      cp.rho += cp.Mrho;
	      cp.impx += cp.Mimpx;
	      cp.impy += cp.Mimpy;
	      cp.impz += cp.Mimpz;
	      cp.rhoE += cp.MrhoE;
	      cp.u = cp.impx/cp.rho;
	      cp.v = cp.impy/cp.rho;
	      cp.w = cp.impz/cp.rho;
	      cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
							
	      cg.rho += cg.Mrho;
	      cg.impx += cg.Mimpx;
	      cg.impy += cg.Mimpy;
	      cg.impz += cg.Mimpz;
	      cg.rhoE += cg.MrhoE;
	      cg.u = cg.impx/cg.rho;
	      cg.v = cg.impy/cg.rho;
	      cg.w = cg.impz/cg.rho;
	      cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
		
	      test = false;
							
	    }
	    else if (grille[i][j+2][k].alpha == 0. && grille[i][j+2][k].p>0. && grille[i][j+2][k].rho>0.&& j+2<Ny+marge && !grille[i][j+2][k].vide)
	    {
	      Cellule& cg = grille[i][j+2][k];
							
	      cp.Mrho =  (cg.rho - cp.rho)/(2. - cp.alpha) ;
	      cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
	      cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
	      cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
	      cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
							
	      cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
	      cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
	      cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
	      cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
	      cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
							
	      cp.rho += cp.Mrho;
	      cp.impx += cp.Mimpx;
	      cp.impy += cp.Mimpy;
	      cp.impz += cp.Mimpz;
	      cp.rhoE += cp.MrhoE;
	      cp.u = cp.impx/cp.rho;
	      cp.v = cp.impy/cp.rho;
	      cp.w = cp.impz/cp.rho;
	      cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
							
	      cg.rho += cg.Mrho;
	      cg.impx += cg.Mimpx;
	      cg.impy += cg.Mimpy;
	      cg.impz += cg.Mimpz;
	      cg.rhoE += cg.MrhoE;
	      cg.u = cg.impx/cg.rho;
	      cg.v = cg.impy/cg.rho;
	      cg.w = cg.impz/cg.rho;
	      cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
		
	      test = false;
	    }
	    else if (grille[i][j][k-2].alpha == 0. && grille[i][j][k-2].p>0. && grille[i][j][k-2].rho>0.&& k-2>=marge && !grille[i][j][k-2].vide)
	    {
	      Cellule& cg = grille[i][j][k-2];
							
	      cp.Mrho =  (cg.rho - cp.rho)/(2. - cp.alpha) ;
	      cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
	      cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
	      cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
	      cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
							
	      cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
	      cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
	      cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
	      cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
	      cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
							
							
	      cp.rho += cp.Mrho;
	      cp.impx += cp.Mimpx;
	      cp.impy += cp.Mimpy;
	      cp.impz += cp.Mimpz;
	      cp.rhoE += cp.MrhoE;
	      cp.u = cp.impx/cp.rho;
	      cp.v = cp.impy/cp.rho;
	      cp.w = cp.impz/cp.rho;
	      cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
							
	      cg.rho += cg.Mrho;
	      cg.impx += cg.Mimpx;
	      cg.impy += cg.Mimpy;
	      cg.impz += cg.Mimpz;
	      cg.rhoE += cg.MrhoE;
	      cg.u = cg.impx/cg.rho;
	      cg.v = cg.impy/cg.rho;
	      cg.w = cg.impz/cg.rho;
	      cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
		
	      test = false;
	    }
	    else if(grille[i][j][k+2].alpha == 0. && grille[i][j][k+2].p>0. && grille[i][j][k+2].rho>0. && k+2 < Nz+marge && !grille[i][j][k+2].vide)
	    {
	      Cellule& cg = grille[i][j][k+2];
							
	      cp.Mrho =  (cg.rho - cp.rho)/(2. - cp.alpha) ;
	      cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
	      cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
	      cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
	      cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
							
	      cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
	      cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
	      cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
	      cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
	      cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
							
							
	      cp.rho += cp.Mrho;
	      cp.impx += cp.Mimpx;
	      cp.impy += cp.Mimpy;
	      cp.impz += cp.Mimpz;
	      cp.rhoE += cp.MrhoE;
	      cp.u = cp.impx/cp.rho;
	      cp.v = cp.impy/cp.rho;
	      cp.w = cp.impz/cp.rho;
	      cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
							
	      cg.rho += cg.Mrho;
	      cg.impx += cg.Mimpx;
	      cg.impy += cg.Mimpy;
	      cg.impz += cg.Mimpz;
	      cg.rhoE += cg.MrhoE;
	      cg.u = cg.impx/cg.rho;
	      cg.v = cg.impy/cg.rho;
	      cg.w = cg.impz/cg.rho;
	      cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
		
	      test = false;
	    }
						
	  }
	  else if(test){
	    std::cout<<"No mixing target cell"<<std::endl; 
	    std::cout<< "Position of the cell center: "<<grille[i][j][k].x << " "<<grille[i][j][k].y << " "<<grille[i][j][k].z << " "<< " rho "<<grille[i][j][k].rho  << " p "<<grille[i][j][k].p <<" alpha " << grille[i][j][k].alpha<<std::endl;
	    std::cout<<"Neighbouring cells: "<<std::endl;
	    for(int ii=-1; ii<=1 && test; ii++){
	      for(int jj=-1; jj<=1 && test; jj++){
		for(int kk=-1; kk<=1 && test; kk++){
		  std::cout<<"alpha "<<grille[i+ii][j+jj][k+kk].alpha<< "  "<< " rho "<<grille[i+ii][j+jj][k+kk].rho << "p "<< grille[i+ii][j+jj][k+kk].p<<std::endl; 
		}
	      }
	    }
	  }
	  if(grille[i][j][k].p<0. || grille[i][j][k].rho<0.){
	    test_fini = false;
	  }
	}
      } 
    }
  }
	
  if(!test_fini){
    cout<<" Mixing did not complete "<<endl;
    Mixage();
  }
}

/*!\brief Resolution of the fluid equations.
   \details Directional (Strang) splitting at each time-step.
   \param t Current simulation time
   \param dt Time-step
   \param n  index of the time iterations
   \return void
*/
void Grille::Solve(const double dt, double t, int n, Solide& S){
    
  for(int i=0;i<Nx+2*marge;i++){
    for(int j=0;j<Ny+2*marge;j++){
      for(int k=0;k<Nz+2*marge;k++){
	Cellule  c = grille[i][j][k];
	c.rho0 = c.rho;
	c.impx0 = c.impx;
	c.impy0 = c.impy;
	c.impz0 = c.impz;
	c.rhoE0 = c.rhoE;
	c.p1=c.p;
	c.alpha0=c.alpha;
	c.kappai0 = c.kappai; c.kappaj0 = c.kappaj; c.kappak0 = c.kappak;
	grille[i][j][k] = c;
      }
    }
  }
   
  //Directional splitting
	
  if(n%6==0){

    fnumx(dt/dx,t);     
    solve_fluidx(dt);  
    BC();           
    Fill_cel(S);
      
    fnumy(dt/dy,t);     
    solve_fluidy(dt);  
    BC();           
    Fill_cel(S);
      
      
    fnumz(dt/dz,t);     
    solve_fluidz(dt);  
    BC();           
    Fill_cel(S);
      
  } 
  else if(n%6==2){
    fnumx(dt/dx,t);     
    solve_fluidx(dt);  
    BC();           
    Fill_cel(S);
      
    fnumz(dt/dz,t);     
    solve_fluidz(dt);  
    BC();           
    Fill_cel(S);
         
    fnumy(dt/dy,t);     
    solve_fluidy(dt);  
    BC();           
    Fill_cel(S);
          
  } 
    
  else if(n%6==1){
      
    fnumy(dt/dy,t);     
    solve_fluidy(dt);  
    BC();           
    Fill_cel(S);
          
    fnumx(dt/dx,t);     
    solve_fluidx(dt);  
    BC();           
    Fill_cel(S);
	    
    fnumz(dt/dz,t);     
    solve_fluidz(dt);  
    BC();           
    Fill_cel(S);
	
  } 
    
  else if(n%6==3){
        
    fnumy(dt/dy,t);     
    solve_fluidy(dt);  
    BC();           
    Fill_cel(S);
		      
    fnumz(dt/dz,t);     
    solve_fluidz(dt);  
    BC();           
    Fill_cel(S);
		      
    fnumx(dt/dx,t);     
    solve_fluidx(dt);  
    BC();           
    Fill_cel(S);
		      
  }
  else if(n%6==4){
        
    fnumz(dt/dz,t);     
    solve_fluidz(dt);  
    BC();           
    Fill_cel(S);
		      
    fnumx(dt/dx,t);     
    solve_fluidx(dt);  
    BC();           
    Fill_cel(S);
		      
    fnumy(dt/dy,t);     
    solve_fluidy(dt);  
    BC();           
    Fill_cel(S);
	   
        
  }
  else if(n%6==5){
        
    fnumz(dt/dz,t);     
    solve_fluidz(dt);  
    BC();           
    Fill_cel(S);
		     
    fnumy(dt/dy,t);     
    solve_fluidy(dt);  
    BC();           
    Fill_cel(S);
		      
    fnumx(dt/dx,t);     
    solve_fluidx(dt);  
    BC();           
    Fill_cel(S);
		      
  }
  
    
}


/*!\brief Filling of the ghost cells (\a alpha = 1)
   \details In order to compute the fluxes at the fluid-solid interface, we define a fictitious state in the cells fully occupied by solid (\a alpha = 1), which is taken equal to the state in the mirror cell with regards to the boundary. \n
 Algorithm: search for the interface closest to the center of the cell (loop on all solid faces) and compute the projection of the cell center on this interface using function <b> CGAL::projection(Point_3) </b>.
 	\param S  Solid 
 	\warning <b> Specific coupling procedure ! </b>
 	\return void
 */

void Grille::Fill_cel(Solide& S){
	
  int nb_part = S.size();
  double dist[6*nb_part];
  double dist_min = 100;
  int poz=0;
  double x_min=0., y_min=0., z_min = 0., x_max = 0., y_max=0., z_max=0.;
  Bbox Fluide(X0,Y0,Z0,X0+domainex,Y0+domainey,Z0+domainez);

  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      for(int k=marge;k<Nz+marge;k++){
	Triangle_3 Tri;
	Cellule& c = grille[i][j][k];
	if((std::abs(c.alpha-1.)<eps) ){
	  Point_3 center_cell(c.x, c.y, c.z);
	  int nbx=0, nby=0,nbz=0;
	  Point_3 projete(0.,0.,0.); //Projection on the closest face
	  Vector_3 V_f(0.,0.,0.); //Velocity of the solid boundary at the projected point
	  double dist_min = 10000000.;
	  bool fluide = false;
	  int cas = 0;
	  Point_3 triangle1;
	  Point_3 triangle2;
	  Point_3 triangle3;
	  for(int iter=0; iter<nb_part; iter++){
	    for(int it=0;it<S.solide[iter].triangles.size();it++){
	      if(S.solide[iter].fluide[it]){
		Plane_3 P(S.solide[iter].triangles[it].operator[](0),S.solide[iter].triangles[it].operator[](1),S.solide[iter].triangles[it].operator[](2));
		Point_3 xP = P.projection(center_cell);
		//Test whether the projection is on the face
		bool test = true;
		for(int k=0;k<3 && test;k++){
		  int kp = (k+1)%3;
		  Point_3 x1 = S.solide[iter].triangles[it].operator[](k);
		  Point_3 x2 = S.solide[iter].triangles[it].operator[](kp);
		  Vector_3 vect1(xP,x1);
		  Vector_3 vect2(xP,x2);
		  if(CGAL::to_double(CGAL::cross_product(vect1,vect2)*S.solide[iter].normales[it])<0.){
		    test = false;
		  }
		}
		//First case: the point is on the face
		if(test){
		  double d = sqrt(CGAL::to_double(CGAL::squared_distance(center_cell,xP)));
		  if(d<dist_min && inside_box(Fluide,xP)){
		    dist_min = d;
		    projete = xP;
		    V_f = S.solide[iter].vitesse_parois(xP);
		    fluide = S.solide[iter].fluide[it];
		    cas = 1;
		    triangle1=S.solide[iter].triangles[it].operator[](0);
		    triangle2=S.solide[iter].triangles[it].operator[](1);
		    triangle3=S.solide[iter].triangles[it].operator[](2);
		  }
		}
		//Second case: the point is out of the face
		else{
		  //Search for the closest point on all edges
		  for(int k=0;k<3;k++){
		    int kp = (k+1)%3;
		    Point_3 x1 = S.solide[iter].triangles[it].operator[](k);
		    Point_3 x2 = S.solide[iter].triangles[it].operator[](kp);
		    double d1 = sqrt(CGAL::to_double(CGAL::squared_distance(center_cell,x1)));
		    double d2 = sqrt(CGAL::to_double(CGAL::squared_distance(center_cell,x2)));
		    double d12 = sqrt(CGAL::to_double(CGAL::squared_distance(x1,x2)));
		    //First subcase: the closest point is x1
		    if(d1*d1+d12*d12<d2*d2){
		      if(d1<dist_min && inside_box(Fluide,x1)){
			dist_min = d1;
			projete = x1;
			V_f = S.solide[iter].vitesse_parois(x1);
			fluide = S.solide[iter].fluide[it];
			cas = 2;
		      }
		    }
		    //Second subcase: the closest point is x2
		    else if(d2*d2+d12*d12<d1*d1){
		      if(d2<dist_min && inside_box(Fluide,x2)){
			dist_min = d2;
			projete = x2;
			V_f = S.solide[iter].vitesse_parois(x2);
			fluide = S.solide[iter].fluide[it];
			cas = 3;
		      }
		    }
		    //Third subcase: take the projection on (x1,x2)
		    else {
		      Line_3 L(x1,x2);
		      double d = sqrt(CGAL::to_double(CGAL::squared_distance(center_cell,L)));
		      Point_3 proj = L.projection(center_cell);
		      if(d<dist_min && inside_box(Fluide,proj)){
			dist_min = d;
			projete = proj;
			V_f = S.solide[iter].vitesse_parois(proj);
			fluide = S.solide[iter].fluide[it];
			cas = 4;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  //Computation of the symmetric point with regards to the plan defined by centre_face and normale_face
	  Point_3 symm_center = center_cell + Vector_3(center_cell,projete)*2;
	  Vector_3 normale(center_cell,projete);
	  double norme = sqrt(CGAL::to_double(normale.squared_length()));
	  assert(norme!= 0.);
	  normale = normale*1./norme;
	  const Cellule& cm = in_cell(symm_center);
	  Vector_3 vit_m(cm.u,cm.v,cm.w); //Velocity at the mirror point
	  Vector_3 vit = vit_m - normale*2.*((vit_m-V_f)*normale);
	  if(abs(cm.alpha-1.)<eps){
	    cout << "solid target cell: original=" << c.x << " " << c.y << " " << c.z << " target=" << cm.x << " " << cm.y << " " << cm.z << " projection=" << projete.x() << " " << projete.y() << " " << projete.z() << " fluid=" << fluide << " case=" <<  cas << " triangle=" << triangle1.x() << " " << triangle1.y() << " " << triangle1.z() << " " << triangle2.x() << " " << triangle2.y() << " " << triangle2.z() << " " << triangle3.x() << " " << triangle3.y() << " " << triangle3.z() << " " << endl;
	  }
		
	  c.rho = cm.rho;
	  c.u = CGAL::to_double(vit.operator[](0));
	  c.v = CGAL::to_double(vit.operator[](1));
	  c.w = CGAL::to_double(vit.operator[](2));
	  c.p = cm.p;
	  c.impx = c.rho*c.u;
	  c.impy = c.rho*c.v;
	  c.impz = c.rho*c.w;
	  if(std::abs(2.*(c.u*c.u+c.v*c.v+c.w*c.w)+c.p/(gam-1.)) > eps_vide){
	    c.rhoE = c.rho/2.*(c.u*c.u+c.v*c.v+c.w*c.w)+c.p/(gam-1.);
	  }
	  if( (std::abs(c.rho) <= eps_vide ) || (std::abs(c.p)<= eps_vide) ){
	    c.vide = true;
	  }
	  else{c.vide = false;}
	}
      }
    }
  }
}

/*!\brief Computation of the signed volume of a prism.
  \details The signed volume of the prism with bases the triangles \f$ T1(A_1,B_1,C_1)  \f$ and \f$ T2(A_2,B_2,C_2) \f$ is given by: \n
 \f{eqnarray*}{
 {\Vert A_1 B_1 C_1 A_2 B_2 C_2 \Vert}_P  = \frac{1}{36}  \left( 2 \vec{A_1 B_1} \wedge \vec{A_1 C_1} + 2 \vec{A_2 B_2} \wedge \vec{A_2 C_2} +  \vec{A_1 B_1}\wedge \vec{A_2 C_2}  +  \vec{A_2 B_2}\wedge \vec{A_1 C_1} \right) \cdot
 \f} 
 \f{eqnarray*}{  \left( \vec{A_1 A_2}  + \vec{B_1 B_2} + \vec{C_1 C_2}\right)\f}
 	\param T1 Triangle_3  basis of the prism
 	\param T2 Triangle_3  basis of prism
 	\warning <b> Specific coupling procedure ! </b>
 	\return double
 */
double volume_prisme(const Triangle_3& T1,const Triangle_3& T2){
	
  double volume=0.;
	
  Vector_3 V = 2.*cross_product( Vector_3(T1.operator[](0),T1.operator[](1)), Vector_3(T1.operator[](0),T1.operator[](2)) )
    + 2*cross_product( Vector_3(T2.operator[](0),T2.operator[](1)), Vector_3(T2.operator[](0),T2.operator[](2)) )
    + cross_product( Vector_3(T1.operator[](0),T1.operator[](1)), Vector_3(T2.operator[](0),T2.operator[](2)) )
    + cross_product( Vector_3(T2.operator[](0),T2.operator[](1)), Vector_3(T1.operator[](0),T1.operator[](2)) );
	
  volume = CGAL::to_double((Vector_3(T1.operator[](0),T2.operator[](0)) + Vector_3(T1.operator[](1),T2.operator[](1)) + 
			    Vector_3(T1.operator[](2),T2.operator[](2)))*V);
	
  volume /=36;
  return volume;
}

/*!\brief Computation of the signed colume of a tetrahedron.
  \details The signed volume of a tetrahedron \f$ T(A,B,C,D)  \f$ is given by: \n
 \f{eqnarray*}{
 {\Vert A B C D \Vert}_{sign} = \frac{1}{6} \vec{A D} \cdot \left(  \vec{A B} \wedge \vec{A C}  \right) 
 \f}
 	\param Tet Tetrahedron
 	\warning <b> Specific coupling procedure ! </b>
 	\return double
 */
double volume_tetra(const Tetrahedron& Tet){
  return CGAL::to_double(Tet.volume());
}

/*!\brief Barycentric transformation of point Xn.
  \details Let Xn a point belonging to triangle Tn(A_1,B_1,C_1). The barycentric transform of Xn is the point Xn1 (belonging to triangle Tn1(A_2,B_2,C_2) ) given by: 
 \f{eqnarray*}{
 \lambda = \frac{\Vert  \vec{C_1 Xn} \wedge \vec{C_1 B_1}  \Vert }{\Vert  \vec{C_1 A_1} \wedge \vec{C_1 B_1}  \Vert} 
 \f}
 \f{eqnarray*}{
 \mu = \frac{\Vert  \vec{C_1 Xn} \wedge \vec{C_1 A_1}  \Vert }{\Vert  \vec{C_1 B_1} \wedge \vec{C_1 A_1}  \Vert} 
 \f}
 \f{eqnarray*}{
 Xn1 = \lambda A_2 + \mu B_2 + (1-\lambda -\mu)C_2 
 \f}
 	\param Tn Triangle_3
 \param Tn1 Triangle_3 
 \param Xn  Point_3
 	\warning <b> Specific coupling procedure ! </b>
 	\return Point_3
 */
Point_3 tr(const Triangle_3& Tn, const Triangle_3& Tn1, const Point_3& Xn){
  
  
  double lambda = 0., mu = 0.;
	
  double AC2 = CGAL::to_double(Vector_3(Tn.vertex(0),Tn.vertex(2)).squared_length());
  double BC2 = CGAL::to_double(Vector_3(Tn.vertex(1),Tn.vertex(2)).squared_length());
  double AC_BC = CGAL::to_double(Vector_3(Tn.vertex(0),Tn.vertex(2))*Vector_3(Tn.vertex(1),Tn.vertex(2)));
  double XC_AC = CGAL::to_double(Vector_3(Xn,Tn.vertex(2))*Vector_3(Tn.vertex(0),Tn.vertex(2)));
  double XC_BC = CGAL::to_double(Vector_3(Xn,Tn.vertex(2))*Vector_3(Tn.vertex(1),Tn.vertex(2)));
  double dom = AC2*BC2-AC_BC*AC_BC;
  double num1 = BC2*XC_AC-AC_BC*XC_BC;
  double num2 = AC2*XC_BC-AC_BC*XC_AC;


  lambda =  num1/dom;
  mu = num2/dom;
	
  double x = CGAL::to_double(lambda * Tn1.operator[](0).operator[](0) + mu*Tn1.operator[](1).operator[](0) 
			     + (1- lambda- mu)*Tn1.operator[](2).operator[](0));
	
  double y = CGAL::to_double(lambda * Tn1.operator[](0).operator[](1) + mu*Tn1.operator[](1).operator[](1) 
			     + (1- lambda- mu)*Tn1.operator[](2).operator[](1));
		 
  double z = CGAL::to_double(lambda * Tn1.operator[](0).operator[](2) + mu*Tn1.operator[](1).operator[](2) 
			     + (1- lambda- mu)*Tn1.operator[](2).operator[](2));

  return Point_3(x, y, z);
}

/*!\brief Barycentric transformation of Triangle T.
 \details Calls function \a tr(Triangle_3, Triangle_3, Point_3) for each triangle vertex.
 \param Tn Triangle_3
 \param Tn1 Triangle_3 
 \param T  Triangle_3 
 \warning <b> Specific coupling procedure ! </b>
 \return Triangle_3 
 */
Triangle_3 tr(const Triangle_3& Tn, const Triangle_3& Tn1, const Triangle_3& T){
	
  const Point_3& s = tr( Tn,Tn1, T.operator[](0) );
  const Point_3& r = tr( Tn,Tn1, T.operator[](1) );
  const Point_3& v = tr( Tn,Tn1, T.operator[](2) );
  
  return Triangle_3(s, r, v);
}


/*!\brief Transformation of a Point_3 into a Point_2.
  \details Let Triangle Tn1(A,B,C), the barycentric transform of Point_3 Xn (belonging to Tn1) into a Point_2 is given by: 

 \f{eqnarray*}{
 \lambda = \frac{\Vert  \vec{C Xn} \wedge \vec{C B}  \Vert }{\Vert  \vec{C A} \wedge \vec{C B}  \Vert} 
 \f}
 \f{eqnarray*}{
 \mu = \frac{\Vert  \vec{C A} \wedge \vec{C Xn}  \Vert }{\Vert  \vec{C A} \wedge \vec{C B}  \Vert} 
 \f}
 \f{eqnarray*}{
 X_{2d} = (\mu, (1-\lambda-\mu))
 \f}
 \param Tn1 Triangle_3 
 \param Xn  Point_3
 \warning <b> Specific coupling procedure ! </b>
 \return Point_2
 */
Point_2 tr(const Triangle_3& Tn1, const Point_3& Xn){
		
  double lambda = 0., mu = 0.;

  double AC2 = CGAL::to_double(Vector_3(Tn1.vertex(0),Tn1.vertex(2)).squared_length());
  double BC2 = CGAL::to_double(Vector_3(Tn1.vertex(1),Tn1.vertex(2)).squared_length());
  double AC_BC = CGAL::to_double(Vector_3(Tn1.vertex(0),Tn1.vertex(2))*Vector_3(Tn1.vertex(1),Tn1.vertex(2)));
  double XC_AC = CGAL::to_double(Vector_3(Xn,Tn1.vertex(2))*Vector_3(Tn1.vertex(0),Tn1.vertex(2)));
  double XC_BC = CGAL::to_double(Vector_3(Xn,Tn1.vertex(2))*Vector_3(Tn1.vertex(1),Tn1.vertex(2)));
  double dom = AC2*BC2-AC_BC*AC_BC;
  double num1 = BC2*XC_AC-AC_BC*XC_BC;
  double num2 = AC2*XC_BC-AC_BC*XC_AC;

												 
 
  lambda =  num1/dom;
  mu = num2/dom;

  Point_2 M(mu, (1-lambda-mu));
	
  return M;
}	

/*!\brief Transformation of a Triangle_3 into a Triangle_2 
 \details Calls function \a tr(Triangle_3, Point_3) for each triangle vertex.
 \param Tn1 Triangle_3 
 \param T  Triangle_3 
 \warning <b> Specific coupling procedure ! </b>
 \return Triangle_2 
 */
Triangle_2 tr(const Triangle_3& Tn1, const Triangle_3& T){
	
  const Point_2& s = tr( Tn1, T.operator[](0) );
  const Point_2& r = tr( Tn1, T.operator[](1) );
  const Point_2& v = tr( Tn1, T.operator[](2) );
	
  return Triangle_2(s, r, v);
}


/*!\brief Barycentric transformation of a Point_2 into a Point_3.
  \details  Let triangle Tn1(A,B,C), the transform of Point_2 Xn(X,Y) into a point_3 (belonging to Tn1) is given by: 
 \f{eqnarray*}{
 \lambda = 1- X - Y
 \f}
 \f{eqnarray*}{
 \mu = X
 \f}
 \f{eqnarray*}{
 X_{3d} = \lambda A + \mu B + (1-\lambda -\mu)C 
 \f}
 \param Tn1 Triangle_3 
 \param Xn  Point_2
 \warning <b> Specific coupling procedure ! </b>
 \return Point_3
 */
Point_3 tr(const Triangle_3& Tn1, const Point_2& Xn){
  
  double lambda = CGAL::to_double(1.- Xn.operator[](0) -  Xn.operator[](1));
  double mu = CGAL::to_double(Xn.operator[](0));
  
  double x = CGAL::to_double(lambda * Tn1.operator[](0).operator[](0) + mu*Tn1.operator[](1).operator[](0) 
			     + (1- lambda- mu)*Tn1.operator[](2).operator[](0));
  
  double y = CGAL::to_double(lambda * Tn1.operator[](0).operator[](1) + mu*Tn1.operator[](1).operator[](1) 
			     + (1- lambda- mu)*Tn1.operator[](2).operator[](1));
  double z = CGAL::to_double(lambda * Tn1.operator[](0).operator[](2) + mu*Tn1.operator[](1).operator[](2) 
			     + (1- lambda- mu)*Tn1.operator[](2).operator[](2));		
  
  
  return Point_3(x,y,z);
}	
/*!\brief  Barycentric transformation of a Triangle_2 into a Triangle_3
 \details Calls function tr(Triangle_3, Point_2) for each triangle vertex.
 \param Tn1 Triangle_3 
 \param T  Triangle_2 
 \warning <b> Specific coupling procedure ! </b>
 \return Triangle_3 
 */
Triangle_3 tr(const Triangle_3& Tn1, const Triangle_2& T){
  
  const Point_3& s = tr( Tn1, T.operator[](0) );
  const Point_3& r = tr( Tn1, T.operator[](1) );
  const Point_3& v = tr( Tn1, T.operator[](2) );
  
  return Triangle_3(s, r, v);
}



/*!\brief List of fluid cells intersected by an interface triangle between times t-dt and t.
 \details Let
 \f$ C_1= grille[in][jn][kn]\f$ the cell where the interface triangle lies at time t and \f$ C_2= grille[in1][jn1][kn1]\f$ the cell where it lies at time t-dt. \n
 Due to the fluid CFL condition, \f$ max(|in-in1 |, |jn-jn1 |, |kn-kn1 |)<=1 \f$. 
 There are therefore at most 8 cells intersected by the prism with bases the interface triangle at times t-dt (\a Particule.triangles_prev)  and t (\a Particule.triangles) ):\n
 -\f$ grille[in][jn][kn]\f$ \n
 -\f$ grille[in1][jn][kn]\f$ \n
 -\f$ grille[in][jn1][kn]\f$ \n
 -\f$ grille[in][jn1][kn1]\f$ \n
 -\f$ grille[in1][jn1][kn]\f$ \n
 -\f$ grille[in1][jn][kn1]\f$ \n
 -\f$ grille[in][jn1][kn1]\f$ \n
 -\f$ grille[in1][jn1][kn1]\f$ \n

 \param (in,jn,kn) index of a Cellule (cell where the interface triangle lies at time t-dt (\a Particule.triangles_prev))
 \param (in1,jn1,kn1) index of a Cellule (cell where the interface triangle lies at time t (\a Particule.triangles))
 \param box_cells vector of Box_3
 \param Cells vector of Cellule
 \warning <b> Specific coupling procedure ! </b>
 \return void
 */
void Grille::cells_intersection_face(int& in,int& jn,int& kn,int& in1,int& jn1,int& kn1, std::vector<Bbox>& box_cells, std::vector<Cellule>& Cells){
	
  if((in!=in1 && jn==jn1 && kn==kn1)|| (in==in1 && jn!=jn1 && kn==kn1) || (in==in1 && jn==jn1 && kn!=kn1))
  {
    Cellule& c0 = grille[in][jn][kn];
    Cellule& c1 = grille[in1][jn1][kn1];
    Cells.push_back(c0); Cells.push_back(c1);
    box_cells.push_back(Bbox(c0.x -c0.dx/2.,c0.y -c0.dy/2.,c0.z -c0.dz/2.,
			     c0.x +c0.dx/2.,c0.y +c0.dy/2.,c0.z + c0.dz/2.));
    box_cells.push_back(Bbox(c1.x -c1.dx/2.,c1.y -c1.dy/2.,c1.z -c1.dz/2.,
			     c1.x +c1.dx/2.,c1.y +c1.dy/2.,c1.z + c1.dz/2.));
  }
  else if(in!=in1 && jn!=jn1 && kn==kn1){
		
    Cellule& c0 = grille[in][jn][kn];
    Cellule& c1 = grille[in1][jn1][kn1];
    Cellule& c2 = grille[in1][jn][kn];
    Cellule& c3 = grille[in][jn1][kn1];
		
    Cells.push_back(c0); Cells.push_back(c1); Cells.push_back(c2); Cells.push_back(c3);
		
    box_cells.push_back(Bbox(c0.x -c0.dx/2.,c0.y -c0.dy/2.,c0.z -c0.dz/2.,
			     c0.x +c0.dx/2.,c0.y +c0.dy/2.,c0.z + c0.dz/2.));
    box_cells.push_back(Bbox(c1.x -c1.dx/2.,c1.y -c1.dy/2.,c1.z -c1.dz/2.,
			     c1.x +c1.dx/2.,c1.y +c1.dy/2.,c1.z + c1.dz/2.));
    box_cells.push_back(Bbox(c2.x -c2.dx/2.,c2.y -c2.dy/2.,c2.z -c2.dz/2.,
			     c2.x +c2.dx/2.,c2.y +c2.dy/2.,c2.z + c2.dz/2.));
    box_cells.push_back(Bbox(c3.x -c3.dx/2.,c3.y -c3.dy/2.,c3.z -c3.dz/2.,
			     c3.x +c3.dx/2.,c3.y +c3.dy/2.,c3.z + c3.dz/2.));
  }
  else if(in!=in1 && jn==jn1 && kn!=kn1){
		
    Cellule& c0 = grille[in][jn][kn];
    Cellule& c1 = grille[in1][jn1][kn1];
    Cellule& c2 = grille[in][jn][kn1];
    Cellule& c3 = grille[in1][jn1][kn];
		
    Cells.push_back(c0); Cells.push_back(c1); Cells.push_back(c2); Cells.push_back(c3);
		
    box_cells.push_back(Bbox(c0.x -c0.dx/2.,c0.y -c0.dy/2.,c0.z -c0.dz/2.,
			     c0.x +c0.dx/2.,c0.y +c0.dy/2.,c0.z + c0.dz/2.));
    box_cells.push_back(Bbox(c1.x -c1.dx/2.,c1.y -c1.dy/2.,c1.z -c1.dz/2.,
			     c1.x +c1.dx/2.,c1.y +c1.dy/2.,c1.z + c1.dz/2.));
    box_cells.push_back(Bbox(c2.x -c2.dx/2.,c2.y -c2.dy/2.,c2.z -c2.dz/2.,
			     c2.x +c2.dx/2.,c2.y +c2.dy/2.,c2.z + c2.dz/2.));
    box_cells.push_back(Bbox(c3.x -c3.dx/2.,c3.y -c3.dy/2.,c3.z -c3.dz/2.,
			     c3.x +c3.dx/2.,c3.y +c3.dy/2.,c3.z + c3.dz/2.));
  }
  else if(in==in1 && jn!=jn1 && kn!=kn1){
		
    Cellule& c0 = grille[in][jn][kn];
    Cellule& c1 = grille[in1][jn1][kn1];
    Cellule& c2 = grille[in][jn1][kn];
    Cellule& c3 = grille[in1][jn][kn1];
		
    Cells.push_back(c0); Cells.push_back(c1); Cells.push_back(c2); Cells.push_back(c3);
		
    box_cells.push_back(Bbox(c0.x -c0.dx/2.,c0.y -c0.dy/2.,c0.z -c0.dz/2.,
			     c0.x +c0.dx/2.,c0.y +c0.dy/2.,c0.z + c0.dz/2.));
    box_cells.push_back(Bbox(c1.x -c1.dx/2.,c1.y -c1.dy/2.,c1.z -c1.dz/2.,
			     c1.x +c1.dx/2.,c1.y +c1.dy/2.,c1.z + c1.dz/2.));
    box_cells.push_back(Bbox(c2.x -c2.dx/2.,c2.y -c2.dy/2.,c2.z -c2.dz/2.,
			     c2.x +c2.dx/2.,c2.y +c2.dy/2.,c2.z + c2.dz/2.));
    box_cells.push_back(Bbox(c3.x -c3.dx/2.,c3.y -c3.dy/2.,c3.z -c3.dz/2.,
			     c3.x +c3.dx/2.,c3.y +c3.dy/2.,c3.z + c3.dz/2.));
  }
  else{
	
    Cellule& c0 = grille[in][jn][kn];
    Cellule& c1 = grille[in1][jn1][kn1];
    Cellule& c2 = grille[in1][jn][kn];
    Cellule& c3 = grille[in][jn1][kn];
    Cellule& c4 = grille[in][jn][kn1];
    Cellule& c5 = grille[in1][jn1][kn];
    Cellule& c6 = grille[in1][jn][kn1];
    Cellule& c7 = grille[in][jn1][kn1];
		
		
    Cells.push_back(c0); Cells.push_back(c1); Cells.push_back(c2); Cells.push_back(c3);
    Cells.push_back(c4); Cells.push_back(c5); Cells.push_back(c6); Cells.push_back(c7);
		
    box_cells.push_back(Bbox(c0.x -c0.dx/2.,c0.y -c0.dy/2.,c0.z -c0.dz/2.,
			     c0.x +c0.dx/2.,c0.y +c0.dy/2.,c0.z + c0.dz/2.));
    box_cells.push_back(Bbox(c1.x -c1.dx/2.,c1.y -c1.dy/2.,c1.z -c1.dz/2.,
			     c1.x +c1.dx/2.,c1.y +c1.dy/2.,c1.z + c1.dz/2.));
    box_cells.push_back(Bbox(c2.x -c2.dx/2.,c2.y -c2.dy/2.,c2.z -c2.dz/2.,
			     c2.x +c2.dx/2.,c2.y +c2.dy/2.,c2.z + c2.dz/2.));
    box_cells.push_back(Bbox(c3.x -c3.dx/2.,c3.y -c3.dy/2.,c3.z -c3.dz/2.,
			     c3.x +c3.dx/2.,c3.y +c3.dy/2.,c3.z + c3.dz/2.));
    box_cells.push_back(Bbox(c4.x -c4.dx/2.,c4.y -c4.dy/2.,c4.z -c4.dz/2.,
			     c4.x +c4.dx/2.,c4.y +c4.dy/2.,c4.z + c4.dz/2.));
    box_cells.push_back(Bbox(c5.x -c5.dx/2.,c5.y -c5.dy/2.,c5.z - c5.dz/2.,
			     c5.x +c5.dx/2.,c5.y +c5.dy/2.,c5.z + c5.dz/2.));
    box_cells.push_back(Bbox(c6.x -c6.dx/2.,c6.y -c6.dy/2.,c6.z -c6.dz/2.,
			     c6.x +c6.dx/2.,c6.y +c6.dy/2.,c6.z + c6.dz/2.));
    box_cells.push_back(Bbox(c7.x -c7.dx/2.,c7.y -c7.dy/2.,c7.z -c7.dz/2.,
			     c7.x +c7.dx/2.,c7.y +c7.dy/2.,c7.z + c7.dz/2.));
  }

}
/*!\brief Computation of the quantity of fluid swept by an interface triangle between t-dt and t, and computation of the boundary flux.
 \details Algorithm:\n
 - Construction of the Box vector containing the prisms with bases T3d_prev and T3d_n. Loop on the prisms obtained:
 - We look for the indices of the cells containing T3d_prev (constructed to be fully contained in one cell) and T3d_n (idem).
 - If the prism is contained in one single cell, compute the volume of the prism using function \a volume_prisme(const Triangle_3&,const Triangle_3&), and the swept quantity for the interface triangle \a Particule.triangles is given by: \f$  volume\_prisme*U^n/volume\_cellule \f$. Otherwise,
 - List the fluid cells intersected by the prism using function \a cells_intersection_face(int& ,int& ,int& ,int& ,int& ,int& , std::vector<Bbox>& s, std::vector<Cellule>& s).
 - Split the prism into tetrahedra: let  \f$ T1(A_1,B_1,C_1)\f$  and \f$ T2(A_2,B_2,C_2)\f$  the prism bases, define points: \f$ A = \frac{1}{4}(B_1 + B_2 + C_1 +C_2) \f$ , \f$ B = \frac{1}{4}(A_1 + A_2 + C_1 +C_2) \f$ and \f$ C = \frac{1}{4}(A_1 + A_2 + B_1 + B_2 ) \f$. The prism is split into the tetrahedra \f$ A_1,B_1,C_1 A_2,B_2,C_2 \f$ sont: \f$ A_1 A_2 C B \f$, \f$ B_1 B_2 A C \f$, \f$ C_1 C_2 B A \f$, \f$ A_1 C C_1 B \f$, \f$ B_1 A C_1 C \f$, \f$ A C B C_1 \f$, \f$ A B C C_2 \f$, \f$ A B_2 C_2 C \f$, \f$ A_1 B_1 C_1 C \f$, \f$ A_2 C_2 C B \f$, \f$ A_2 B_2 C C_2. \f$
 - Intersect these tetrahedra with the fluid cells intersected by the prism using function \a intersect_cube_tetrahedron(Bbox&, Tetrahedron&). The swept quantity for the interface triangle is given by the sum of the following quantities: \f$  volume\_{intersection\_cellule\_tetrahedre}*U^n/volume\_cellule. \f$ \n

 Computation of the boundary fluxes: let \a f an interface triangle, the boundary flux is given by:
 \f{eqnarray*}{
 \Phi_f  =  \left(0,  p^x \, A_f n^{x}_f, \, p^y \,A_f n^{y}_f, \, p^z \,A_f n^{z}_f, V_f \cdot \left( p^x \, A_f n^{x}_f,p^y \,A_f n^{y}_f,p^z \,A_f n^{z}_f \right)^t \right)^t
 \f} \n
 where \f$ A_f \f$ is the area of f,  \f$ n_f \f$ is the exterior normal to f, \f$ V_f \f$ is the velocity at the center of f computed with function \a vitesse_parois(Point_3& ) and \f$ p^x, p^y, p^z \f$ are the effective pressures in the x, y and z direction during the time-step (\a Cellule.pdtx, \a  Cellule.pdty and \a Cellule.pdtz).

 \param T3d_prev Triangles_3 (interface triangles at time t: \a Particule.Triangles_interface)
 \param T3d_n    Triangles_3 (interface triangles at time t-dt: \a Particule.Triangles_interface_prev)
 \param dt Time-step
 \param P Particule 
 \warning <b> Specific coupling procedure ! </b>
 \return void
 */
void Grille::swap_face(const Triangles& T3d_prev, const Triangles& T3d_n, const double dt, Particule & P, double & volume_test){
	
  std::vector<Bbox> box_prismes(T3d_prev.size());
  for (int i=0; i< T3d_prev.size(); i++){
    const Bbox& box_triangles_prev = T3d_prev[i].bbox();
    const Bbox& box_triangles_n = T3d_n[i].bbox();
    box_prismes[i]= box_triangles_prev + box_triangles_n;
  } 
	
  for (int i=0; i< box_prismes.size(); i++){
    //double vol_test=0.; 
    int in=0, jn=0, kn=0, in1=0, jn1=0, kn1=0;
    bool interieur = true;
    const Point_3& center_prev= centroid(T3d_prev[i].operator[](0),T3d_prev[i].operator[](1),T3d_prev[i].operator[](2));
    const Point_3& center_n= centroid(T3d_n[i].operator[](0),T3d_n[i].operator[](1),T3d_n[i].operator[](2));
    in_cell(center_prev, in, jn, kn, interieur);
    in_cell(center_n, in1, jn1, kn1, interieur);
	
    if((std::abs(grille[in1][jn1][kn1].alpha -1.)<eps)  && (interieur==true)){
      double x= CGAL::to_double(center_n.operator[](0));
      double y= CGAL::to_double(center_n.operator[](1));
      double z= CGAL::to_double(center_n.operator[](2));
      Cellule& cd= grille[in1+1][jn1][kn1];
      if (cd.is_in_cell(x,y,z) ) {in1=in1+1;}
      else{
	Cellule& cg= grille[in1-1][jn1][kn1];
	if (cg.is_in_cell(x,y,z)) {in1=in1-1;}
	else{
	  Cellule& ch= grille[in1][jn1+1][kn1];
	  if (ch.is_in_cell(x,y,z)) {jn1=jn1+1;}
	  else{
	    Cellule& cb= grille[in1][jn1-1][kn1];
	    if (cb.is_in_cell(x,y,z)) {jn1=jn1-1;}
	    else{
	      Cellule& cd= grille[in1][jn1][kn1+1];
	      if (cd.is_in_cell(x,y,z)) {kn1=kn1+1;}
	      else{
		Cellule& cder= grille[in1][jn1][kn1-1];
		if (cder.is_in_cell(x,y,z)) {kn1=kn1-1;}
	      }
	    }
	  }
	}
      }
    } 
		  
    Cellule& c= grille[in1][jn1][kn1];
    double volume_cel = c.dx*c.dy*c.dz;  
    if ( (in==in1) && (jn==jn1) && (kn==kn1) && (interieur==true)){
      //The prism is contained in one single cell
      double volume_p=volume_prisme(T3d_prev[i],T3d_n[i]);
      //Computation of the volume
      if( (std::abs(volume_p)>eps) && (std::abs(1.-c.alpha)>eps)){
	c.delta_w[0] += volume_p*c.rho0/volume_cel; 
	c.delta_w[1] += volume_p*c.impx0/volume_cel;
	c.delta_w[2] += volume_p*c.impy0/volume_cel; 
	c.delta_w[3] += volume_p*c.impz0/volume_cel; 
	c.delta_w[4] += volume_p*c.rhoE0/volume_cel;
      }
      volume_test += volume_p;
    }	
    else if((std::abs(volume_prisme(T3d_prev[i],T3d_n[i])) >eps)  && (interieur==true) && (std::abs(1.-c.alpha)>eps)) {
      std::vector<Bbox> box_cells;
      std::vector<Cellule> Cells ;
      cells_intersection_face(in,jn,kn,in1,jn1,kn1,box_cells,Cells);

      //Definition of the tetrahedra 
      std::vector<Tetrahedron> vect_Tet;
      std::vector<Bbox> box_Tet;
      Point_3 e = centroid(T3d_prev[i].operator[](1),T3d_n[i].operator[](1), T3d_prev[i].operator[](2),T3d_n[i].operator[](2));
      Point_3 f = centroid(T3d_prev[i].operator[](0),T3d_n[i].operator[](0), T3d_prev[i].operator[](2),T3d_n[i].operator[](2));
      Point_3 g = centroid(T3d_prev[i].operator[](0),T3d_n[i].operator[](0), T3d_prev[i].operator[](1),T3d_n[i].operator[](1));

      Tetrahedron tet0 (T3d_prev[i].operator[](0),T3d_n[i].operator[](0), g, f);
      if(abs(tet0.volume ())>eps){
	vect_Tet.push_back(tet0);
	box_Tet.push_back(tet0.bbox());
      }
			
      Tetrahedron tet1 (T3d_prev[i].operator[](1),T3d_n[i].operator[](1), e, g);
      if(abs(tet1.volume ())>eps){
	vect_Tet.push_back(tet1);
	box_Tet.push_back(tet1.bbox());
      }

      Tetrahedron tet2 (T3d_prev[i].operator[](2),T3d_n[i].operator[](2), f, e);
      if(abs(tet2.volume ())>eps){
	vect_Tet.push_back(tet2);
	box_Tet.push_back(tet2.bbox());
      }

      Tetrahedron tet3 (T3d_prev[i].operator[](0),T3d_prev[i].operator[](1), T3d_prev[i].operator[](2), g);
      if(abs(tet3.volume ())>eps){
	vect_Tet.push_back(tet3);
	box_Tet.push_back(tet3.bbox());
      }

      Tetrahedron tet4 (T3d_prev[i].operator[](0), f, g, T3d_prev[i].operator[](2));
      if(abs(tet4.volume ())>eps){
	vect_Tet.push_back(tet4);
	box_Tet.push_back(tet4.bbox());
      }
      Tetrahedron tet5 (T3d_prev[i].operator[](1), g, e, T3d_prev[i].operator[](2));
      if(abs(tet5.volume ())>eps){
	vect_Tet.push_back(tet5);
	box_Tet.push_back(tet5.bbox());
      }

      Tetrahedron tet6 (e,g,f, T3d_prev[i].operator[](2));
      if(abs(tet6.volume ())>eps){
	vect_Tet.push_back(tet6);
	box_Tet.push_back(tet6.bbox());
      }

      Tetrahedron tet7 (e,f,g, T3d_n[i].operator[](2));
      if(abs(tet7.volume ())>eps){
	vect_Tet.push_back(tet7);
	box_Tet.push_back(tet7.bbox());
      }
      Tetrahedron tet8 (e, T3d_n[i].operator[](1),T3d_n[i].operator[](2), g);
      if(abs(tet8.volume ())>eps){
	vect_Tet.push_back(tet8);
	box_Tet.push_back(tet8.bbox());
      }
      Tetrahedron tet9 (T3d_n[i].operator[](0),T3d_n[i].operator[](2),g,f);
      if(abs(tet9.volume ())>eps){
	vect_Tet.push_back(tet9);
	box_Tet.push_back(tet9.bbox());
      }

      Tetrahedron tet10 (T3d_n[i].operator[](0),T3d_n[i].operator[](1), g, T3d_n[i].operator[](2));
      if(abs(tet10.volume ())>eps){
	vect_Tet.push_back(tet10);
	box_Tet.push_back(tet10.bbox());
      }
		
      for(int iter=0; iter<box_cells.size(); iter++){
			
	double volume = 0.;
			
	if (CGAL::do_intersect(box_prismes[i], box_cells[iter]) ) {
	  //Intersection test of the bounding boxes of the tetrahedra with the 8 neighbouring cells
	  for(int it=0; it<box_Tet.size(); it++){
	    if (CGAL::do_intersect(box_Tet[it],box_cells[iter]) ) {
	      //Test to check whethet vect_Tet[it] is contained in cell "iter"
	      if (inside_box(box_cells[iter], vect_Tet[it].operator[](0)) &&  inside_box(box_cells[iter], vect_Tet[it].operator[](1))
		  && inside_box(box_cells[iter], vect_Tet[it].operator[](2)) && inside_box(box_cells[iter], vect_Tet[it].operator[](3))){
		volume += volume_tetra(vect_Tet[it]); 
	      }
	      else {
		//Computation of the intersection volume
		double temps_intersections=0.,temps_triangulation=0.;
		volume += (intersect_cube_tetrahedron_bis(box_cells[iter], vect_Tet[it],temps_intersections,temps_triangulation) * sign(volume_tetra(vect_Tet[it])) );
	      }
	    } 
	  } 			
	}
			
	c.delta_w[0] += volume*Cells[iter].rho0/volume_cel; 
	c.delta_w[1] += volume*Cells[iter].impx0/volume_cel;
	c.delta_w[2] += volume*Cells[iter].impy0/volume_cel; 
	c.delta_w[3] += volume*Cells[iter].impz0/volume_cel; 
	c.delta_w[4] += volume*Cells[iter].rhoE0/volume_cel;
	grille[in1][jn1][kn1] = c;
				
	volume_test += volume;
      } 
    }
    if (explicite){
      Vector_3 norm_prev= orthogonal_vector(T3d_prev[i].operator[](0),T3d_prev[i].operator[](1),T3d_prev[i].operator[](2));
      double norm2_prev= sqrt(CGAL::to_double(norm_prev*norm_prev));
      if(norm2_prev>eps){ 
	Cellule c_prev= grille[in][jn][kn];
	Vector_3 n_prev = norm_prev/norm2_prev;
	double aire_prev = sqrt(CGAL::to_double(T3d_prev[i].squared_area()));
	c.phi_x += c_prev.pdtx * aire_prev *( CGAL::to_double(n_prev.x()))/volume_cel;
	c.phi_y += c_prev.pdty * aire_prev *( CGAL::to_double(n_prev.y()))/volume_cel;
	c.phi_z += c_prev.pdtz * aire_prev *( CGAL::to_double(n_prev.z()))/volume_cel;
			
	Vector_3 V_f = P.vitesse_parois_prev(center_prev);
	c.phi_v += aire_prev * (CGAL::to_double(c.pdtx*n_prev.x()*V_f.x()  + c.pdty*n_prev.y()*V_f.y()+
						c.pdtz*n_prev.z()*V_f.z()))/volume_cel;
	grille[in1][jn1][kn1] = c;
      }
    }
    else {
      Vector_3 norm= orthogonal_vector(T3d_n[i].operator[](0),T3d_n[i].operator[](1),T3d_n[i].operator[](2));
      double norm2= sqrt(CGAL::to_double(norm*norm));
      if(norm2>eps){ 
	Vector_3 n = norm/norm2;
	double aire = sqrt(CGAL::to_double(T3d_n[i].squared_area()));
	c.phi_x += c.pdtx * aire *( CGAL::to_double(n.x()))/(c.dx*c.dy*c.dz);
	c.phi_y += c.pdty * aire *( CGAL::to_double(n.y()))/(c.dx*c.dy*c.dz);
	c.phi_z += c.pdtz * aire *( CGAL::to_double(n.z()))/(c.dx*c.dy*c.dz);
	Vector_3 V_f = P.vitesse_parois(center_n);
	c.phi_v += aire * (CGAL::to_double(c.pdtx*n.x()*V_f.x()  + c.pdty*n.y()*V_f.y() + c.pdtz*n.z()*V_f.z()))/(c.dx*c.dy*c.dz);
	grille[in1][jn1][kn1] = c;
      }
    } 
  } 
}	

/*!\brief Computation of the swept quantity for an interface triangle between times t-dt and t. Computation of the boundary flux.
 \details Algorithm:\n
 - Construction of the vector of the bounding boxes of the prisms with bases T3d_prev and T3d_n. Loop on the prisms obtained:
 - Look for the indices of the cells containing T3d_prev (due to the construction, the triangle is fully contained in one single cell) and T3d_n (idem).
 - If the prism is contained in one single cell, compute the prism volume using function \a volume_prisme(const Triangle_3&,const Triangle_3&), and the quantity swept by the interface triangle (\a Particule.triangles) is given by: \f$  volume\_prisme*U^n/volume\_cellule \f$. Otherwise,
 - List fluid cells intersecting the prism using function \a cells_intersection_face(int& ,int& ,int& ,int& ,int& ,int& , std::vector<Bbox>& s, std::vector<Cellule>& s).
 - Split the prism into tetrahedra: let  \f$ T1(A_1,B_1,C_1)\f$  and \f$ T2(A_2,B_2,C_2)\f$  the prism bases, and define points: \f$ A = \frac{1}{4}(B_1 + B_2 + C_1 +C_2) \f$ , \f$ B = \frac{1}{4}(A_1 + A_2 + C_1 +C_2) \f$ et \f$ C = \frac{1}{4}(A_1 + A_2 + B_1 + B_2 ) \f$. The prism is split into the prisms \f$ A_1,B_1,C_1 A_2,B_2,C_2 \f$ sont: \f$ A_1 A_2 C B \f$, \f$ B_1 B_2 A C \f$, \f$ C_1 C_2 B A \f$, \f$ A_1 C C_1 B \f$, \f$ B_1 A C_1 C \f$, \f$ A C B C_1 \f$, \f$ A B C C_2 \f$, \f$ A B_2 C_2 C \f$, \f$ A_1 B_1 C_1 C \f$, \f$ A_2 C_2 C B \f$, \f$ A_2 B_2 C C_2. \f$
 - Intersection of these tetrahedra with the fluid cells intersected by the prism using function \a intersect_cube_tetrahedron(Bbox&, Tetrahedron&). The quantity swept by the interface triangle is given by the sum of the following terms: \f$  volume\_{intersection\_cellule\_tetrahedre}*U^n/volume\_cellule. \f$ \n

 Computation of the boundary flux: let \a f an interface triangle, the boundary flux is given by:
 \f{eqnarray*}{
 \Phi_f  =  \left(0,  p^x \, A_f n^{x}_f, \, p^y \,A_f n^{y}_f, \, p^z \,A_f n^{z}_f, V_f \cdot \left( p^x \, A_f n^{x}_f,p^y \,A_f n^{y}_f,p^z \,A_f n^{z}_f \right)^t \right)^t
 \f} \n
 where \f$ A_f \f$ is the area of f,  \f$ n_f \f$ is the exterior normal to f, \f$ V_f \f$ is the velocity at the center of \a f computed using function \a vitesse_parois(Point_3& ) and \f$ p^x, p^y, p^z \f$ are the effective pressures in the directions x, y and z during the time-step (\a Cellule.pdtx, \a  Cellule.pdty et \a Cellule.pdtz).

 \param T3d_prev Triangles_3 (interface triangles at time t: \a Particule.Triangles_interface)
 \param T3d_n    Triangles_3 (interface triangles at time t-dt: \a Particule.Triangles_interface_prev)
 \param dt Time-step
 \param P Particule 
 \warning <b> Specific coupling procedure ! </b>
 \return void
 */
void Grille::swap_face_inexact(const Triangle_3& Tr_prev, const Triangle_3& Tr, const Triangles& T3d_prev, const Triangles& T3d_n, const double dt, Particule & P, double & volume_test){
  CGAL::Timer delta_time,total_time,boucle1_time,boucle2_time;
  delta_time.start();total_time.start();boucle1_time.start();boucle2_time.start();
  double temps_delta=0.,temps_total=0.,temps_boucle1=0.,temps_boucle2=0.,temps_intersections=0.,temps_triangulation=0.;
  delta_time.reset();
  Bbox box_prisme = Tr_prev.bbox()+Tr.bbox();

  //Definition tetrahedra 
  Tetrahedron Tet[11];
  Point_3 e = centroid(Tr_prev.operator[](1),Tr.operator[](1), Tr_prev.operator[](2),Tr.operator[](2));
  Point_3 f = centroid(Tr_prev.operator[](0),Tr.operator[](0), Tr_prev.operator[](2),Tr.operator[](2));
  Point_3 g = centroid(Tr_prev.operator[](0),Tr.operator[](0), Tr_prev.operator[](1),Tr.operator[](1));
  
  Tet[0] = Tetrahedron(Tr_prev.operator[](0),Tr.operator[](0), g, f);
  Tet[1] = Tetrahedron(Tr_prev.operator[](1),Tr.operator[](1), e, g);
  Tet[2] = Tetrahedron(Tr_prev.operator[](2),Tr.operator[](2), f, e);
  Tet[3] = Tetrahedron(Tr_prev.operator[](0),Tr_prev.operator[](1), Tr_prev.operator[](2), g);
  Tet[4] = Tetrahedron(Tr_prev.operator[](0), f, g, Tr_prev.operator[](2));
  Tet[5] = Tetrahedron(Tr_prev.operator[](1), g, e, Tr_prev.operator[](2));
  Tet[6] = Tetrahedron(e,g,f, Tr_prev.operator[](2));
  Tet[7] = Tetrahedron(e,f,g, Tr.operator[](2));
  Tet[8] = Tetrahedron(e, Tr.operator[](1),Tr.operator[](2), g);
  Tet[9] = Tetrahedron(Tr.operator[](0),Tr.operator[](2),g,f);
  Tet[10]= Tetrahedron(Tr.operator[](0),Tr.operator[](1), g, Tr.operator[](2));

  double test=volume_prisme(Tr_prev,Tr);
  double test2=0.;
  for(int l=0;l<11;l++){
    test2 += CGAL::to_double(Tet[l].volume());
  }
  
  double delta_w_tot[5];
  delta_w_tot[0] = delta_w_tot[1] = delta_w_tot[2] = delta_w_tot[3] = delta_w_tot[4] =0.;
  double volume_tot = 0.;
  //Computation of the quantity swept by the face
  for(int t=0;t<11;t++){
    double volume_tet = 0.;
    for(int i=0;i<Nx+2*marge;i++){
      for(int j=0;j<Ny+2*marge;j++){
	for(int k=0;k<Nz+2*marge;k++){
	  Cellule& c= grille[i][j][k];
	  Bbox box_cell(c.x -c.dx/2.,c.y -c.dy/2.,c.z -c.dz/2.,c.x +c.dx/2.,c.y +c.dy/2.,c.z + c.dz/2.);
	  
	  if (CGAL::do_overlap(box_prisme, box_cell) ) {
	    if(CGAL::do_overlap(Tet[t].bbox(), box_cell)){
	      double volume = (intersect_cube_tetrahedron(box_cell, Tet[t],temps_intersections,temps_triangulation) * sign(Tet[t].volume()) );
	      
	      volume_test += volume;
	      volume_tot += volume;
	      volume_tet += volume;
	      delta_w_tot[0] += volume*c.rho0; 
	      delta_w_tot[1] += volume*c.impx0;
	      delta_w_tot[2] += volume*c.impy0; 
	      delta_w_tot[3] += volume*c.impz0; 
	      delta_w_tot[4] += volume*c.rhoE0;
	    }
	  }
	} 
      } 
    }
  }
  temps_delta += delta_time.time();
  
  //Preliminary loop on the interface triangles
  //Rough evaluation of the swept quantity on each interface triangle
  boucle1_time.reset();
  double volume_eval = 0.;
  std::vector<Bbox> box_prismes(T3d_prev.size());
  for (int i=0; i< T3d_prev.size(); i++){
    const Bbox& box_triangles_prev = T3d_prev[i].bbox();
    const Bbox& box_triangles_n = T3d_n[i].bbox();
    box_prismes[i]= box_triangles_prev + box_triangles_n;
  } 
  for (int i=0; i< box_prismes.size(); i++){
    int in=0, jn=0, kn=0, in1=0, jn1=0, kn1=0;
    bool interieur = true;
    const Point_3& center_prev= centroid(T3d_prev[i].operator[](0),T3d_prev[i].operator[](1),T3d_prev[i].operator[](2));
    const Point_3& center_n= centroid(T3d_n[i].operator[](0),T3d_n[i].operator[](1),T3d_n[i].operator[](2));
    in_cell(center_prev, in, jn, kn, interieur);
    in_cell(center_n, in1, jn1, kn1, interieur);

    
    if((std::abs(grille[in1][jn1][kn1].alpha -1.)<eps)  && (interieur==true)){
      double x= CGAL::to_double(center_n.operator[](0));
      double y= CGAL::to_double(center_n.operator[](1));
      double z= CGAL::to_double(center_n.operator[](2));
      Cellule& cd= grille[in1+1][jn1][kn1];
      if (cd.is_in_cell(x,y,z) ) {in1=in1+1;}
      else{
	Cellule& cg= grille[in1-1][jn1][kn1];
	if (cg.is_in_cell(x,y,z)) {in1=in1-1;}
	else{
	  Cellule& ch= grille[in1][jn1+1][kn1];
	  if (ch.is_in_cell(x,y,z)) {jn1=jn1+1;}
	  else{
	    Cellule& cb= grille[in1][jn1-1][kn1];
	    if (cb.is_in_cell(x,y,z)) {jn1=jn1-1;}
	    else{
	      Cellule& cd= grille[in1][jn1][kn1+1];
	      if (cd.is_in_cell(x,y,z)) {kn1=kn1+1;}
	      else{
		Cellule& cder= grille[in1][jn1][kn1-1];
		if (cder.is_in_cell(x,y,z)) {kn1=kn1-1;}
	      }
	    }
	  }
	}
      }
    }
		  
    Cellule& c= grille[in1][jn1][kn1];
    Cellule& c_prev= grille[in][jn][kn];
    double volume_cel = c.dx*c.dy*c.dz;  
    if ((interieur==true)){ 
      double volume_p=volume_prisme(T3d_prev[i],T3d_n[i]);
      //Evaluation of the swept quantity as the product of the volume of the prism by the value of the fluid in the cell
      if( (std::abs(volume_p)>eps) && (std::abs(1.-c.alpha)>eps)){
	c.delta_w[0] += volume_p*c_prev.rho0/volume_cel; 
	c.delta_w[1] += volume_p*c_prev.impx0/volume_cel;
	c.delta_w[2] += volume_p*c_prev.impy0/volume_cel; 
	c.delta_w[3] += volume_p*c_prev.impz0/volume_cel; 
	c.delta_w[4] += volume_p*c_prev.rhoE0/volume_cel;
	volume_eval += abs(volume_p);
	delta_w_tot[0] -= volume_p*c_prev.rho0; 
	delta_w_tot[1] -= volume_p*c_prev.impx0;
	delta_w_tot[2] -= volume_p*c_prev.impy0; 
	delta_w_tot[3] -= volume_p*c_prev.impz0; 
	delta_w_tot[4] -= volume_p*c_prev.rhoE0;
      }
    }

    if (explicite){
      Vector_3 norm_prev= orthogonal_vector(T3d_prev[i].operator[](0),T3d_prev[i].operator[](1),T3d_prev[i].operator[](2));
      double norm2_prev= sqrt(CGAL::to_double(norm_prev*norm_prev));
      if(norm2_prev>eps){
	Vector_3 n_prev = norm_prev/norm2_prev;
	double aire_prev = sqrt(CGAL::to_double(T3d_prev[i].squared_area()));
	c.phi_x += c_prev.pdtx * aire_prev *( CGAL::to_double(n_prev.x()))/volume_cel;
	c.phi_y += c_prev.pdty * aire_prev *( CGAL::to_double(n_prev.y()))/volume_cel;
	c.phi_z += c_prev.pdtz * aire_prev *( CGAL::to_double(n_prev.z()))/volume_cel;
			
	Vector_3 V_f = P.vitesse_parois_prev(center_prev);
	c.phi_v += aire_prev * (CGAL::to_double(c.pdtx*n_prev.x()*V_f.x()  + c.pdty*n_prev.y()*V_f.y()+
						c.pdtz*n_prev.z()*V_f.z()))/volume_cel;
      }
    }
    else {
      Vector_3 norm= orthogonal_vector(T3d_n[i].operator[](0),T3d_n[i].operator[](1),T3d_n[i].operator[](2));
      double norm2= sqrt(CGAL::to_double(norm*norm));
      if(norm2>eps){ 
	Vector_3 n = norm/norm2;
	double aire = sqrt(CGAL::to_double(T3d_n[i].squared_area()));
	c.phi_x += c.pdtx * aire *( CGAL::to_double(n.x()))/(c.dx*c.dy*c.dz);
	c.phi_y += c.pdty * aire *( CGAL::to_double(n.y()))/(c.dx*c.dy*c.dz);
	c.phi_z += c.pdtz * aire *( CGAL::to_double(n.z()))/(c.dx*c.dy*c.dz);
	Vector_3 V_f = P.vitesse_parois(center_n);
	c.phi_v += aire * (CGAL::to_double(c.pdtx*n.x()*V_f.x()  + c.pdty*n.y()*V_f.y() + c.pdtz*n.z()*V_f.z()))/(c.dx*c.dy*c.dz);
	
      }
    } 
  }
  temps_boucle1 += boucle1_time.time();

  //Second loop on interface triangles: smoothes the error on swept quantity among cells according to their volume
  boucle2_time.reset();
  for (int i=0; i< box_prismes.size(); i++){
    int in=0, jn=0, kn=0, in1=0, jn1=0, kn1=0;
    bool interieur = true;
    const Point_3& center_prev= centroid(T3d_prev[i].operator[](0),T3d_prev[i].operator[](1),T3d_prev[i].operator[](2));
    const Point_3& center_n= centroid(T3d_n[i].operator[](0),T3d_n[i].operator[](1),T3d_n[i].operator[](2));
    in_cell(center_prev, in, jn, kn, interieur);
    in_cell(center_n, in1, jn1, kn1, interieur);
		
    if((std::abs(grille[in1][jn1][kn1].alpha -1.)<eps)  && (interieur==true)){
      double x= CGAL::to_double(center_n.operator[](0));
      double y= CGAL::to_double(center_n.operator[](1));
      double z= CGAL::to_double(center_n.operator[](2));
      Cellule& cd= grille[in1+1][jn1][kn1];
      if (cd.is_in_cell(x,y,z) ) {in1=in1+1;}
      else{
	Cellule& cg= grille[in1-1][jn1][kn1];
	if (cg.is_in_cell(x,y,z)) {in1=in1-1;}
	else{
	  Cellule& ch= grille[in1][jn1+1][kn1];
	  if (ch.is_in_cell(x,y,z)) {jn1=jn1+1;}
	  else{
	    Cellule& cb= grille[in1][jn1-1][kn1];
	    if (cb.is_in_cell(x,y,z)) {jn1=jn1-1;}
	    else{
	      Cellule& cd= grille[in1][jn1][kn1+1];
	      if (cd.is_in_cell(x,y,z)) {kn1=kn1+1;}
	      else{
		Cellule& cder= grille[in1][jn1][kn1-1];
		if (cder.is_in_cell(x,y,z)) {kn1=kn1-1;}
	      }
	    }
	  }
	}
      }
    }
		  
    Cellule& c= grille[in1][jn1][kn1];
    double volume_cel = c.dx*c.dy*c.dz;  
    if ((interieur==true)){ 
      double volume_p=volume_prisme(T3d_prev[i],T3d_n[i]);
      //Evaluation of the swept quantity as the product of the volume of the prism by the value in the cell
      if( (std::abs(volume_p)>eps) && (std::abs(1.-c.alpha)>eps)){
	for(int l=0;l<5;l++){
	  c.delta_w[l] += abs(volume_p)/volume_eval*delta_w_tot[l]/volume_cel;
	}
      }
    }
  }
  temps_boucle2 += boucle2_time.time();
}	

/*!\brief Construction of the 2d submesh of two meshes on a triangular 2d face of the solid.
   \details The goal is to split the solid face into interface triangles such that each triangle is fully enclosed in one single cell at times t-dt and t (npt necessarily the same cell).\n
   Algorithm:\n
   - Construct the 2d bounding box for each triangle.  \n
   - Loop on the bounding boxes.
   - Test the intersection of the bounding boxes using function <b> CGAL::do_overlap(Bbox_2, Bbox_2)</b>. If the boxes intersect: \n
   - Test the intersection of the triangles using function <b> CGAL::do_intersect(Triangle_2,Triangle_2)</b>. If the triangles intersect:  \n
   - Compute the intersections between the two triangles using function <b>CGAL::intersection(Triangle_2,Triangle_2)</b>. \n
   - If the intersection is a triangle, add it to the submesh. If the intersection is a polygon, triangulate it using class <b>CGAL::Triangulation</b> and function <b> CGAL::insert </b> in the class. If the intersection is a point or a segment, do nothing since the swept volume is null.
   \warning <b> Specific coupling procedure ! </b>
   \param Tn Triangles_2 
   \param Tn1 Triangles_2
   \param tri2 Triangles_2
   \return void
*/
void Sous_Maillage_2d(const Triangles_2& Tn, const Triangles_2& Tn1, Triangles_2& tri2){
	
  Triangle_2 tri;
  std::vector<Point_2> vPoints; 
  int k=0;
  CGAL::Timer user_time;
  user_time.start();
  for(Triangles_2::const_iterator t=Tn.begin(); t!=Tn.end(); t++){ 
    for(Triangles_2::const_iterator t1=Tn1.begin(); t1!=Tn1.end(); t1++){
      if (CGAL::do_overlap((*t).bbox(),(*t1).bbox())){ //test the intersection of Bbox 
	if (CGAL::do_intersect(*t,*t1)){ // test the intersection of the triangles
	  const CGAL::Object& result = CGAL::intersection(*t,*t1); //Compute the intersection between two triangles
	  if(CGAL::assign(tri,result)){ tri2.push_back(tri); }
	  else if(CGAL::assign(vPoints,result)){
	    Triangulation_2 T;
	    T.insert(vPoints.begin(), vPoints.end());
	    for (Triangulation_2::Finite_faces_iterator fit=T.finite_faces_begin(); fit!=T.finite_faces_end();++fit)
	    { 
	      if(T.triangle(fit).area()>eps){
		tri2.push_back(T.triangle(fit));
	      }
	    }
	  }
	}
      }
    }
  }
}

/*!\brief Split the face into a submesh of the intersections of the solid face with the fluid grid at times t-dt and t
   \details From the position of the interface at time t (\a Tn) and t-dt (\a Tn1), split the face into triangles fully contained in one single cell at times t and t-dt (not necessarily the same cell). \n
   Algorithm:\n
   - Barycentric transformation of \a Tn using function \a tr(Triangle_3, Triangle_3, Triangle_3) and transformation of the resulting triangles into 2d triangles using function \a tr(Triangle_3, Triangle_3).
   - Transformation of tn1 into 2d triangles 2d using function \a tr(Triangle_3, Triangle_3).
   - Construction of the 2d submesh of the face using function \a Sous_Maillage_2d(const Triangles_2&, const Triangles_2&, Triangles_2&).
   - Transformation of the 2d submesh 2d into a 3d submesh of the face using function \a tr(Triangle_3, Triangle_2).

   \warning <b> Specific coupling procedure !</b> 
   \param Tn const Triangle_3: Interface triangle at time t-dt (\a Particule.triangles_prev)
   \param tn const vector<Triangle_3>: Triangulation of face Tn (\a Particule.Triangles_interface_prev)
   \param Tn1 const Triangle_3: Interface at time t (\a Particule.triangles)
   \param tn1 const vector<Triangle_3>: Triangulation of face Tn1 (\a Particule.Triangles_interface)
   \param N   const Vector_3: exterior normal vector to Tn1 (\a Particule.normales)
   \param T3d_n vector of Triangle_3: Triangular submesh of face \a Particule.Triangles_interface_prev at time t-dt 
   \param T3d_n1 vector of Triangle_3: Triangular submesh of face \a Particule.Triangles_interface at time t 
   \return void
*/
void sous_maillage_faceTn_faceTn1_2d(const Triangle_3& Tn, const Triangles& tn, const Triangle_3& Tn1, const Triangles& tn1, const Vector_3& N,Triangles& T3d_n,Triangles& T3d_n1){
	
  CGAL::Timer total_time,bary1_time,bary2_time,bary3_time,sous_maillage_time,time_2d_3d,time_2d_3d_bis;
  double temps_total=0.,temps_bary1=0.,temps_bary2=0.,temps_bary3=0.,temps_sous_maillage=0.,temps_2d_3d=0.,temps_2d_3d_bis=0.;	
  total_time.start();bary1_time.start();bary2_time.start();bary3_time.start();sous_maillage_time.start();time_2d_3d.start();time_2d_3d_bis.start();
  
  bary2_time.reset();
  Triangles_2 Tn_2(tn.size());
  for(int i=0; i<tn.size(); i++){
    Tn_2[i] = tr(Tn, tn[i]);
  }
  temps_bary2 += bary2_time.time();
  
  bary3_time.reset();
  Triangles_2 Tn1_2(tn1.size());
  for(int i=0; i<tn1.size(); i++){
    Tn1_2[i] =tr(Tn1, tn1[i]);
  }
  temps_bary3 += bary3_time.time();
  
  sous_maillage_time.reset();
  Triangles_2 tri2;
  Sous_Maillage_2d(Tn1_2, Tn_2, tri2);
  temps_sous_maillage += sous_maillage_time.time();
  
  time_2d_3d.reset();
  T3d_n1.clear();
  for(Triangles_2::const_iterator t2=tri2.begin();t2!=tri2.end();t2++){
    const Triangle_3& Tri = tr(Tn1,*t2);
    Vector_3 vect0(Tri.operator[](0),Tri.operator[](1));
    Vector_3 vect1(Tri.operator[](0),Tri.operator[](2));
    Vector_3 normale = CGAL::cross_product(vect0,vect1);
    if (normale*N > 0.){ T3d_n1.push_back(Tri); }
    else{ T3d_n1.push_back(Triangle_3(Tri.operator[](0),Tri.operator[](2),Tri.operator[](1)));}
  }
  temps_2d_3d += time_2d_3d.time();
  time_2d_3d_bis.reset();
  T3d_n.clear();
  for(Triangles_2::const_iterator t2=tri2.begin();t2!=tri2.end();t2++){
    const Triangle_3& Tri = tr(Tn,*t2);
    Vector_3 vect0(Tri.operator[](0),Tri.operator[](1));
    Vector_3 vect1(Tri.operator[](0),Tri.operator[](2));
    Vector_3 normale = CGAL::cross_product(vect0,vect1);
    if (normale*N > 0.){ T3d_n.push_back(Tri); }
    else{ T3d_n.push_back(Triangle_3(Tri.operator[](0),Tri.operator[](2),Tri.operator[](1)));}
  }
  temps_2d_3d_bis += time_2d_3d_bis.time();
}

/*!\brief Computation of the quantity of fluid swept by the solid between times t-dt and t.
   \details Algorithm:\n
   - Split the solid faces (\a Particule.triangles and \a Particule.triangles_prev) into triangles fully contained in one single cell at times t-dt and t (not necessarily the same cell) using function \a sous_maillage_faceTn_faceTn1_2d(Triangle_3&, Triangles&, Triangle_3&, Triangles&, Vector_3& ,Triangles& ,Triangles&).\n
   - Compute the swept quantity and the boundary flux using function \a swap_face(Triangles&, Triangles&, const double ,  Particule &).
   \warning <b> Specific coupling procedure ! </b> 
   \param S Solide
   \param dt Time-step
   \return void
*/

void Grille::Swap_2d(const double dt, Solide& S){
	
  CGAL::Timer swap_face_time,total_time,sous_maillage_time;
  swap_face_time.start();total_time.start();sous_maillage_time.start();
  double temps_swap_face=0.,temps_total=0.,temps_sous_maillage=0.,nb=0.;
  double volume_test=0.;
  for(int i=0;i<S.solide.size();i++){
    for (int j=0; j<S.solide[i].triangles.size(); j++){
      if (S.solide[i].fluide[j]){
	nb+=1.;
	Triangles T3d_n,T3d_n1;
	sous_maillage_time.reset();
	sous_maillage_faceTn_faceTn1_2d(S.solide[i].triangles_prev[j], S.solide[i].Triangles_interface_prev[j], S.solide[i].triangles[j], S.solide[i].Triangles_interface[j], S.solide[i].normales[j], T3d_n, T3d_n1);
	temps_sous_maillage += sous_maillage_time.time();
	swap_face_time.reset();
	if(exact_swap){
	  //Exact swept quantity
	  swap_face(T3d_n,T3d_n1,dt, S.solide[i],volume_test );
	} else {
	  //Inexact swept quantity: compute the swept quatity inexactly and distribute the default of fluid on the interface elements
	  swap_face_inexact(S.solide[i].triangles_prev[j],S.solide[i].triangles[j],T3d_n,T3d_n1,dt,S.solide[i],volume_test);
	}
	temps_swap_face += swap_face_time.time();
      }
    }
  }
  temps_total += total_time.time();
  cout << "################# COUT SWAP ##############" << endl;
  cout << "sous_maillage=" << 100*temps_sous_maillage/temps_total << "%     t_moy=" << temps_sous_maillage/nb << " nb=" << nb << endl;
  cout << "swap_face=" << 100*temps_swap_face/temps_total << "%     t_moy=" << temps_swap_face/nb << " nb=" << nb << endl;
  cout << "Reste=" << 100*(1.-(temps_sous_maillage+temps_swap_face)/temps_total) << "%" << endl;
  cout << "##########################################" << endl;
  cout<<"volume balayee = "<< volume_test<<endl;
}


/*!\brief Conservative mixing of small cut-cells.
   \details We define a small cut-cell as a cell such that \f$ alpha > epsa \f$ (\a Cellule.alpha si the solid occupancy ratio in the cell, and \a epsa is the maximum value accepted for CFL stability reasons and is set in \a parametres.hpp ). In order not to modify the time-step while keeping the original CFL condition, the small cut-cells are merged with their neighbours.
   \warning <b> Specific coupling procedure ! </b>
   \return void
 */
void Grille::Mixage_cible(){
	
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){ 
      for(int k=marge;k<Nz+marge;k++){
				
	grille[i][j][k].cible_alpha = 0.;
	grille[i][j][k].cible_rho = 0.;
	grille[i][j][k].cible_impx = 0.;
	grille[i][j][k].cible_impy = 0.;
	grille[i][j][k].cible_impz = 0.;
	grille[i][j][k].cible_rhoE = 0.;
      }
    }
  }
	
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){ 
      for(int k=marge;k<Nz+marge;k++){
	Cellule& cp = grille[i][j][k];
	if((cp.alpha>epsa || cp.p<0. || cp.rho<0.) && abs(cp.alpha-1.)>eps && !cp.vide){
	  std::vector< std::vector<int> > tab_cible; 
	  std::vector<int> poz(3); poz[0]= i; poz[1] = j; poz[2] = k; tab_cible.push_back(poz);


          Cellule cg = cible(grille[i][j][k], tab_cible);

	  cg.cible_alpha += (1.-cp.alpha);
	  cg.cible_rho  += (1.-cp.alpha)*cp.rho;
	  cg.cible_impx += (1.-cp.alpha)*cp.impx;
	  cg.cible_impy += (1.-cp.alpha)*cp.impy;
	  cg.cible_impz += (1.-cp.alpha)*cp.impz;
	  cg.cible_rhoE += (1.-cp.alpha)*cp.rhoE;
					
	  cp.cible_i= cg.i;
	  cp.cible_j = cg.j;
	  cp.cible_k = cg.k; 
	
	}
	else{
	  grille[i][j][k].cible_i = i;
	  grille[i][j][k].cible_j = j;
	  grille[i][j][k].cible_k = k;
	}
      }
    }
  }
	
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){ 
      for(int k=marge;k<Nz+marge;k++){
	Cellule& cp = grille[i][j][k];

	if(std::abs(cp.cible_alpha)>0. && !cp.vide){
	  cp.rho = ((1.-cp.alpha)*cp.rho + cp.cible_rho)/((1.-cp.alpha) + cp.cible_alpha);
	  cp.impx = ((1.-cp.alpha)*cp.impx + cp.cible_impx)/((1.-cp.alpha) + cp.cible_alpha);
	  cp.impy = ((1.-cp.alpha)*cp.impy + cp.cible_impy)/((1.-cp.alpha) + cp.cible_alpha);
	  cp.impz = ((1.-cp.alpha)*cp.impz + cp.cible_impz)/((1.-cp.alpha) + cp.cible_alpha);
	  cp.rhoE = ((1.-cp.alpha)*cp.rhoE + cp.cible_rhoE)/((1.-cp.alpha) + cp.cible_alpha);
	  if(std::abs(cp.rho) > eps_vide){
	    cp.u = cp.impx/cp.rho;
	    cp.v = cp.impy/cp.rho;
	    cp.w = cp.impz/cp.rho;
	    cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
	    if(std::abs(cp.p) < eps_vide){
	      cp.vide=true;
	    }
	  }
	  else{
	    cp.u = 0.;
	    cp.v = 0.;
	    cp.w = 0.;
	    cp.p = 0.;
	    cp.vide=true;
	  }
	}
      }
    }
  }
  bool test_fini = true;
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){ 
      for(int k=marge;k<Nz+marge;k++){
	Cellule& cp = grille[i][j][k];
	Cellule& cible=grille[cp.cible_i][cp.cible_j][cp.cible_k];
	cp.rho = cible.rho;
	cp.impx = cible.impx;
	cp.impy = cible.impy;
	cp.impz = cible.impz;
	cp.rhoE = cible.rhoE;
	if(std::abs(cp.rho) > eps_vide){
	  cp.u = cp.impx/cp.rho;
	  cp.v = cp.impy/cp.rho;
	  cp.w = cp.impz/cp.rho;
	  cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
	  if(std::abs(cp.p) < eps_vide){
	    cp.vide=true;
	  }
	}
	else{
	  cp.u = 0.;
	  cp.v = 0.;
	  cp.w = 0.;
	  cp.p = 0.;
	  cp.vide=true;
	}
					
	if((grille[i][j][k].p<0. || grille[i][j][k].rho<0.) && !cp.vide){
	  cout << "Unfinished test x=" << grille[i][j][k].x << " y=" <<  grille[i][j][k].y << " z=" <<  grille[i][j][k].z << " p=" <<  grille[i][j][k].p << " rho=" <<  grille[i][j][k].rho << " target x=" << cible.x << " y=" << cible.y << " z=" << cible.z << endl;
	  test_fini = false;
	}
      }
    }
  }
	
  if(!test_fini){
    cout<<" Mixing is not completed "<<endl;
    Mixage_cible();
  }
}

/*!\brief Conservative mixing of small cut-cells.
   \details We define a small cut-cell as a cell such that \f$ alpha > epsa \f$ (\a Cellule.alpha si the solid occupancy ratio in the cell, and \a epsa is the maximum value accepted for CFL stability reasons and is set in \a parametres.hpp ). In order not to modify the time-step while keeping the original CFL condition, the small cut-cells are merged with their neighbours. This is an alternative version of \a Grille::Mixage_cible()
   \warning <b> Specific coupling procedure ! </b>
   \return void
 */
bool Grille::Mixage_cible2(){
  //Test whether there are still negative densities or pressure after mixing
  bool test_fini = true;
  cout << "Mixage_cible2" << endl;
  //Step 0: initialize
  for(int i=0;i<Nx+2*marge;i++){
    for(int j=0;j<Ny+2*marge;j++){ 
      for(int k=0;k<Nz+2*marge;k++){
	grille[i][j][k].cible_alpha = 0.;
	grille[i][j][k].cible_rho = 0.;
	grille[i][j][k].cible_impx = 0.;
	grille[i][j][k].cible_impy = 0.;
	grille[i][j][k].cible_impz = 0.;
	grille[i][j][k].cible_rhoE = 0.;
	grille[i][j][k].cible_i = i;
	grille[i][j][k].cible_j = j;
	grille[i][j][k].cible_k = k;
      }
    }
  }
  //Step 1: Define the target in a neighbourhood of each cell (keep the cell itself as target if it has no issue). We use to that end functions \a voisin_fluide (fully fluid neighbour, p>0, rho>0, minimal kappa), voisin_mixt (alpha_cible<alpha, p>0, rho>0, minimal kappa) and voisin (p and rho >0 if possible, minimal kappa)
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      for(int k=marge;k<Nz+marge;k++){
	Cellule& cp = grille[i][j][k];
	if((cp.alpha>epsa || cp.p<0. || cp.rho<0.) && abs(cp.alpha-1.)>eps && !cp.vide){
	  //Search for a candidate target cell with voisin_fluide
	  bool target = true;
	  Cellule cell_cible;
	  cell_cible = voisin_fluide(cp, target);
	  if(target){ 
	    cp.cible_i = cell_cible.i;
	    cp.cible_j = cell_cible.j;
	    cp.cible_k = cell_cible.k;
	  } else {
	    //Search for a candidate target cell with voisin_mixt
	    target = true;
	    cell_cible= voisin_mixt(cp,target);
	    if(target){
	      cp.cible_i = cell_cible.i;
	      cp.cible_j = cell_cible.j;
	      cp.cible_k = cell_cible.k;
	    } else {
	      //Search for a candidate target cell with voisin
	      cell_cible= voisin(cp);
	      cp.cible_i = cell_cible.i;
	      cp.cible_j = cell_cible.j;
	      cp.cible_k = cell_cible.k;
	    }
	  }
	}
      }
    }
  }
  //Step 2: Determine the end target of each cell by following recursively the target of the target. Of course, stop as soon as there is no movement or a cycle. Then, update the equivalent quantities for the target cell after mixing
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      for(int k=marge;k<Nz+marge;k++){
	Cellule& cp = grille[i][j][k];
	if((cp.alpha>epsa || cp.p<0. || cp.rho<0.) && abs(cp.alpha-1.)>eps && !cp.vide){
	  //List of traveled points
	  std::vector< std::vector<int> > tab_cible; 
	  std::vector<int> poz(3); poz[0]= i; poz[1] = j; poz[2] = k; tab_cible.push_back(poz);
	  int l=cp.cible_i;
	  int m=cp.cible_j;
	  int n=cp.cible_k;
	  bool test=true;
	  for(int count=1;test;count++){
	    Cellule& cible = grille[l][m][n];
	    poz[0] = l;
	    poz[1] = m;
	    poz[2] = n;
	    //Test to check whether the target points to an already visited cell
	    for(int iter=0;iter<tab_cible.size() && test;iter++){
	      //If a cell has already been visited, put the test to false
	      if(abs(tab_cible[iter][0]-l)+abs(tab_cible[iter][1]-m)+abs(tab_cible[iter][2]-n)<eps){
		test = false;
	      }
	    }
	    //If the cell has never been visited, take it as a target cell
	    if(test){
	      tab_cible.push_back(poz);
	      l = cible.cible_i;
	      m = cible.cible_j;
	      n = cible.cible_k;
	    }
	    //Otherwise, take it as the target cell for all previously visited cells in the loop
	    else {
	      for(int iter=0;iter<tab_cible.size() && test;iter++){
		grille[tab_cible[iter][0]][tab_cible[iter][1]][tab_cible[iter][2]].cible_i = l;
		grille[tab_cible[iter][0]][tab_cible[iter][1]][tab_cible[iter][2]].cible_j = m;
		grille[tab_cible[iter][0]][tab_cible[iter][1]][tab_cible[iter][2]].cible_k = n;
	      }
	      grille[l][m][n].cible_i = l;
	      grille[l][m][n].cible_j = m;
	      grille[l][m][n].cible_k = n;
	    }
	  }
	  
	  cp.cible_i = l;
	  cp.cible_j = m;
	  cp.cible_k = n;
	  
	}
      }
    }
  }
  //Step 3: Update the values of the target cells (which point to themselves)
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      for(int k=marge;k<Nz+marge;k++){
	Cellule& cp =grille[i][j][k];
	Cellule& cg = grille[cp.cible_i][cp.cible_j][cp.cible_k];
	cg.cible_alpha += (1.-cp.alpha);
	cg.cible_rho  += (1.-cp.alpha)*cp.rho;
	cg.cible_impx += (1.-cp.alpha)*cp.impx;
	cg.cible_impy += (1.-cp.alpha)*cp.impy;
	cg.cible_impz += (1.-cp.alpha)*cp.impz;
	cg.cible_rhoE += (1.-cp.alpha)*cp.rhoE;      }
    }
  }
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      for(int k=marge;k<Nz+marge;k++){
	Cellule& cp =grille[i][j][k];
	if(abs(cp.cible_i-i)+abs(cp.cible_j-j)+abs(cp.cible_k-k)<eps && abs(cp.cible_alpha)>eps && !cp.vide){
	  cp.rho = cp.cible_rho/cp.cible_alpha;
	  cp.impx = cp.cible_impx/cp.cible_alpha;
	  cp.impy = cp.cible_impy/cp.cible_alpha;
	  cp.impz = cp.cible_impz/cp.cible_alpha;
	  cp.rhoE = cp.cible_rhoE/cp.cible_alpha;
	  if(std::abs(cp.rho) > eps_vide){
	    cp.u = cp.impx/cp.rho;
	    cp.v = cp.impy/cp.rho;
	    cp.w = cp.impz/cp.rho;
	    cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
	    if(std::abs(cp.p) < eps_vide){
	      cp.vide=true;
	    }
	  }
	  else{
	    cp.u = 0.;
	    cp.v = 0.;
	    cp.w = 0.;
	    cp.p = 0.;
	    cp.vide=true;
	  }
	  if(cp.rho<0. || cp.p<0.){
	    Cellule& cible = grille[cp.cible_i][cp.cible_j][cp.cible_k];
	    cout << "Unfinished test x=" << cp.x << " y=" <<  cp.y << " z=" <<  cp.z << " p=" <<  cp.p << " rho=" <<  cp.rho << " alpha=" << cp.alpha << " target x=" << cible.x << " y=" << cible.y << " z=" << cible.z << " p=" << cible.p << " rho=" << cible.rho << " alpha=" << cible.alpha << " cible_alpha=" << cp.cible_alpha << " cible_rho=" << cp.cible_rho << " cible_rhoE" << cp.cible_rhoE << endl;
	    test_fini = false;
	    //Search for a possible target cell with voisin_fluide
	    bool target = true;
	    Cellule cell_cible;
	    cell_cible = voisin_fluide(cp, target);
	    if(target){ 
	      cout << "voisin_fluide x=" << cell_cible.x << " y=" << cell_cible.y << " z=" << cell_cible.z << " rho=" << cell_cible.rho << " p=" << cell_cible.p << " alpha=" << cell_cible.alpha << endl;
	    } else {
	      //Search for a possible target cell with voisin_mixt
	      target = true;
	      cell_cible= voisin_mixt(cp,target);
	      if(target){
		cout << "voisin_mixt x=" << cell_cible.x << " y=" << cell_cible.y << " z=" << cell_cible.z << " rho=" << cell_cible.rho << " p=" << cell_cible.p << " alpha=" << cell_cible.alpha << endl;
	      } else {
		//Search for a possible target cell with voisin
		cell_cible= voisin(cp);
		cout << "voisin x=" << cell_cible.x << " y=" << cell_cible.y << " z=" << cell_cible.z << " rho=" << cell_cible.rho << " p=" << cell_cible.p << " alpha=" << cell_cible.alpha << endl;
	      }
	    }
	  }
	}
      }
    }
  }
  //Step 4: put the value of the target in each cell
  for(int i=marge;i<Nx+marge;i++){
    for(int j=marge;j<Ny+marge;j++){
      for(int k=marge;k<Nz+marge;k++){
	Cellule& cp = grille[i][j][k];
	if(abs(cp.cible_i-i)+abs(cp.cible_j-j)+abs(cp.cible_k-k)>eps){
	  Cellule& cg = grille[cp.cible_i][cp.cible_j][cp.cible_k];
	  cp.rho = cg.rho;
	  cp.impx = cg.impx;
	  cp.impy = cg.impy;
	  cp.impz = cg.impz;
	  cp.rhoE = cg.rhoE;
	  cp.u = cg.u;
	  cp.v = cg.v;
	  cp.w = cg.w;
	  cp.p = cg.p;
	  cp.vide = cg.vide;
	}
      }
    }
  }
  return test_fini;
}
