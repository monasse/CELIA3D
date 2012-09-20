#include "fluide.hpp"
#include "intersections.cpp"


void Grille::modif_fnum(const double dt){
	
	Cellule c,ci,cj,ck; 
	for(int i=marge;i<Nx+marge;i++){
		for(int j=marge;j<Ny+marge;j++){ 
			for(int k=marge;k<Nz+marge;k++){
				c = grille[i][j][k];
				if(std::abs(c.alpha-1.)>eps){
					ci = grille[i-1][j][k];    //Cellule  i-1
					cj = grille[i][j-1][k];    //Cellule  j-1
					ck = grille[i][j][k-1];    //Cellule  k-1
           
          c.flux_modif[0] = 0.;
					c.flux_modif[1] = c.pdtx * c.phi_x;
					c.flux_modif[2] = c.pdty * c.phi_y;
					c.flux_modif[3] = c.pdtz * c.phi_z;
					c.flux_modif[4] = 0.;

					
					for(int l=0.; l<5; l++){  
						c.flux_modif[l] -= (1.-c.kappai)*c.dtfxi[l] - (1.-ci.kappai)*ci.dtfxi[l]
						+ (1.-c.kappaj)*c.dtfyj[l] - (1.-cj.kappaj)*cj.dtfyj[l]
						+ (1.-c.kappak)*c.dtfzk[l] - (1.-ck.kappak)*ck.dtfzk[l];
						if(abs(c.flux_modif[l])>eps){
						c.flux_modif[l] /= (1.-c.alpha);
						}
						else {c.flux_modif[l] = 0.;}
					}		
					//Mise a jour des valeurs dans le cellules
					c.rho = c.rho0  +  c.flux_modif[0];
					c.impx = c.impx0 + c.flux_modif[1];
					c.impy = c.impy0 + c.flux_modif[2];
					c.impz = c.impz0 + c.flux_modif[3];
					c.rhoE = c.rhoE0 + c.flux_modif[4];
					c.u = c.impx/c.rho;
					c.v = c.impy/c.rho;
					c.w = c.impz/c.rho;
					c.p = (gam-1.)*(c.rhoE-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2. - c.rho*c.w*c.w/2.);
					
				}
				
				else {
					c.rho = 0.;
					c.impx = 0.;
					c.impy = 0.;
					c.impz = 0.;
					c.rhoE = 0.;
					c.u = 0.;
					c.v = 0.;
					c.w = 0.;
					c.p = 0.;;
				}
				grille[i][j][k] = c;      
			}
		}
	}
}


void Grille:: mixage(){
	
	Cellule cp, cg;
	bool test=true;
	for(int i=marge;i<Nx+marge;i++){
		for(int j=marge;j<Ny+marge;j++){ 
			for(int k=marge;k<Nz+marge;k++){
				cp = grille[i][j][k];
				if((cp.alpha>epsa || cp.p<0. || cp.rho<0.) && abs(cp.alpha-1.)>eps){
					
					if (grille[i-1][j][k].alpha == 0.)
					{
						cg = grille[i-1][j][k];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						 
						grille[i][j][k] = cp;
						grille[i-1][j][k] = cg;
					}
					else if (grille[i+1][j][k].alpha == 0.)
					{
						cg = grille[i+1][j][k];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						
						grille[i][j][k] = cp;
						grille[i+1][j][k] = cg;
					}
					
					else if (grille[i][j-1][k].alpha == 0.)
					{
						cg = grille[i][j-1][k];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						
						grille[i][j][k] = cp;
						grille[i][j-1][k] = cg;
					}
					else if (grille[i][j+1][k].alpha == 0.)
					{
						cg = grille[i][j+1][k];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2.- cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						grille[i][j][k] = cp;
						grille[i][j+1][k] = cg;
					}
					else if (grille[i][j][k-1].alpha == 0.)
					{
						cg = grille[i][j][k-1];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						grille[i][j][k] = cp;
						grille[i][j][k-1] = cg;
					}
					else if(grille[i][j][k+1].alpha == 0.)
					{
						cg = grille[i][j][k+1];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						grille[i][j][k] = cp;
						grille[i][j][k+1] = cg;
					}
					
					else if (grille[i-2][j][k].alpha == 0.)
					{
						cg = grille[i-2][j][k];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						
						grille[i][j][k] = cp;
						grille[i-2][j][k] = cg;
					}
					else if (grille[i+2][j][k].alpha == 0.)
					{
						cg = grille[i+2][j][k];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						
						grille[i][j][k] = cp;
						grille[i+2][j][k] = cg;
					}
					
					else if (grille[i][j-2][k].alpha == 0.)
					{
						cg = grille[i][j-2][k];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						
						grille[i][j][k] = cp;
						grille[i][j-2][k] = cg;
					}
					else if (grille[i][j+2][k].alpha == 0.)
					{
						cg = grille[i][j+2][k];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2.- cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						grille[i][j][k] = cp;
						grille[i][j+2][k] = cg;
					}
					else if (grille[i][j][k-2].alpha == 0.)
					{
						cg = grille[i][j][k-2];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						grille[i][j][k] = cp;
						grille[i][j][k-2] = cg;
					}
					else if(grille[i][j][k+2].alpha == 0.)
					{
						cg = grille[i][j][k+2];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						grille[i][j][k] = cp;
						grille[i][j][k+2] = cg;
					}
					else if(grille[i-1][j-1][k].alpha <= eps_relat)
					{
						cg = grille[i-1][j-1][k];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						grille[i][j][k] = cp;
						grille[i-1][j-1][k] = cg;
					}
					else if(grille[i-1][j+1][k].alpha <= eps_relat)
					{
						cg = grille[i-1][j+1][k];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						grille[i][j][k] = cp;
						grille[i-1][j+1][k] = cg;
					}
					else if(grille[i+1][j-1][k].alpha <= eps_relat)
					{
						cg = grille[i+1][j-1][k];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						grille[i][j][k] = cp;
						grille[i+1][j-1][k] = cg;
					}
					else if(grille[i+1][j+1][k].alpha <= eps_relat)
					{
						cg = grille[i+1][j+1][k];
						
						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
						
						cp.rho += cp.Mrho;
						cp.impx += cp.Mimpx;
						cp.impy += cp.Mimpy;
						cp.impz += cp.Mimpz;
						cp.rhoE += cp.MrhoE;
						cp.u = cp.impx/cp.rho;
						cp.v = cp.impy/cp.rho;
						cp.w = cp.impz/cp.rho;
						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
						
						
						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
						
						cg.rho += cg.Mrho;
						cg.impx += cg.Mimpx;
						cg.impy += cg.Mimpy;
						cg.impz += cg.Mimpz;
						cg.rhoE += cg.MrhoE;
						cg.u = cg.impx/cg.rho;
						cg.v = cg.impy/cg.rho;
						cg.w = cg.impz/cg.rho;
						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
						grille[i][j][k] = cp;
						grille[i+1][j+1][k] = cg;
					}
					else {
						if(test){std::cout<<"Pas de cellule cible pour le mixage"<<std::endl; 
										std::cout<<" grille qui faut melange "<<std::endl;
										cp.Affiche();
										test=false;
						}
					}
					
				}
				
			}
		}
	}
	
} 
void Grille::fill_cel(Solide& S){
	
	Cellule c, cm;
	int nb_part = S.size();
	double dist[6*nb_part];
	double dist_min = 100;
	int poz=0;
	double x_min=0., y_min=0., z_min = 0., x_max = 0., y_max=0., z_max=0.;

	//std::cout<<"center faces number: " <<count<<std::endl;
	for(int i=marge;i<Nx+marge;i++){
		for(int j=marge;j<Ny+marge;j++){ 
			for(int k=marge;k<Nz+marge;k++){
				c = grille[i][j][k];
				if((std::abs(c.alpha-1.)<eps)){
				  Point_3 center_cell(c.x, c.y, c.z);
				  int nbx=0, nby=0,nbz=0;
				  Point_3 projete(0.,0.,0.); //Projet� sur la face la plus proche
				  double dist_min = 10000000.;
				  for(int iter=0; iter<nb_part; iter++){
					Particule P = S.solide[iter];
					for(int it=0;it<P.size();it++){
					  Face F = P.faces[it];
					  if(F.voisin==-1){
						Plane_3 P(F.vertex[0].pos,F.vertex[1].pos,F.vertex[2].pos);
						for(int k=3;k<F.size() && P.is_degenerate();k++){//Test si le plan est degenere
						  P = Plane_3(F.vertex[0].pos,F.vertex[1].pos,F.vertex[k].pos);
						}
						Point_3 xP = P.projection(center_cell);
						//Test pour savoir si le projete est dans la face
						bool test = true;
						for(int k=0;k<F.size()-1 && test;k++){
						  Point_3 x1 = F.vertex[k].pos;
						  Point_3 x2 = F.vertex[k+1].pos;
						  Vector_3 vect1(xP,x1);
						  Vector_3 vect2(xP,x2);
						  if(CGAL::to_double(CGAL::cross_product(vect1,vect2)*F.normale)<0.){
							test = false;
						  }
						}
						//1er cas : on est dans la face
						if(test){
						  double d = sqrt(CGAL::to_double(CGAL::squared_distance(center_cell,xP)));
						  if(d<dist_min){
							dist_min = d;
							projete = xP;
						  }
						}
						//2eme cas : on est hors de la face
						else{
						  //Recherche des deux points les plus proches de xP sur la face
						  Point_3 x1 = F.vertex[0].pos; //Point le plus proche
						  Point_3 x2; // deuxieme point le plus proche
						  double d1 = sqrt(CGAL::to_double(CGAL::squared_distance(xP,x1)));
						  double d2 = 10000000.;
						  for(int k=1;k<F.size();k++){
							Point_3 x = F.vertex[k].pos;
							double d = sqrt(CGAL::to_double(CGAL::squared_distance(xP,x)));
							if(d<d2){
							  if(d<d1){
								x2 = x1;
								d2 = d1;
								x1 = x;
								d1 = d;
							  } else {
								x2 = x;
								d2 = d;
							  }
							}
						  }
						  double d12 = sqrt(CGAL::to_double(CGAL::squared_distance(x1,x2)));
						  if(d1*d1+d12*d12<d2*d2){
							//Le projete est le point le plus proche dans ce cas
							double d = sqrt(CGAL::to_double(CGAL::squared_distance(center_cell,xP)));
							if(d<dist_min){
							  dist_min = d;
							  projete = x1;
							}
						  } else {
							//Le projete est la projection sur (x1,x2)
							Line_3 L(x1,x2);
							double d = sqrt(CGAL::to_double(CGAL::squared_distance(center_cell,L)));
							if(d<dist_min){
							  dist_min = d;
							  projete = L.projection(center_cell);
							}
						  }
						  
						  
						}

					  }
					}
				  }
				  //Calcul du symetrique par rapport au plan defini par centre_face et normale_face
				  Point_3 symm_center = center_cell + Vector_3(center_cell,projete)*2;
				  cm = in_cell(symm_center);
				  c.rho = cm.rho;
				  c.u = cm.u;
				  c.v = cm.v;
				  c.w = cm.w;
				  c.p = cm.p;
				  c.impx = cm.impx;
				  c.impy = cm.impy;
				  c.impz = cm.impz;
				  c.rhoE = cm.rhoE;
				  grille[i][j][k] = c;
				}
			}
		}
	}
}


//Transformation barycentrique du point Xn
Point_3 tr(Triangle_3 Tn, Triangle_3 Tn1, Point_3 Xn){
	
	
	double lambda = 0., mu = 0.;
	
	
	double dom = CGAL::to_double((Tn.operator[](0).operator[](0) - Tn.operator[](2).operator[](0)) * 
	           (Tn.operator[](1).operator[](1) - Tn.operator[](2).operator[](1))-
	           (Tn.operator[](0).operator[](1) - Tn.operator[](2).operator[](1)) * 
	           (Tn.operator[](1).operator[](0) - Tn.operator[](2).operator[](0)));
						 
						 
	double num1 = CGAL::to_double((Xn.operator[](0) - Tn.operator[](2).operator[](0)) * 
						    (Tn.operator[](1).operator[](1) - Tn.operator[](2).operator[](1)) -
						    (Xn.operator[](1) - Tn.operator[](2).operator[](1)) * 
						    (Tn.operator[](1).operator[](0) - Tn.operator[](2).operator[](0)));
								
	double num2 =-1*( CGAL::to_double((Xn.operator[](0) - Tn.operator[](2).operator[](0)) * 
								(Tn.operator[](0).operator[](1) - Tn.operator[](2).operator[](1)) -
								(Xn.operator[](1) - Tn.operator[](2).operator[](1)) * 
								(Tn.operator[](0).operator[](0) - Tn.operator[](2).operator[](0))));
	if(dom>eps){							
	lambda =  num1/dom;
	mu = num2/dom;
	}
	else {std::cout<<"Oupps division par zero dans la fonction tr "<<std::endl;}
	
	double x = CGAL::to_double(lambda * Tn1.operator[](0).operator[](0) + mu*Tn1.operator[](1).operator[](0) 
	                           + (1- lambda- mu)*Tn1.operator[](2).operator[](0));
														 
	double y = CGAL::to_double(lambda * Tn1.operator[](0).operator[](1) + mu*Tn1.operator[](1).operator[](1) 
	                           + (1- lambda- mu)*Tn1.operator[](2).operator[](1));
	double z = CGAL::to_double(lambda * Tn1.operator[](0).operator[](2) + mu*Tn1.operator[](1).operator[](2) 
														 + (1- lambda- mu)*Tn1.operator[](2).operator[](2));		
														 
	
	return Point_3(x, y, z);
}


// Transormation barycentrique du Triangle T 
Triangle_3 tr(Triangle_3 Tn, Triangle_3 Tn1, Triangle_3 T){
	
	Point_3 s = tr( Tn,Tn1, T.operator[](0) );
	Point_3 r = tr( Tn,Tn1, T.operator[](1) );
	Point_3 v = tr( Tn,Tn1, T.operator[](2) );
	
	return Triangle_3(s, r, v);
}
//transformation inverse tr(Tn1,Tn, tr(Tn,Tn1,T))

//Transformation d'un Point_3 en Point_2
Point_2 tr(Triangle_3 Tn1, Point_3 Xn){
		
 double lambda = 0., mu = 0.;
  
 double dom = CGAL::to_double((Tn1.operator[](0).operator[](0) - Tn1.operator[](2).operator[](0)) * 
 (Tn1.operator[](1).operator[](1) - Tn1.operator[](2).operator[](1))-
 (Tn1.operator[](0).operator[](1) - Tn1.operator[](2).operator[](1)) * 
 (Tn1.operator[](1).operator[](0) - Tn1.operator[](2).operator[](0)));
 
 double num1 = CGAL::to_double((Xn.operator[](0) - Tn1.operator[](2).operator[](0)) * 
 (Tn1.operator[](1).operator[](1) - Tn1.operator[](2).operator[](1)) -
 (Xn.operator[](1) - Tn1.operator[](2).operator[](1)) * 
 (Tn1.operator[](1).operator[](0) - Tn1.operator[](2).operator[](0)));
 
 double num2 = -1*(CGAL::to_double((Xn.operator[](0) - Tn1.operator[](2).operator[](0)) * 
 (Tn1.operator[](0).operator[](1) - Tn1.operator[](2).operator[](1)) -
 (Xn.operator[](1) - Tn1.operator[](2).operator[](1)) * 
 (Tn1.operator[](0).operator[](0) - Tn1.operator[](2).operator[](0))));
 if(dom>eps){							
	 lambda =  num1/dom;
	 mu = num2/dom;
 }
 else {std::cout<<"Oupps division par zero dans la fonction tr "<<std::endl;}
	
	Point_2 M(mu, (1-lambda-mu));
	
	return M;
}	

// Transormation d'un Triangle_3 en Triangle_2 
Triangle_2 tr(Triangle_3 Tn1, Triangle_3 T){
	
	Point_2 s = tr( Tn1, T.operator[](0) );
	Point_2 r = tr( Tn1, T.operator[](1) );
	Point_2 v = tr( Tn1, T.operator[](2) );
	
	return Triangle_2(s, r, v);
}



//Transformation d'un Point_2 en Point_3
Point_3 tr(Triangle_3 Tn1, Point_2 Xn){
	
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

// Transormation d'un Triangle_2 en Triangle_3
Triangle_3 tr(Triangle_3 Tn1, Triangle_2 T){
	
	Point_3 s = tr( Tn1, T.operator[](0) );
	Point_3 r = tr( Tn1, T.operator[](1) );
	Point_3 v = tr( Tn1, T.operator[](2) );
	
	return Triangle_3(s, r, v);
}



CDT sous_maillage_face(Triangles_2 Tn1, Triangles_2 Tn_n1, CDT &cdt){
		
	std::vector<Bbox_2> boxesTn1, boxesTn_n1; //tres outil pour les intersection 
	
	for(Triangle2_iterator it= Tn1.begin(); it!= Tn1.end(); ++it){  //on associe a chaque triangle un Box(une boite contenant le triangle)
		boxesTn1.push_back(Bbox_2(it->bbox()));
	}
	
	for(Triangle2_iterator it= Tn_n1.begin(); it!= Tn_n1.end(); ++it){
		boxesTn_n1.push_back(Bbox_2(it->bbox()));
	}
	
	Triangle_2 t;
	Point_2 P;
	Segment_2 seg;
	std::vector<Point_2> vPoints; 
	
	std::vector<Point_2> intPoints; //vector de Point_2 d'intersection
	std::vector<Segment_2> intSeg;  //vector de Segment_2 d'intersection
	
	for(int i=0; i<boxesTn1.size(); i++ ){ 
		for(int j=0; j<boxesTn_n1.size(); j++ ){
			//cout<<"Triangle 1: "<<Tn1[i]<<" Triangle 2: "<<Tn_n1[j]<<endl;
			if (CGAL::do_overlap( boxesTn1[i],boxesTn_n1[j]) ) //test d'intersection des Box 
			{
					if (CGAL::do_intersect(Tn1[i],Tn_n1[j]) ){ // test d'intersection des triangles contenues dans les Box
						
						CGAL::Object result = CGAL::intersection(Tn1[i],Tn_n1[j]); //calcul d'intersection entre les deux triangles
						
						if(CGAL::assign(P,result)){ 
							intPoints.push_back(P);
						}
						else if(CGAL::assign(seg,result)){
							intSeg.push_back(seg);
							
						}
						else if(CGAL::assign(t,result)){
							Segment_2 s1(t.operator[](0), t.operator[](1));
							Segment_2 s2(t.operator[](1), t.operator[](2));
							Segment_2 s3(t.operator[](2), t.operator[](0));
							intSeg.push_back(s1);	
							intSeg.push_back(s2);	
							intSeg.push_back(s3);							
						}
						else if(CGAL::assign(vPoints,result)){ 
							for(int l= 0; l<vPoints.size(); l++)
							{
								intPoints.push_back(vPoints[l]);
							}
							
						}
						else {cout<<"Intersection type: ? "<<endl;
									cout<<"Triangle 1: "<<Tn1[i]<<" Triangle 2: "<<Tn_n1[j]<<endl;
						}
						
					}
				}
			}
		}
	
	 cdt.insert(intPoints.begin(), intPoints.end()); //insertion des points d'intersection dans le maillage
	 
	 for(int i = 0; i<intSeg.size(); i++){
		 //construction du maillage 2d sous la contrainte "intSeg[i] est une arrete dans le maillage"
		 cdt.insert_constraint(intSeg[i].operator[](0), intSeg[i].operator[](1));
	 }
	
	return cdt;
}	


double swap (Triangle_3 Tn, Triangles tn, Triangle_3 Tn1, Triangles tn1){

	double swap=0.;
	// transf barycentrique de tn 
	Triangles tn_n1(tn.size());
	for(int i=0; i<tn.size(); i++){
		tn_n1[i] = tr(Tn, Tn1, tn[i]);
	}
	
	
	Point_2 Ap(0., 0.); 
	Point_2 Bp(1., 0.);
	Point_2 Cp(0., 1.);
	Triangle_2 Ref(Ap,Bp,Cp);
	
	// Transf du Triangles_3  tn_n1 en Triangle_2
	Triangles_2 Tn_n1_2(1+tn_n1.size());
	Tn_n1_2[0] = Ref;
	for(int i=0; i<tn_n1.size(); i++){
		Tn_n1_2[i+1] = tr(Tn1, tn_n1[i]);
	}
	
	
	
	Triangles_2 Tn1_2(1+tn1.size());
	Tn1_2[0] = Ref;
	
	// Transf du Triangles_3  tn1 en Triangle_2
	for(int i=0; i<tn1.size(); i++){
		Tn1_2[i+1] =tr(Tn1, tn1[i]);
	}
	
	
	// sous maillage triangulaire de l'interface
	CDT cdt;
	cdt.insert(Ap); cdt.insert(Bp); cdt.insert(Cp);
	sous_maillage_face(Tn1_2, Tn_n1_2, cdt);
	assert(cdt.is_valid());
	
	Triangles_2 T2d; //recuperation faces du maillage Triangle_2
	
	
	
// 	int count = 0; 
// 	for (CDT::Finite_edges_iterator eit = cdt.finite_edges_begin();
// 	eit != cdt.finite_edges_end();
// 	++eit)
// 	if (cdt.is_constrained(*eit)) ++count;

// Vertex_cdt v;
// for (CDT::Vertex_iterator ver=cdt.vertices_begin();
// ver!=cdt.vertices_end();++ver){
// 	
// 	v= *ver;
// 	std::cout<<v.point()<<std::endl;
// }
	
	int i=0;
	for (CDT::Finite_faces_iterator fit=cdt.finite_faces_begin();
	          fit!=cdt.finite_faces_end();++fit)
	{
		Point_2 s = fit->vertex(0)->point();
		Point_2 v = fit->vertex(1)->point();
		Point_2 r = fit->vertex(2)->point();
		T2d.push_back(Triangle_2(s,v,r));
		//cout<<"triangle intersect "<<T2d[i]<<endl;
		//i++;
	}
	



	//transf des Triangle_2 en Triangle_3
	Triangles T3d(T2d.size());
	for(int i=0; i<T2d.size(); i++){ 
		T3d[i] = tr(Tn1,T2d[i]);
		//cout<<"triangle 3d "<<T3d[i]<<endl;
	}
	
	
	//transf inverse 
	Triangles T3d_n(T3d.size());
	for(int i=0; i<T3d.size(); i++){ 
		T3d_n[i] = tr(Tn1,Tn,T3d[i]);
		//cout<<"triangle 3d n "<<T3d_n[i]<<endl;
	}
	return swap;
}




/*
Triangles_2 tr(Triangle_3 Ref, Triangles T, Triangles_2 &Ttr){
	
	//Ref est le triangle Principal de la surface du solide
	
	Vector_3 AB(Ref.operator[](0), Ref.operator[](1));
	Vector_3 AC(Ref.operator[](0), Ref.operator[](2));
	Vector_3 N = cross_product(AB, AC);
	Vector_3 Np = cross_product(N, AB);
	
	double ABnorm= std::sqrt(CGAL::to_double(AB.squared_length() ));
	double Npnorm= std::sqrt(CGAL::to_double(Np.squared_length() ));
	
	double Hx= (CGAL::to_double(AC*AB))/(ABnorm);
	double Hy= (CGAL::to_double(AC*Np))/(Npnorm);
	
	Point_2 Ap(0., 0.); 
	Point_2 Bp(ABnorm, 0.);
	Point_2 Cp(Hx, Hy);
	
	
	Triangle_2 Refp(Ap,Bp,Cp);
	Ttr.push_back(Refp);
	
	for(int it= 0; it<T.size(); it++){
		
		double Mx= 0., My=0. ;
		double C[3][2];
		for(int j= 0; j<3; j++){
			
			if(abs(T[it].operator[](j).operator[](0) - Ref.operator[](0).operator[](0))<eps &&
				abs(T[it].operator[](j).operator[](1) - Ref.operator[](0).operator[](1))<eps &&
				abs(T[it].operator[](j).operator[](2) - Ref.operator[](0).operator[](2))<eps){
				Mx= 0.;
			My= 0.;
			}
			
			
			else if(abs(T[it].operator[](j).operator[](0) - Ref.operator[](1).operator[](0))<eps &&
				abs(T[it].operator[](j).operator[](1) - Ref.operator[](1).operator[](1))<eps &&
				abs(T[it].operator[](j).operator[](2) - Ref.operator[](1).operator[](2))<eps){
				Mx= ABnorm;
			My= 0.;
			}
			
			else if(abs(T[it].operator[](j).operator[](0) - Ref.operator[](2).operator[](0))<eps &&
				abs(T[it].operator[](j).operator[](1) - Ref.operator[](2).operator[](1))<eps &&
				abs(T[it].operator[](j).operator[](2) - Ref.operator[](2).operator[](2))<eps){
				Mx= Hx;
			My= Hy;
			}
			
			else {
				Vector_3 AM (Ref.operator[](0), T[it].operator[](j)); 
				Mx= (CGAL::to_double(AM*AB))/(ABnorm);
				My= (CGAL::to_double(AM*Np))/(Npnorm);
				
			}
			C[j][0]= Mx; C[j][1]= My; 
		}
		
		Ttr.push_back(Triangle_2(Point_2(C[0][0],C[0][1]),Point_2(C[1][0],C[1][1]),Point_2(C[2][0],C[2][1])));
	}
	
	
	return Ttr;
}	*/


/*Triangulation sous_maillage_face(Triangles Tn1, Triangles Tn_n1){

	CDT cdt;
	
	std::vector<Bbox> boxesTn1(Tn1.size()), boxesTn_n1(Tn_n1.size());
	
	for(Tri_iterator it= Tn1.begin(); it!= Tn1.end(); ++it){
		boxesTn1.push_back(Bbox(it->bbox()));
		}
		
		for(Tri_iterator it= Tn_n1.begin(); it!= Tn_n1.end(); ++it){
			boxesTn_n1.push_back(Bbox(it->bbox()));
		}
		
		Triangle_3 t;
		Point_3 P;
		Segment_3 seg;
		std::vector<Point_3> vPoints; 
		
		std::vector<Point_3> intPoints;
		std::vector<Segment_3> intSeg;
		
		
		for(int i=0; i<boxesTn1.size(); i++ ){
			for(int j=0; j<boxesTn_n1.size(); j++ ){
				if (CGAL::do_intersect( boxesTn1[i],boxesTn_n1[j]) )
				{
					if (CGAL::do_intersect(boxesTn1[i],Tn_n1[j]) ){
						if (CGAL::do_intersect(Tn1[i],Tn_n1[j]) ){
							
							CGAL::Object result = CGAL::intersection(Tn1[i],Tn_n1[j]);
							
							if(CGAL::assign(P,result)){
								intPoints.push_back(P);
							}
							else if(CGAL::assign(seg,result)){
								intSeg.push_back(seg);
								
							}
							else if(CGAL::assign(t,result)){
								Segment_3 s1(t.operator[](0), t.operator[](1));
								Segment_3 s2(t.operator[](1), t.operator[](2));
								Segment_3 s3(t.operator[](2), t.operator[](0));
								intSeg.push_back(s1);	
								intSeg.push_back(s2);	
								intSeg.push_back(s3);							
							}
							else if(CGAL::assign(vPoints,result)){ 
								for(int l= 0; l<vPoints.size(); l++)
								{
									intPoints.push_back(vPoints[l]);
								}
								
						}
						else {cout<<"Intersection type: ? "<<endl;}
						
						}
					}
				}
			}
		}
		
		Triangulation T(intPoints.begin(), intPoints.end());
		//assert( T.dimension() == 2);
		
		if (T.dimension() == 2){
			//cdt.insert();
		}
		else {std::cout<<"il y a un probleme dans la construction du sous-maillage triangulaire de l'interface"<<std::cout;}
		
		return T;
		}*/	