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
				  Point_3 projete(0.,0.,0.); //Projeté sur la face la plus proche
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


