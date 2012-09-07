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
void Grille::fill_cel(std::vector<Solide>& S){
	
	Cellule c, cm;
	int nb_part = S.size();
	double dist[6*nb_part];
	double dist_min = 100;
	int poz=0;
	double x_min=0., y_min=0., z_min = 0., x_max = 0., y_max=0., z_max=0.;

	//Solide solide(x_min,y_min,z_min,x_max,y_max,z_max);
	

	Point_3 center_faces[6][nb_part];
	for(int it=0; it<S.size(); it++){
			center_faces[0][it]= S[it].centre[0]; 
			center_faces[1][it]= S[it].centre[1]; 
			center_faces[2][it]= S[it].centre[2]; 
			center_faces[3][it]= S[it].centre[3]; 
			center_faces[4][it]= S[it].centre[4]; 
			center_faces[5][it]= S[it].centre[5]; 
	}

   Point_3 ref(100., 100., 100.);
		int count=0;
		for(int i=0; i<6; i++){	
			for(int j=0; j<nb_part; j++){
				for(int k=i+1; k<6; k++){
					for(int l=0; l<nb_part; l++){	
					if(CGAL::to_double(squared_distance(center_faces[i][j], center_faces[k][l]))<eps) 
					{ center_faces[i][j] = ref ;
					center_faces[k][l] = ref;}
				 }
				}
				count++;
			}
		}
	//std::cout<<"center faces number: " <<count<<std::endl;
	for(int i=marge;i<Nx+marge;i++){
		for(int j=marge;j<Ny+marge;j++){ 
			for(int k=marge;k<Nz+marge;k++){
				c = grille[i][j][k];
				if((std::abs(c.alpha-1.)<eps))
				{
					Point_3 center_cell(c.x, c.y, c.z);
				  int nbx=0, nby=0,nbz=0;
					double dist_min = 100;
					double nb=0.;
					count = 0;
					for(int it=0; it<6; it++){
						for(int iter=0; iter<nb_part; iter++){
							dist[count] = CGAL::to_double(squared_distance(center_cell, center_faces[it][iter]));
							if(dist[count]< dist_min) {
								dist_min = dist[count];
								poz = it;
							}
							count++;
						}
					}
					if (poz == 0){
						nb = dist_min/c.dx;
						if (nb != (int)(nb)){ nbx= (int)(nb)+1;}
						else {nbx = nb;}
						if(i-2*nbx>0){
						cm = grille[i-2*nbx][j][k]; // a definir la cellule mirroir par rapport à l'interface
						while(cm.alpha>eps && (i-2*nbx)>marge){nbx++; cm = grille[i-2*nbx][j][k]; }
						}
						else {cm = grille[0][j][k]; }
					}
					else if (poz == 1){
						nb = dist_min/c.dx;
						if (nb != (int)(nb)){ nbx= (int)(nb)+1;}
						else {nbx = nb;}
						if(i+2*nbx <Nx+2*marge){
						cm = grille[i+2*nbx][j][k]; // a definir la cellule mirroir par rapport à l'interface
						while(cm.alpha>eps && (i+2*nbx)<Nx+marge){nbx++; cm = grille[i+2*nbx][j][k]; }
						}
						else{cm = grille[Nx+marge][j][k];}
					}
					
					else if (poz == 2){
						nb = dist_min/c.dy;
						if (nb != (int)(nb)){ nby= (int)(nb)+1;}
						else {nby = nb;}
						if(j-2*nby>0){
						cm = grille[i][j-2*nby][k]; // a definir la cellule mirroir par rapport à l'interface
						while(cm.alpha>eps && (j-2*nby)>marge){nby++; cm = grille[i][j-2*nby][k]; }
						}
						else{cm = grille[i][0][k];}
					}
					
					else if (poz == 3){
						nb = dist_min/c.dy;
						if (nb != (int)(nb)){ nby= (int)(nb)+1;}
						else {nby = nb;}
						if(j+2*nby<Ny+2*marge){
						cm = grille[i][j+2*nby][k]; // a definir la cellule mirroir par rapport à l'interface
						while(cm.alpha>eps && (j+2*nby)<Ny+marge){nby++; cm = grille[i][j+2*nby][k]; }
						}
						else{cm = grille[i][Ny+marge][k];}
					}
					else if (poz == 4){
						nb = dist_min/c.dz;
						if (nb != (int)(nb)){ nbz= (int)(nb)+1;}
						else {nbz = nb;}
						if(k-2*nbz>0){
						cm = grille[i][j][k-2*nbz]; // a definir la cellule mirroir par rapport à l'interface
						while(cm.alpha>eps && (k-2*nbz)>marge){nbz++; cm = grille[i][j][k-2*nbz]; }
						}
						else{cm = grille[i][j][0];}
					}
					else{
						nb = dist_min/c.dz;
						if (nb != (int)(nb)){ nbz= (int)(nb)+1;}
						else {nbz = nb;}
						if(k+2*nbz<Nz+2*marge){
						cm = grille[i][j][k+2*nbz]; // a definir la cellule mirroir par rapport à l'interface
						while(cm.alpha>eps && (k+2*nbz)<Nz+marge){nbz++; cm = grille[i][j][k+2*nbz]; }
						}
						else{cm = grille[i][j][Nz+marge];}
					}
					c.rho = cm.rho;
					c.impx = cm.impx;
					c.impy = cm.impy;
					c.impz = cm.impz;
					c.rhoE = cm.rhoE;
					c.u = cm.u ;
					c.v = cm.v;
					c.w = cm.w;
					c.p = cm.p;
// 					{std::cout<<" center cell "<< center_cell<< "distance min "<< dist_min<< "position "<<poz<<std::endl;
// 					std::cout<<"cellule miroir "<<std::endl; cm.Affiche();}
					grille[i][j][k] = c;
				}
			}
		}
	}
}


