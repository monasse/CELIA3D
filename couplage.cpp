#include "fluide.hpp"
#include "intersections.cpp"

//Calcul des Forces fluides et Moments fluides exerces sur le solide	
void Grille::Forces_fluide(Solide& S, const double dt){
	
	//Mise à jour des Forces fluides et Moments fluides exerces sur le solide 
	for(int iter_s=0; iter_s<S.size(); iter_s++){ 
		
		S.solide[iter_s].Ffprev = S.solide[iter_s].Ff;
		S.solide[iter_s].Mfprev = S.solide[iter_s].Mf;
		Point_3 Xn(S.solide[iter_s].x0.operator[](0) + S.solide[iter_s].Dx.operator[](0), S.solide[iter_s].x0.operator[](1) + S.solide[iter_s].Dx.operator[](1),S.solide[iter_s].x0.operator[](2) + S.solide[iter_s].Dx.operator[](2));
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
		
		S.solide[iter_s].Ff = Vector_3(fx,fy,fz);
		S.solide[iter_s].Mf = Vector_3(CGAL::to_double(mx),CGAL::to_double(my),CGAL::to_double(mz)); 
		//cout<<" Ff "<< S.solide[iter_s].Ff <<endl;
		//cout<<" Mf "<< S.solide[iter_s].Mf <<endl;
	} //fin boucle sur les particules
	
}	


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
					c.flux_modif[1] = c.phi_x;
					c.flux_modif[2] = c.phi_y;
					c.flux_modif[3] = c.phi_z;
					c.flux_modif[4] = c.phi_v;
					for(int l=0.; l<5; l++){  
						c.flux_modif[l] -= (1.-c.kappai)*c.dtfxi[l] - (1.-ci.kappai)*ci.dtfxi[l]
						+ (1.-c.kappaj)*c.dtfyj[l] - (1.-cj.kappaj)*cj.dtfyj[l]
						+ (1.-c.kappak)*c.dtfzk[l] - (1.-ck.kappak)*ck.dtfzk[l] - c.delta_w[l];
						if(std::abs(c.flux_modif[l])>eps){
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
	
	for(int i=marge;i<Nx+marge;i++){
		for(int j=marge;j<Ny+marge;j++){ 
			for(int k=marge;k<Nz+marge;k++){
				cp = grille[i][j][k];
				bool test=true;
				if((cp.alpha>epsa || cp.p<0. || cp.rho<0.) && abs(cp.alpha-1.)>eps){
					
					for(int ii=-1; ii<=1 && test; ii++){
						for(int jj=-1; jj<=1 && test; jj++){
							for(int kk=-1; kk<=1 && test; kk++){
								if (grille[i+ii][j+jj][k+kk].alpha <eps && grille[i+ii][j+jj][k+kk].p>0. && grille[i+ii][j+jj][k+kk].rho>0.)
								{
									test=false;
									cg = grille[i+ii][j+jj][k+kk];
									
									cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
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
									
									grille[i][j][k] = cp;
									grille[i+ii][j+jj][k+kk] = cg;
								} //if cg.alpha==0
								
							}
						}
					}
					
					if(test){
						
						if (grille[i-2][j][k].alpha == 0. && grille[i-2][j][k].p>0. && grille[i-2][j][k].rho>0.)
					{
							cg = grille[i-2][j][k];
							
							cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
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
							
							grille[i][j][k] = cp;
							grille[i-2][j][k] = cg;
							test = false;
						}
						else if (grille[i+2][j][k].alpha == 0. && grille[i+2][j][k].p>0. && grille[i+2][j][k].rho>0.)
						{
								cg = grille[i+2][j][k];
								
								cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
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
								
								grille[i][j][k] = cp;
								grille[i+2][j][k] = cg;
								test = false;
							}
							
							else if (grille[i][j-2][k].alpha == 0. && grille[i][j-2][k].p>0. && grille[i][j-2][k].rho>0.)
							{
									cg = grille[i][j-2][k];
									
									cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
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
									
									grille[i][j][k] = cp;
									grille[i][j-2][k] = cg;
									test = false;
								}
								else if (grille[i][j+2][k].alpha == 0. && grille[i][j+2][k].p>0. && grille[i][j+2][k].rho>0.)
								{
										cg = grille[i][j+2][k];
										
										cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
										cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
										cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
										cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
										cp.MrhoE = (cg.rhoE - cp.rhoE)/(2.- cp.alpha);
										
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
										
										grille[i][j][k] = cp;
										grille[i][j+2][k] = cg;
										test = false;
									}
									else if (grille[i][j][k-2].alpha == 0. && grille[i][j][k-2].p>0. && grille[i][j][k-2].rho>0.)
									{
											cg = grille[i][j][k-2];
											
											cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
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
											
											grille[i][j][k] = cp;
											grille[i][j][k-2] = cg;
											test = false;
										}
										else if(grille[i][j][k+2].alpha == 0. && grille[i][j][k+2].p>0. && grille[i][j][k+2].rho>0.)
										{
												cg = grille[i][j][k+2];
												
												cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
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
												
												grille[i][j][k] = cp;
												grille[i][j][k+2] = cg;
												test = false;
											}
						
					}//fin if(test)
					else if(test){
						std::cout<<"Pas de cellule cible pour le mixage"<<std::endl; 
						std::cout<< "position du centre de la cellule : "<<grille[i][j][k].x << " "<<grille[i][j][k].y << " "<<grille[i][j][k].z << " "<< " rho "<<grille[i][j][k].rho  << " p "<<grille[i][j][k].p <<" alpha " << grille[i][j][k].alpha<<std::endl;
						std::cout<<"cellules voisines : "<<std::endl;
						for(int ii=-1; ii<=1 && test; ii++){
							for(int jj=-1; jj<=1 && test; jj++){
								for(int kk=-1; kk<=1 && test; kk++){
									std::cout<<"alpha "<<grille[i+ii][j+jj][k+kk].alpha<< "  "<< " rho "<<grille[i+ii][j+jj][k+kk].rho << "p "<< grille[i+ii][j+jj][k+kk].p<<std::endl; 
								}
							}
						}
					} //fin else if(test)
				}// 0.5<c.alpha<1.
			} //fin boucle sur la grille
		}
	}
	
}





// void Grille:: mixage(){
// 	
// 	Cellule cp, cg;
// 	bool test=true;
// 	for(int i=marge;i<Nx+marge;i++){
// 		for(int j=marge;j<Ny+marge;j++){ 
// 			for(int k=marge;k<Nz+marge;k++){
// 				cp = grille[i][j][k];
// 				if((cp.alpha>epsa || cp.p<0. || cp.rho<0.) && abs(cp.alpha-1.)>eps){
// 					
// 					if (grille[i-1][j][k].alpha == 0.)
// 					{
// 						cg = grille[i-1][j][k];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 						
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						 
// 						grille[i][j][k] = cp;
// 						grille[i-1][j][k] = cg;
// 					}
// 					else if (grille[i+1][j][k].alpha == 0.)
// 					{
// 						cg = grille[i+1][j][k];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 
// 						
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i+1][j][k] = cg;
// 					}
// 					
// 					else if (grille[i][j-1][k].alpha == 0.)
// 					{
// 						cg = grille[i][j-1][k];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i][j-1][k] = cg;
// 					}
// 					else if (grille[i][j+1][k].alpha == 0.)
// 					{
// 						cg = grille[i][j+1][k];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2.- cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i][j+1][k] = cg;
// 					}
// 					else if (grille[i][j][k-1].alpha == 0.)
// 					{
// 						cg = grille[i][j][k-1];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 						
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i][j][k-1] = cg;
// 					}
// 					else if(grille[i][j][k+1].alpha == 0.)
// 					{
// 						cg = grille[i][j][k+1];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 	
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i][j][k+1] = cg;
// 					}
// 					
// 					else if (grille[i-2][j][k].alpha == 0.)
// 					{
// 						cg = grille[i-2][j][k];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i-2][j][k] = cg;
// 					}
// 					else if (grille[i+2][j][k].alpha == 0.)
// 					{
// 						cg = grille[i+2][j][k];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 						
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i+2][j][k] = cg;
// 					}
// 					
// 					else if (grille[i][j-2][k].alpha == 0.)
// 					{
// 						cg = grille[i][j-2][k];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i][j-2][k] = cg;
// 					}
// 					else if (grille[i][j+2][k].alpha == 0.)
// 					{
// 						cg = grille[i][j+2][k];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2.- cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 						
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i][j+2][k] = cg;
// 					}
// 					else if (grille[i][j][k-2].alpha == 0.)
// 					{
// 						cg = grille[i][j][k-2];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i][j][k-2] = cg;
// 					}
// 					else if(grille[i][j][k+2].alpha == 0.)
// 					{
// 						cg = grille[i][j][k+2];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i][j][k+2] = cg;
// 					}
// 					else if(grille[i-1][j-1][k].alpha <= eps_relat)
// 					{
// 						cg = grille[i-1][j-1][k];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 						
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i-1][j-1][k] = cg;
// 					}
// 					else if(grille[i-1][j+1][k].alpha <= eps_relat)
// 					{
// 						cg = grille[i-1][j+1][k];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 						
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i-1][j+1][k] = cg;
// 					}
// 					else if(grille[i+1][j-1][k].alpha <= eps_relat)
// 					{
// 						cg = grille[i+1][j-1][k];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
// 						
// 						
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i+1][j-1][k] = cg;
// 					}
// 					else if(grille[i+1][j+1][k].alpha <= eps_relat)
// 					{
// 						cg = grille[i+1][j+1][k];
// 						
// 						cp.Mrho = (cg.rho - cp.rho)/(2. - cp.alpha) ;
// 						cp.Mimpx = (cg.impx - cp.impx)/(2. - cp.alpha);
// 						cp.Mimpy = (cg.impy - cp.impy)/(2. - cp.alpha);
// 						cp.Mimpz = (cg.impz - cp.impz)/(2. - cp.alpha);
// 						cp.MrhoE = (cg.rhoE - cp.rhoE)/(2. - cp.alpha);
// 						
// 						cg.Mrho = (1.-cp.alpha)*(cp.rho - cg.rho)/(2. - cp.alpha) ;
// 						cg.Mimpx = (1.-cp.alpha)*(cp.impx - cg.impx)/(2. - cp.alpha);
// 						cg.Mimpy = (1.-cp.alpha)*(cp.impy - cg.impy)/(2. - cp.alpha);
// 						cg.Mimpz = (1.-cp.alpha)*(cp.impz - cg.impz)/(2. - cp.alpha);
// 						cg.MrhoE = (1.-cp.alpha)*(cp.rhoE - cg.rhoE)/(2. - cp.alpha);
// 						
// 						cp.rho += cp.Mrho;
// 						cp.impx += cp.Mimpx;
// 						cp.impy += cp.Mimpy;
// 						cp.impz += cp.Mimpz;
// 						cp.rhoE += cp.MrhoE;
// 						cp.u = cp.impx/cp.rho;
// 						cp.v = cp.impy/cp.rho;
// 						cp.w = cp.impz/cp.rho;
// 						cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);						
// 
// 						cg.rho += cg.Mrho;
// 						cg.impx += cg.Mimpx;
// 						cg.impy += cg.Mimpy;
// 						cg.impz += cg.Mimpz;
// 						cg.rhoE += cg.MrhoE;
// 						cg.u = cg.impx/cg.rho;
// 						cg.v = cg.impy/cg.rho;
// 						cg.w = cg.impz/cg.rho;
// 						cg.p = (gam-1.)*(cg.rhoE-cg.rho*cg.u*cg.u/2.-cg.rho*cg.v*cg.v/2. - cg.rho*cg.w*cg.w/2.);
// 						
// 						grille[i][j][k] = cp;
// 						grille[i+1][j+1][k] = cg;
// 					}
// 					else {
// 						if(test){std::cout<<"Pas de cellule cible pour le mixage"<<std::endl; 
// 										std::cout<<" grille qui faut melanger "<<std::endl;
// 										cp.Affiche();
// 										test=false;
// 						}
// 					}
// 					
// 				}
// 				
// 			}
// 		}
// 	}
// 	
// } 


void Grille::fill_cel(Solide& S){
	
	Cellule c, cm;
	int nb_part = S.size();
	double dist[6*nb_part];
	double dist_min = 100;
	int poz=0;
	double x_min=0., y_min=0., z_min = 0., x_max = 0., y_max=0., z_max=0.;
	Bbox Fluide(X0,Y0,Z0,X0+domainex,Y0+domainey,Z0+domainez);

	//std::cout<<"center faces number: " <<count<<std::endl;
	for(int i=marge;i<Nx+marge;i++){
		for(int j=marge;j<Ny+marge;j++){
		  for(int k=marge;k<Nz+marge;k++){
				Triangle_3 Tri;
				c = grille[i][j][k];
				if((std::abs(c.alpha-1.)<eps)){
				  Point_3 center_cell(c.x, c.y, c.z);
				  int nbx=0, nby=0,nbz=0;
				  Point_3 projete(0.,0.,0.); //Projete sur la face la plus proche
					Vector_3 V_f(0.,0.,0.); //Vitesse de la paroi au point projete
				  double dist_min = 10000000.;
				  for(int iter=0; iter<nb_part; iter++){
					for(int it=0;it<S.solide[iter].size();it++){
					  if(S.solide[iter].fluide[it]){
							Plane_3 P(S.solide[iter].triangles[it].operator[](0),S.solide[iter].triangles[it].operator[](1),S.solide[iter].triangles[it].operator[](2));
						for(int k=3;k<S.solide[iter].triangles.size() && P.is_degenerate();k++){//Test si le plan est degenere
							P = Plane_3(S.solide[iter].triangles[it].operator[](0),S.solide[iter].triangles[it].operator[](1),S.solide[iter].triangles[it].operator[](k));
						}
						Point_3 xP = P.projection(center_cell);
						//Test pour savoir si le projete est dans la face
						bool test = true;
						for(int k=0;k<2 && test;k++){
							Point_3 x1 = S.solide[iter].triangles[it].operator[](k);
							Point_3 x2 = S.solide[iter].triangles[it].operator[](k+1);
						  Vector_3 vect1(xP,x1);
						  Vector_3 vect2(xP,x2);
						  if(CGAL::to_double(CGAL::cross_product(vect1,vect2)*S.solide[iter].normales[it])<0.){
							test = false;
						  }
						}
						//1er cas : on est dans la face
						if(test){
						  double d = sqrt(CGAL::to_double(CGAL::squared_distance(center_cell,xP)));
						  if(d<dist_min && inside_box(Fluide,xP)){
							dist_min = d;
							projete = xP;
							V_f = S.solide[iter].vitesse_parois(xP);
						  }
						}
						//2eme cas : on est hors de la face
						else{
						  //Recherche du point le plus proche sur toutes les aretes
						  for(int k=0;k<3;k++){
							int kp = (k+1)%3;
							Point_3 x1 = S.solide[iter].triangles[it].operator[](k);
							Point_3 x2 = S.solide[iter].triangles[it].operator[](kp);
							double d1 = sqrt(CGAL::to_double(CGAL::squared_distance(center_cell,x1)));
							double d2 = sqrt(CGAL::to_double(CGAL::squared_distance(center_cell,x2)));
							double d12 = sqrt(CGAL::to_double(CGAL::squared_distance(x1,x2)));
							//1er sous-cas : on est plus proche du point x1
							if(d1*d1+d12*d12<d2*d2){
							  if(d1<dist_min && inside_box(Fluide,x1)){
								dist_min = d1;
								projete = x1;
								V_f = S.solide[iter].vitesse_parois(x1);
							  }
							}
							//2eme sous-cas : on est plus proche du point x2
							else if(d2*d2+d12*d12<d1*d1){
							  if(d2<dist_min && inside_box(Fluide,x2)){
								dist_min = d2;
								projete = x2;
								V_f = S.solide[iter].vitesse_parois(x2);
							  }
							}
							//3eme sous-cas : on prend le projete sur (x1,x2)
							else {
							  Line_3 L(x1,x2);
							  double d = sqrt(CGAL::to_double(CGAL::squared_distance(center_cell,L)));
							  Point_3 proj = L.projection(center_cell);
							  if(d<dist_min && inside_box(Fluide,proj)){
								dist_min = d;
								projete = proj;
								V_f = S.solide[iter].vitesse_parois(proj);
							  }
							}
							}//Recherche du point le plus proche sur toutes les aretes
						}
					  }
					}
				  }
				  //Calcul du symetrique par rapport au plan defini par centre_face et normale_face
				  Point_3 symm_center = center_cell + Vector_3(center_cell,projete)*2;
				  Vector_3 normale(center_cell,projete);
				  double norme = sqrt(CGAL::to_double(normale.squared_length()));
          assert(norme!= 0.);
				  normale = normale*1./norme;
				  cm = in_cell(symm_center);
				  Vector_3 vit_m(cm.u,cm.v,cm.w); //Vitesse au point miroir
					Vector_3 vit = vit_m - normale*2.*((vit_m-V_f)*normale);
				  c.rho = cm.rho;
				  c.u = CGAL::to_double(vit.operator[](0));
				  c.v = CGAL::to_double(vit.operator[](1));
				  c.w = CGAL::to_double(vit.operator[](2));
				  c.p = cm.p;
				  c.impx = c.rho*c.u;
				  c.impy = c.rho*c.v;
				  c.impz = c.rho*c.w;
				  c.rhoE = c.rho/2.*(c.u*c.u+c.v*c.v+c.w*c.w)+c.p/(gam-1.);
				  grille[i][j][k] = c;
				}
			}
		}
	}
}


 
 //Transformation barycentrique du point Xn
 Point_3 tr2(Triangle_3 Tn, Triangle_3 Tn1, Point_3 Xn){
 	
 	
 	double lambda = 0., mu = 0.;
 	
	double dom =CGAL::to_double(cross_product(Vector_3(Tn.operator[](2),Tn.operator[](0)),Vector_3(Tn.operator[](2),
	                           Tn.operator[](1))).operator[](0));
	double num1 =CGAL::to_double(cross_product(Vector_3(Tn.operator[](2),Xn),Vector_3(Tn.operator[](2),
	Tn.operator[](1))).operator[](0));
	double num2 =CGAL::to_double(cross_product(Vector_3(Tn.operator[](2),Tn.operator[](0)),
	Vector_3(Tn.operator[](2),Xn)).operator[](0));
	if(std::abs(dom)<eps){
		dom =CGAL::to_double(cross_product(Vector_3(Tn.operator[](2),Tn.operator[](0)),Vector_3(Tn.operator[](2),
		Tn.operator[](1))).operator[](1));
		num1 =CGAL::to_double(cross_product(Vector_3(Tn.operator[](2),Xn),Vector_3(Tn.operator[](2),
		Tn.operator[](1))).operator[](1));
		num2 =CGAL::to_double(cross_product(Vector_3(Tn.operator[](2),Tn.operator[](0)),
		Vector_3(Tn.operator[](2),Xn)).operator[](1));
		
		if(std::abs(dom)<eps){
			dom =CGAL::to_double(cross_product(Vector_3(Tn.operator[](2),Tn.operator[](0)),Vector_3(Tn.operator[](2),
			Tn.operator[](1))).operator[](2));
			num1 =CGAL::to_double(cross_product(Vector_3(Tn.operator[](2),Xn),Vector_3(Tn.operator[](2),
			Tn.operator[](1))).operator[](2));
			num2 =CGAL::to_double(cross_product(Vector_3(Tn.operator[](2),Tn.operator[](0)),
			Vector_3(Tn.operator[](2),Xn)).operator[](2));
			}
		}

 	lambda =  num1/dom;
 	 mu = num2/dom;


 	double x =CGAL::to_double( lambda * Tn1.operator[](0).operator[](0) + mu*Tn1.operator[](1).operator[](0) 
 	                           + (1- lambda- mu)*Tn1.operator[](2).operator[](0));
 														 
 	double y = CGAL::to_double(lambda * Tn1.operator[](0).operator[](1) + mu*Tn1.operator[](1).operator[](1) 
 	                           + (1- lambda- mu)*Tn1.operator[](2).operator[](1));
 	double z = CGAL::to_double(lambda * Tn1.operator[](0).operator[](2) + mu*Tn1.operator[](1).operator[](2) 
 														 + (1- lambda- mu)*Tn1.operator[](2).operator[](2));		
 														 
 	
 	return Point_3(x, y, z);
 }

// Transormation barycentrique du Triangle T 
Triangle_3 tr2(Triangle_3 Tn, Triangle_3 Tn1, Triangle_3 T){
	
	Point_3 s = tr2( Tn,Tn1, T.operator[](0) );
	Point_3 r = tr2( Tn,Tn1, T.operator[](1) );
	Point_3 v = tr2( Tn,Tn1, T.operator[](2) );
	
	return Triangle_3(s, r, v);
}

//Transformation barycentrique du point Xn
Point_3 tr(Triangle_3 Tn, Triangle_3 Tn1, Point_3 Xn){
	
	
	  double lambda = 0., mu = 0.;
	
		double dom = std::sqrt(CGAL::to_double(cross_product(Vector_3(Tn.operator[](2),Tn.operator[](0)),
											Vector_3(Tn.operator[](2),Tn.operator[](1))).squared_length() ));
		double num1 = std::sqrt(CGAL::to_double(cross_product(Vector_3(Tn.operator[](2),Xn),
											 Vector_3(Tn.operator[](2),Tn.operator[](1))).squared_length() ));
		double num2 = std::sqrt(CGAL::to_double(cross_product(Vector_3(Tn.operator[](2),Tn.operator[](0)),
											 Vector_3(Tn.operator[](2),Xn)).squared_length() )); 

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
  
//  double dom = CGAL::to_double((Tn1.operator[](0).operator[](0) - Tn1.operator[](2).operator[](0)) * 
//  (Tn1.operator[](1).operator[](1) - Tn1.operator[](2).operator[](1))-
//  (Tn1.operator[](0).operator[](1) - Tn1.operator[](2).operator[](1)) * 
//  (Tn1.operator[](1).operator[](0) - Tn1.operator[](2).operator[](0)));
//  
//  double num1 = CGAL::to_double((Xn.operator[](0) - Tn1.operator[](2).operator[](0)) * 
//  (Tn1.operator[](1).operator[](1) - Tn1.operator[](2).operator[](1)) -
//  (Xn.operator[](1) - Tn1.operator[](2).operator[](1)) * 
//  (Tn1.operator[](1).operator[](0) - Tn1.operator[](2).operator[](0)));
//  
//  double num2 = -1*(CGAL::to_double((Xn.operator[](0) - Tn1.operator[](2).operator[](0)) * 
//  (Tn1.operator[](0).operator[](1) - Tn1.operator[](2).operator[](1)) -
//  (Xn.operator[](1) - Tn1.operator[](2).operator[](1)) * 
//  (Tn1.operator[](0).operator[](0) - Tn1.operator[](2).operator[](0))));
 
 	double dom =std::sqrt(CGAL::to_double(cross_product(Vector_3(Tn1.operator[](2),Tn1.operator[](0)),
										    Vector_3(Tn1.operator[](2),Tn1.operator[](1))).squared_length() ));
	double num1 =std::sqrt(CGAL::to_double(cross_product(Vector_3(Tn1.operator[](2),Xn),
										     Vector_3(Tn1.operator[](2),Tn1.operator[](1))).squared_length() ));
	double num2 =std::sqrt(CGAL::to_double(cross_product(Vector_3(Tn1.operator[](2),Tn1.operator[](0)),
												 Vector_3(Tn1.operator[](2),Xn)).squared_length() )); 
												 
 
	 lambda =  num1/dom;
	 mu = num2/dom;

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



 CDT sous_maillage_face(Triangles_2& Tn1, Triangles_2& Tn_n1, CDT &cdt){
 	
 	std::vector<Bbox_2> boxesTn1, boxesTn_n1; //tres outil pour les intersections 
 	
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
 	std::vector<Segment_2> intSeg;  //vector de Segment_2 d'intersection a conserver dans la triangularisation de la face 
 	CGAL::Timer user_time, user_time2;
 	user_time.start();
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
 						else {cout<<"Intersection type: ? in sous_maillage_face"<<endl;
 									cout<<"Triangle 1: "<<Tn1[i]<<" Triangle 2: "<<Tn_n1[j]<<endl;
 						}
 						
 					}
 				}
 			}
 		}
 		//cout << "Intersection triangles 2d pour une face time is: " << user_time.time() << " seconds." <<"nb des triangles " <<boxesTn1.size()+boxesTn_n1.size() <<endl;
 		
 		user_time2.start();
 	  cdt.insert(intPoints.begin(), intPoints.end()); //insertion des points d'intersection dans le maillage
 	 
 	 for(int i = 0; i<intSeg.size(); i++){
 		 //construction du maillage 2d sous la contrainte "intSeg[i] est une arrete dans le maillage"
 		 cdt.insert_constraint(intSeg[i].operator[](0), intSeg[i].operator[](1));
 	 }
 	 //cout << "construction sous-maillage sous contrainte pour une face time is: " << user_time2.time() << " seconds." <<endl;
 	return cdt;
 }	
 
 


double volume_prisme(const Triangle_3& T1,const Triangle_3& T2){
	
	double volume=0.;
	
	Vector_3 V = 2.*cross_product( Vector_3(T1.operator[](0),T1.operator[](1)), Vector_3(T1.operator[](0),T1.operator[](2)) )
	             +2*cross_product( Vector_3(T2.operator[](0),T2.operator[](1)), Vector_3(T2.operator[](0),T2.operator[](2)) )
	             +cross_product( Vector_3(T1.operator[](0),T1.operator[](1)), Vector_3(T2.operator[](0),T2.operator[](2)) )
	             +cross_product( Vector_3(T2.operator[](0),T2.operator[](1)), Vector_3(T1.operator[](0),T1.operator[](2)) );
	
	volume = CGAL::to_double((Vector_3(T1.operator[](0),T2.operator[](0)) + Vector_3(T1.operator[](1),T2.operator[](1)) + 
	                          Vector_3(T1.operator[](2),T2.operator[](2)))*V);
	
	volume /=36;
	return volume;
}

double volume_tetra(const Tetrahedron& Tet){
	
	double volume=0.;
	
	Vector_3 V = cross_product( Vector_3(Tet.operator[](0),Tet.operator[](1)), Vector_3(Tet.operator[](0),Tet.operator[](2)) );
	
	
	volume = CGAL::to_double( (Vector_3(Tet.operator[](0),Tet.operator[](3)) *V) );
	volume /= 6.;
	
	return volume;
}


void Grille::cells_intersection_face(int& in,int& jn,int& kn,int& in1,int& jn1,int& kn1, std::vector<Bbox>& box_cells, std::vector<Cellule>& Cells){
	
	if((in!=in1 && jn==jn1 && kn==kn1)|| (in==in1 && jn!=jn1 && kn==kn1) || (in==in1 && jn==jn1 && kn!=kn1))
	{
		Cellule c0 = grille[in][jn][kn];
		Cellule c1 = grille[in1][jn1][kn1];
		Cells.push_back(c0); Cells.push_back(c1);
		box_cells.push_back(Bbox(c0.x -c0.dx/2.,c0.y -c0.dy/2.,c0.z -c0.dz/2.,
											c0.x +c0.dx/2.,c0.y +c0.dy/2.,c0.z + c0.dz/2.));
		box_cells.push_back(Bbox(c1.x -c1.dx/2.,c1.y -c1.dy/2.,c1.z -c1.dz/2.,
											c1.x +c1.dx/2.,c1.y +c1.dy/2.,c1.z + c1.dz/2.));
	}
	else if(in!=in1 && jn!=jn1 && kn==kn1){
		
		Cellule c0 = grille[in][jn][kn];
		Cellule c1 = grille[in1][jn1][kn1];
		Cellule c2 = grille[in1][jn][kn];
		Cellule c3 = grille[in][jn1][kn1];
		
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
		
		Cellule c0 = grille[in][jn][kn];
		Cellule c1 = grille[in1][jn1][kn1];
		Cellule c2 = grille[in][jn][kn1];
		Cellule c3 = grille[in1][jn1][kn];
		
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
		
		Cellule c0 = grille[in][jn][kn];
		Cellule c1 = grille[in1][jn1][kn1];
		Cellule c2 = grille[in][jn1][kn];
		Cellule c3 = grille[in1][jn][kn1];
		
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
	
		Cellule c0 = grille[in][jn][kn];
		Cellule c1 = grille[in1][jn1][kn1];
		Cellule c2 = grille[in1][jn][kn];
		Cellule c3 = grille[in][jn1][kn];
		Cellule c4 = grille[in][jn][kn1];
		Cellule c5 = grille[in1][jn1][kn];
		Cellule c6 = grille[in1][jn][kn1];
		Cellule c7 = grille[in][jn1][kn1];
		
		
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

void Grille::swap_face(Triangles& T3d_prev, Triangles& T3d_n, const double dt){
	
	CGAL::Timer user_time, user_time2;
	double time=0.;
	std::vector<Bbox> box_prismes(T3d_prev.size());
	for (int i=0; i< T3d_prev.size(); i++){
		Bbox box_triangles_prev = T3d_prev[i].bbox();
		Bbox box_triangles_n = T3d_n[i].bbox();
		box_prismes[i]= box_triangles_prev.operator+(box_triangles_n);
	} 
	int nb=0, nb2=0;  double diff=0.;
	for (int i=0; i< box_prismes.size(); i++){
		double vol_test=0.; 
		int in=0, jn=0, kn=0, in1=0, jn1=0, kn1=0;
		bool interieur = true;
		Point_3 center_prev= centroid(T3d_prev[i].operator[](0),T3d_prev[i].operator[](1),T3d_prev[i].operator[](2));
		Point_3 center_n= centroid(T3d_n[i].operator[](0),T3d_n[i].operator[](1),T3d_n[i].operator[](2));
		in_cell(center_prev, in, jn, kn, interieur);
		in_cell(center_n, in1, jn1, kn1, interieur);
		
		Cellule c_cur= grille[in1][jn1][kn1];
		if((std::abs(c_cur.alpha -1.)<eps)  && (interieur==true)){
			double x= CGAL::to_double(center_n.operator[](0));
			double y= CGAL::to_double(center_n.operator[](1));
			double z= CGAL::to_double(center_n.operator[](2));
			Cellule cd= grille[in1+1][jn1][kn1];
			if (cd.is_in_cell(x,y,z) ) {in1=in1+1;}
			else{
				Cellule cg= grille[in1-1][jn1][kn1];
				if (cg.is_in_cell(x,y,z)) {in1=in1-1;}
				else{
					Cellule ch= grille[in1][jn1+1][kn1];
					if (ch.is_in_cell(x,y,z)) {jn1=jn1+1;}
					else{
						Cellule cb= grille[in1][jn1-1][kn1];
						if (cb.is_in_cell(x,y,z)) {jn1=jn1-1;}
						else{
							Cellule cd= grille[in1][jn1][kn1+1];
							if (cd.is_in_cell(x,y,z)) {kn1=kn1+1;}
							else{
								Cellule cder= grille[in1][jn1][kn1-1];
								if (cder.is_in_cell(x,y,z)) {kn1=kn1-1;}
							}
						}
					}
				}
			}
		} // end if alpha==1
		  
    Cellule c= grille[in1][jn1][kn1];

		if ( (in==in1) && (jn==jn1) && (kn==kn1) && (interieur==true)){
			// le prisme est contenu dans une seule cellule 
			double volume_p=volume_prisme(T3d_prev[i],T3d_n[i]);
			//calcul du volume
			if( (std::abs(volume_p)>eps) && (std::abs(1.-c.alpha)>eps)){
			c.delta_w[0] += volume_p*c.rho0/(c.dx*c.dy*c.dz); 
			c.delta_w[1] += volume_p*c.impx0/(c.dx*c.dy*c.dz);
			c.delta_w[2] += volume_p*c.impy0/(c.dx*c.dy*c.dz); 
			c.delta_w[3] += volume_p*c.impz0/(c.dx*c.dy*c.dz); 
			c.delta_w[4] += volume_p*c.rhoE0/(c.dx*c.dy*c.dz);
			grille[in1][jn1][kn1] = c;
			}
			//vol_test += volume_p;
		}	
		else if((std::abs(volume_prisme(T3d_prev[i],T3d_n[i])) >eps)  && (interieur==true) && (std::abs(1.-c.alpha)>eps)) {
		std::vector<Bbox> box_cells;
		std::vector<Cellule> Cells ;
		cells_intersection_face(in,jn,kn,in1,jn1,kn1,box_cells,Cells);

		//definition Tetraedres 
		std::vector<Tetrahedron> vect_Tet;
		std::vector<Bbox> box_Tet;
		Point_3 e = centroid(T3d_prev[i].operator[](1),T3d_n[i].operator[](1), T3d_prev[i].operator[](2),T3d_n[i].operator[](2));
		Point_3 f = centroid(T3d_prev[i].operator[](0),T3d_n[i].operator[](0), T3d_prev[i].operator[](2),T3d_n[i].operator[](2));
		Point_3 g = centroid(T3d_prev[i].operator[](0),T3d_n[i].operator[](0), T3d_prev[i].operator[](1),T3d_n[i].operator[](1));

		Tetrahedron tet0 (T3d_prev[i].operator[](0),T3d_n[i].operator[](0), g, f);
		//if(!tet0.is_degenerate ()){
			if(abs(tet0.volume ())>eps){
				vect_Tet.push_back(tet0);
				box_Tet.push_back(tet0.bbox());
			}
			
		Tetrahedron tet1 (T3d_prev[i].operator[](1),T3d_n[i].operator[](1), e, g);
		//if(!tet1.is_degenerate ()){
		if(abs(tet1.volume ())>eps){
		vect_Tet.push_back(tet1);
		box_Tet.push_back(tet1.bbox());
		}

		Tetrahedron tet2 (T3d_prev[i].operator[](2),T3d_n[i].operator[](2), f, e);
		//if(!tet2.is_degenerate ()){
		if(abs(tet2.volume ())>eps){
			vect_Tet.push_back(tet2);
			box_Tet.push_back(tet2.bbox());
		}

		Tetrahedron tet3 (T3d_prev[i].operator[](0),T3d_prev[i].operator[](1), T3d_prev[i].operator[](2), g);
		if(abs(tet3.volume ())>eps){
			vect_Tet.push_back(tet3);
			box_Tet.push_back(tet3.bbox());
		}

		Tetrahedron tet4 (T3d_prev[i].operator[](0),g, T3d_prev[i].operator[](2), f);
		// 			if(!tet4.is_degenerate ()){
		if(abs(tet4.volume ())>eps){
		vect_Tet.push_back(tet4);
		box_Tet.push_back(tet4.bbox());
		}
		Tetrahedron tet5 (T3d_prev[i].operator[](1),e, T3d_prev[i].operator[](2), g);
		//if(!tet5.is_degenerate ()){
		if(abs(tet5.volume ())>eps){
			vect_Tet.push_back(tet5);
			box_Tet.push_back(tet5.bbox());
		}

		Tetrahedron tet6 (e,g,f, T3d_prev[i].operator[](2));
		//if(!tet6.is_degenerate ()){
		if(abs(tet6.volume ())>eps){
			vect_Tet.push_back(tet6);
			box_Tet.push_back(tet6.bbox());
		}

		Tetrahedron tet7 (e,f,g, T3d_n[i].operator[](2));
		//if(!tet7.is_degenerate ()){
		if(abs(tet7.volume ())>eps){
			vect_Tet.push_back(tet7);
			box_Tet.push_back(tet7.bbox());
		}
		Tetrahedron tet8 (e, T3d_n[i].operator[](1),T3d_n[i].operator[](2),g);
		//if(!tet8.is_degenerate ()){
		if(abs(tet8.volume ())>eps){
			vect_Tet.push_back(tet8);
			box_Tet.push_back(tet8.bbox());
		}
		Tetrahedron tet9 (T3d_n[i].operator[](0),T3d_n[i].operator[](2),g,f);
		//if(!tet9.is_degenerate ()){
		if(abs(tet9.volume ())>eps){
			vect_Tet.push_back(tet9);
			box_Tet.push_back(tet9.bbox());
		}

		Tetrahedron tet10 (T3d_n[i].operator[](0),T3d_n[i].operator[](1),g,T3d_n[i].operator[](2));
		//if(!tet10.is_degenerate ()){
		if(abs(tet10.volume ())>eps){
			vect_Tet.push_back(tet10);
			box_Tet.push_back(tet10.bbox());
		}
		
		for(int iter=0; iter<box_cells.size(); iter++){
			
			double volume = 0.;
			
			if (CGAL::do_intersect(box_prismes[i], box_cells[iter]) ) {
			// test d'intersection des box_tetraedre avec les 8 cellules autour
			for(int it=0; it<box_Tet.size(); it++){
					if (CGAL::do_intersect(box_Tet[it],box_cells[iter]) ) {
						// test pour verifier si vect_Tet[it] est contenu dans la cellule "iter"
						if (inside_box(box_cells[iter], vect_Tet[it].operator[](0)) &&  inside_box(box_cells[iter], vect_Tet[it].operator[](1))
							&& inside_box(box_cells[iter], vect_Tet[it].operator[](2)) && inside_box(box_cells[iter], vect_Tet[it].operator[](3))){
						   volume += volume_tetra(vect_Tet[it]); 
						}
						else {
							//calcul volume intersection
							user_time.start();
							volume += (intersect_cube_tetrahedron(box_cells[iter], vect_Tet[it]) * sign(volume_tetra(vect_Tet[it])) );
							time+=user_time.time();
							user_time.reset();
						}
					} //if intersect Box_Tetra avec Box_Cell
				} // boucle sur tetra			
			}//if inter box_cell inter box_prisme
			
			if(std::abs(volume)>eps){
				c.delta_w[0] += volume*Cells[iter].rho0/(c.dx*c.dy*c.dz); 
				c.delta_w[1] += volume*Cells[iter].impx0/(c.dx*c.dy*c.dz);
				c.delta_w[2] += volume*Cells[iter].impy0/(c.dx*c.dy*c.dz); 
				c.delta_w[3] += volume*Cells[iter].impz0/(c.dx*c.dy*c.dz); 
				c.delta_w[4] += volume*Cells[iter].rhoE0/(c.dx*c.dy*c.dz);
				grille[in1][jn1][kn1] = c;
			}
			//vol_test += volume;
		 } // boucle sur les box_cells
		}//end else 
	} //end boucle sur les prismes
}	


void sous_maillage_faceTn_faceTn1_3d(Triangle_3& Tn, Triangles& tn, Triangle_3& Tn1, Triangles& tn1, Vector_3& N,Triangles& T3d_n,Triangles& T3d_n1){
	
	CGAL::Timer user_time, user_time2, user_time3, user_time4 ;
	double time=0.;	
	user_time.start();
	// transf barycentrique de tn 
	Triangles tn_n1(tn.size());
	for(int i=0; i<tn.size(); i++){		
		tn_n1[i] = tr(Tn, Tn1, tn[i]);
	}
	//cout << "Mapping Tn vers Tn1 pour une face time is: " << user_time.time() << " seconds." << endl;
	user_time.reset();
	
	user_time2.start();
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
	//cout << "Passage 3d-2d pour une face time is: " << user_time2.time() << " seconds." << endl;
	user_time2.reset();
	
	user_time3.start();
	// sous maillage triangulaire de l'interface
	CDT cdt;
	cdt.insert(Ap); cdt.insert(Bp); cdt.insert(Cp);
	sous_maillage_face(Tn1_2, Tn_n1_2, cdt);
	assert(cdt.is_valid());
	//	cout << "Sous-maillage en 2d pour une face time is: " << user_time3.time() << " seconds." << endl;
	user_time3.reset();
	
	Triangles_2 T2d; //recuperation faces du maillage Triangle_2
	
	
	for (CDT::Finite_faces_iterator fit=cdt.finite_faces_begin();
	fit!=cdt.finite_faces_end();++fit)
	{ 
		Point_2 s = fit->vertex(0)->point();
		Point_2 v = fit->vertex(1)->point();
		Point_2 r = fit->vertex(2)->point();
		if(Triangle_2(s,v,r).area()>eps){
			T2d.push_back(Triangle_2(s,v,r));
		}
	}
	
	user_time4.start();
	//transf des Triangle_2 en Triangle_3
	T3d_n1.resize(T2d.size());
	for(int i=0; i<T2d.size(); i++){
		Triangle_3 Tri = tr(Tn1,T2d[i]);
		Vector_3 vect0(Tri.operator[](0),Tri.operator[](1));
		Vector_3 vect1(Tri.operator[](0),Tri.operator[](2));
		Vector_3 normale = CGAL::cross_product(vect0,vect1);
		if (normale*N > 0.){ T3d_n1[i] =  Tri; }
		else{ T3d_n1[i] = Triangle_3(Tri.operator[](0),Tri.operator[](2),Tri.operator[](1));} 
	}
	//cout << "Passage 2d-3d pour une face time is: " << user_time4.time() << " seconds." << endl;
	user_time4.reset();
	
	//transf inverse 
	T3d_n.resize(T3d_n1.size());
	for(int i=0; i<T3d_n1.size(); i++){ 
		T3d_n[i] = tr(Tn1,Tn,T3d_n1[i]);
	}
	
}


void Grille::swap_3d(const double dt, Solide& S, int& n, int &n1, int& m){
	
	CGAL::Timer user_time, user_time2;
	double time_1=0., time_2=0.;
	for(int i=0;i<S.solide.size();i++){
		for (int j=0; j<S.solide[i].triangles.size(); j++){
			Triangles T3d_n,T3d_n1;
			user_time.start();
			n+= S.solide[i].Triangles_interface_prev[j].size();
			n1+=S.solide[i].Triangles_interface[j].size() ;
			sous_maillage_faceTn_faceTn1_3d(S.solide[i].triangles_prev[j], S.solide[i].Triangles_interface_prev[j] ,
										 S.solide[i].triangles[j], S.solide[i].Triangles_interface[j],
										 S.solide[i].normales[j], T3d_n,T3d_n1);
			m+= T3d_n.size();
			time_1+=CGAL::to_double(user_time.time());
		  user_time.reset();
		  user_time2.start();
		  swap_face(T3d_n,T3d_n1,dt);
			time_2+=CGAL::to_double(user_time2.time());
		  user_time2.reset();
		}
	}
	//cout << "Temps sous-maillage face: " << time_1 << " seconds." << endl;
	//cout << "Temps swap_modification_flux: " << time_2 << " seconds." << endl;
	
}


void Sous_Maillage_2d(const Triangles_2& Tn, const Triangles_2& Tn1, Triangles_2& tri2){
	
	std::vector<Bbox_2> boxesTn(Tn.size()), boxesTn1(Tn1.size()); //tres outil pour les intersections 
	//on associe a chaque triangle un Box (une boite contenant le triangle)
	for(int it=0; it< Tn.size(); it++){
		boxesTn[it] = Tn[it].bbox();
	}
	for(int iter=0; iter< Tn1.size(); iter++){
		boxesTn1[iter] = Tn1[iter].bbox();
	}
	Triangle_2 tri;
	std::vector<Point_2> vPoints; 
	int k=0;
	CGAL::Timer user_time;
	user_time.start();
	for(int i=0; i<boxesTn.size(); i++){ 
		for(int j=0; j<boxesTn1.size(); j++){
			if (CGAL::do_overlap(boxesTn[i],boxesTn1[j])){ //test d'intersection des Box 
				if (CGAL::do_intersect(Tn[i],Tn1[j])){ // test d'intersection des triangles contenues dans les Box
					CGAL::Object result = CGAL::intersection(Tn[i],Tn1[j]); //calcul d'intersection entre les deux triangles
					if(CGAL::assign(tri,result)){ tri2.push_back(tri); }
					else if(CGAL::assign(vPoints,result)){
						Triangulation_2 T;
						T.insert(vPoints.begin(), vPoints.end());
						if( (T.is_valid() ) && (T.dimension() == 2)){
							for (Triangulation_2::Finite_faces_iterator fit=T.finite_faces_begin(); fit!=T.finite_faces_end();++fit)
							{ 
								Point_2 s = fit->vertex(0)->point();
								Point_2 v = fit->vertex(1)->point();
								Point_2 r = fit->vertex(2)->point();
								if(Triangle_2(s,v,r).area()>eps){ tri2.push_back(Triangle_2(s,v,r));}
							}
						}
					}
				}
			}
		}
	}
}	
void sous_maillage_faceTn_faceTn1_2d(Triangle_3& Tn, Triangles& tn, Triangle_3& Tn1, Triangles& tn1, Vector_3& N,Triangles& T3d_n,Triangles& T3d_n1){
	
	CGAL::Timer user_time, user_time2, user_time3, user_time4 ;
	double time=0.;	
	user_time.start();
	// transf barycentrique de tn 
	Triangles tn_n1(tn.size());
	for(int i=0; i<tn.size(); i++){		
		tn_n1[i] = tr(Tn, Tn1, tn[i]);
	}
	user_time.reset();
	
	Triangles_2 Tn_n1_2_test(tn_n1.size());
	for(int i=0; i<tn_n1.size(); i++){
		Tn_n1_2_test[i] = tr(Tn1, tn_n1[i]);
	}
	
	Triangles_2 Tn1_2_test(tn1.size());
	for(int i=0; i<tn1.size(); i++){
		Tn1_2_test[i] =tr(Tn1, tn1[i]);
	}
	Triangles_2 tri2;
	Sous_Maillage_2d(Tn1_2_test, Tn_n1_2_test, tri2);
	//cout<<" size "<<tri2.size()<<endl;
	T3d_n1.resize(tri2.size());
	for(int i=0; i<T3d_n1.size(); i++){
		Triangle_3 Tri = tr(Tn1,tri2[i]);
		Vector_3 vect0(Tri.operator[](0),Tri.operator[](1));
		Vector_3 vect1(Tri.operator[](0),Tri.operator[](2));
		Vector_3 normale = CGAL::cross_product(vect0,vect1);
		if (normale*N > 0.){ T3d_n1[i] =  Tri; }
		else{ T3d_n1[i] = Triangle_3(Tri.operator[](0),Tri.operator[](2),Tri.operator[](1));} 
	}
	T3d_n.resize(T3d_n1.size());
	for(int i=0; i<T3d_n1.size(); i++){ 
		T3d_n[i] = tr(Tn1,Tn,T3d_n1[i]);
	}
	
}


void Grille::swap_2d(const double dt, Solide& S, int& n, int &n1, int& m){
	
	CGAL::Timer user_time, user_time2;
	double time_1=0., time_2=0.;

// 	double var=0.;
// 	vector< vector< vector<double > > > Test(Nx+2*marge, vector< vector<double> >(Ny+2*marge, vector<double>(Nz+2*marge,0.)) );
	for(int i=0;i<S.solide.size();i++){
		for (int j=0; j<S.solide[i].triangles.size(); j++){
			Triangles T3d_n,T3d_n1;
			user_time.start();
			n+= S.solide[i].Triangles_interface_prev[j].size();
			n1+=S.solide[i].Triangles_interface[j].size() ;
			sous_maillage_faceTn_faceTn1_2d(S.solide[i].triangles_prev[j], S.solide[i].Triangles_interface_prev[j] ,
										 S.solide[i].triangles[j], S.solide[i].Triangles_interface[j],
										 S.solide[i].normales[j], T3d_n,T3d_n1);
	     //cout << "Temps sous-maillage face: " << user_time.time() << " seconds." << endl;
			 m+= T3d_n.size();


		     //test 26 nov calcul aire faces
					 // cout<< " aire face is : "<< std::sqrt(CGAL::to_double(S.solide[i].triangles_prev[j].squared_area () ))<<endl;
					 double aire_f= std::sqrt(CGAL::to_double(S.solide[i].triangles_prev[j].squared_area () ));
					 double volume_f= volume_prisme(S.solide[i].triangles_prev[j],S.solide[i].triangles[j]);
						double a_n=0., a_n1=0., vol=0.;
						for(int iter=0; iter<T3d_n1.size(); iter++){
							a_n += std::sqrt(CGAL::to_double(T3d_n1[iter].squared_area () )); 
							//a_n1 += std::sqrt(CGAL::to_double(Test2[iter].squared_area () )); 
							vol += volume_prisme(T3d_n[iter], T3d_n1[iter]);
							}
							if (std::abs(vol)<eps){vol=0.;}
							if (std::abs(a_n -aire_f)>0.0001){cout<< " aire n is : "<< a_n<< "aire direct "<<aire_f<< endl;}
							if(std::abs(vol- volume_f)>0.00001){cout<<"volume balayee "<<vol<< "volume balayee n-n1 calcul direct "<<volume_f<<endl;}
				//fin test 26 nov	calcul aire faces	 OK
				
			time_1+=CGAL::to_double(user_time.time());
		  user_time.reset();
		  user_time2.start();
			
// 			//test 5 nov calcul volume
// 			cout<< "volume balayee n-n1 calcul direct: "<<volume_prisme(S.solide[i].triangles_prev[j],S.solide[i].triangles[j])<<endl;
// 			double vol=0., vol_tetra=0.;
// 			for(int iter=0; iter<T3d_n1.size(); iter++){
// 				//double vol=0., vol_tetra=0.;
// 				vol += volume_prisme(T3d_n[iter],T3d_n1[iter]);
// 				//test 6 nov
// 				Point_3 a = centroid(T3d_n[iter].operator[](1),T3d_n1[iter].operator[](1), T3d_n[iter].operator[](2),T3d_n1[iter].operator[](2));
// 				Point_3 b = centroid(T3d_n[iter].operator[](0),T3d_n1[iter].operator[](0), T3d_n[iter].operator[](2),T3d_n1[iter].operator[](2));
// 				Point_3 c = centroid(T3d_n[iter].operator[](0),T3d_n1[iter].operator[](0), T3d_n[iter].operator[](1),T3d_n1[iter].operator[](1));
// 				
// 				Tetrahedron tet0 (T3d_n[iter].operator[](0),T3d_n1[iter].operator[](0), c, b);
// 				vol_tetra+= volume_tetra(tet0);
// 				Tetrahedron tet1 (T3d_n[iter].operator[](1),T3d_n1[iter].operator[](1), a, c);
// 				vol_tetra+= volume_tetra(tet1);
// 				Tetrahedron tet2 (T3d_n[iter].operator[](2),T3d_n1[iter].operator[](2), b, a);
// 				vol_tetra+= volume_tetra(tet2);
// 				Tetrahedron tet3 (T3d_n[iter].operator[](0),T3d_n[iter].operator[](1), T3d_n[iter].operator[](2), c);
// 				vol_tetra+= volume_tetra(tet3);
// 				Tetrahedron tet4 (T3d_n[iter].operator[](0),c, T3d_n[iter].operator[](2), b);
// 				vol_tetra+= volume_tetra(tet4);
// 				Tetrahedron tet5 (T3d_n[iter].operator[](1),a, T3d_n[iter].operator[](2), c);
// 				vol_tetra+= volume_tetra(tet5);
// 				Tetrahedron tet6 (a,c,b, T3d_n[iter].operator[](2));
// 				vol_tetra+= volume_tetra(tet6);
// 				Tetrahedron tet7 (a,b,c, T3d_n1[iter].operator[](2));
// 				vol_tetra+= volume_tetra(tet7);
// 				Tetrahedron tet8 (a, T3d_n1[iter].operator[](1),T3d_n1[iter].operator[](2),c);
// 				vol_tetra+= volume_tetra(tet8);
// 				Tetrahedron tet9 (T3d_n1[iter].operator[](0),T3d_n1[iter].operator[](2),c,b);
// 				vol_tetra+= volume_tetra(tet9);
// 				Tetrahedron tet10 (T3d_n1[iter].operator[](0),T3d_n1[iter].operator[](1),c,T3d_n1[iter].operator[](2));
// 				vol_tetra+= volume_tetra(tet10);
// 				
// // 				if( std::abs(vol - vol_tetra)>eps ){ 
// // 					 cout<<"volume prisme "<< vol<<" i= "<< i<<endl; 
// // 					 cout<<"volume  tetra"<< vol_tetra<<endl;
// // 				}
// 			}
// 			if(std::abs(vol)<eps){vol=0.;}
// 			cout<< "volume balayee n-n1:               "<<vol<<endl;
// 			//fin test 5 nov	calcul volume	 OK
// 
//        if(std::abs(vol_tetra)<eps){vol_tetra=0.;}
//        cout<< "volume tetra:                      "<<vol_tetra<<endl;
// 			 //fin test 6 nov ok
		   swap_face(T3d_n,T3d_n1,dt);
		  //cout << "Temps swap_modification_flux: " << user_time2.time() << " seconds." << endl;
			//cout << "nb sous triang face: " << T3d_n.size()<< endl;
			time_2+=CGAL::to_double(user_time2.time());
		  user_time2.reset();
		}
	}
	//cout << "Temps sous-maillage face: " << time_1 << " seconds." << endl;
	//cout << "Temps swap_modification_flux: " << time_2 << " seconds." << endl;

}

