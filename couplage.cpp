/*!
   \file couplage.cpp
   \brief D&eacute;finitions des fonctions sp&eacute;cifiques au couplage.
   \details Calcul des forces et moments fluides exerc&eacute;s sur le solide, modifications flux fluide, remplissage de cellules fant&ocirc;mes, calcul de la quantit&eacute; balay&eacute;e.
   \warning  <b> Proc&eacute;dures sp&eacute;cifique au couplage! </b>
 */


#include "fluide.hpp"
#include "intersections.cpp"

/*!
* \fn void Grille::Forces_fluide(Solide& S, const double dt)
* \brief Calcul des Forces (\a Particule.Ff) et Moments fluides (\a Particule.Mf) exerc&eacute;s sur le Solide.
* \details Soit \a f un morceau d'interface, la force de pression exerc&eacute;e par le fluide sur l'interface \a f est donn&eacute;e par :
\f{eqnarray*}{
	F_f  =  (- p^x \, A_f n^{x}_f, \,- p^y \,A_f n^{y}_f, \,- p^z \,A_f n^{z}_f )^t
\f} \n
o&ugrave;  \f$ A_f \f$ l'aire de l'interface f,  \f$ n_f \f$ la normale sortante &agrave; l'interface f, et \f$ p^x, p^y, p^z \f$ les pressions efficaces selon les directions x, y et z pendant le pas de temps (\a Cellule.pdtx, 
\a  Cellule.pdty et \a Cellule.pdtz).
Le moment fluide exerc&eacute; sur f est donn&eacute; par :
\f{eqnarray*}{
	M_f = F_f  \wedge (X_f - X_I),
\f}
o&ugrave; \f$ X_f \f$ centre de l'interface f et \f$  X_I \f$ centre de la particule qui contient f.
Ces forces vont &ecirc;tre transmises au Solide comme des forces exerc&eacute;es par le fluide sur le solide pendant le pas de temps.
* \param S Solide
* \param dt pas de temps
* \warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
* \return void
*/
void Grille::Forces_fluide(Solide& S, const double dt){
	
	Vector_3 Ffluide(0.,0.,0.);
	//Mise &agrave; jour des Forces fluides et Moments fluides exerces sur le solide 
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
		Ffluide = Ffluide + S.solide[iter_s].Ff;
	} //fin boucle sur les particules
	cout<<"forces fluide "<<Ffluide<<endl;	
}	

/*!
* \fn void Grille::Modif_fnum(const double dt)
*  \brief Modification des flux fluide et bilan discret sur la cellule (m&eacute;thode de fronti&egrave;res immerg&eacute;es).
*  \details C'est &agrave; cette &eacute;tape que le fluide "voie" la pr&eacute;sence du solide. On calcule la valeur finale de l'&eacute;tat \f$ U^{n+1}_{i, j, k}\f$  dans la cellule en utilisant:
\f{eqnarray*}{
	\left( 1-  \Lambda_{i,j,k}^{n+1} \right) U^{n+1}_{i,j,k}   = \left( 1-  \Lambda_{i,j,k} ^{n+1}\right) U^n_{i,j,k}  + \Delta t \, \left(   \frac{(1-\lambda_{i-1/2,j,k}^{n+1} )}{\Delta x_{ i,j,k}} F_{i-1/2, j, k}^{n+1/2} -\frac{(1-\lambda_{i+1/2,j,k}^{n+1} )}{\Delta x_{ i,j,k}} F_{i+1/2, j, k}^{n+1/2}  + ...\right)	\f}
	\f{eqnarray*}{+  \frac{\Delta t}{V_{i,j,k}} \sum_{f \in C_{i,j,k}}{A_{f}} {\phi}_{f_{ i,j,k}}  +   \sum_{f \in C_{i,j,k}} \Delta U^{n, n+1}_{f_{ i,j,k}}  .
	\f}
	o&ugrave;	\f$ \Lambda_{i,j,k}^{n+1} \f$: \a Cellule.alpha (fraction occup&eacute;e par du solide dans la cellule(i,j,k)),\n
	\f$ \lambda_{i+1/2,j,k}^{n+1} \f$: \a Cellule.kappai (fraction occup&eacute;e par du solide sur les faces de la cellule(i,j,k)),\n
	\f$ F_{i+1/2, j, k}^{n+1/2} \f$:   \a  Cellule.fluxi (flux &agrave; droite dans la cellule(i,j,k)), \n
	\f$ V_{i,j,k} \f$: volume de la cellule(i,j,k), \n
	\f{eqnarray*}{ \phi_{f_{ i,j,k}} = (0, \Pi_x, \Pi_y,\Pi_z, \Pi_v)^t. 	\f}
	\f$ \Pi_x, \Pi_y,\Pi_z, \Pi_v\f$: \a Cellule.phi_x, \a Cellule.phi_y, Cellule.phi_z, \a Cellule.phi_v (les flux &agrave; la parois),\n
	\f$ \Delta U^{n, n+1}_{f_{ i,j,k}} 	\f$:  \a Cellule.delta_w (quantit&eacute;e balay&eacute;e).
*	\param dt pas de temps
*	\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
*	\return void
*/
void Grille::Modif_fnum(const double dt){
	
	Cellule c,ci,cj,ck; 
	double phi_x=0., phi_y=0., phi_z=0.;
	double vol=deltax*deltay*deltaz;
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
						c.flux_modif[l] /= (1.-c.alpha);
					}		
					//Mise a jour des valeurs dans les cellules
					c.rho = c.rho0  +  c.flux_modif[0];
					c.impx = c.impx0 + c.flux_modif[1];
					c.impy = c.impy0 + c.flux_modif[2];
					c.impz = c.impz0 + c.flux_modif[3];
					c.rhoE = c.rhoE0 + c.flux_modif[4];
					c.u = c.impx/c.rho;
					c.v = c.impy/c.rho;
					c.w = c.impz/c.rho;
					c.p = (gam-1.)*(c.rhoE-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2. - c.rho*c.w*c.w/2.);
					phi_x+=c.phi_x*vol/dt; phi_y+=c.phi_y*vol/dt; phi_z+=c.phi_z*vol/dt;
					
				}
				grille[i][j][k] = c;      
			}
		}
	}
	Vector_3 Phi(phi_x, phi_y, phi_z); cout<<" Flux a la parois "<<Phi<<endl;
}

/*!
* \fn void Grille:: Mixage()
*  \brief M&eacute;lange conservatif de petites cellules coup&eacute;es.
*  \details On d&eacute;finit une petite cellule tel que \f$ alpha > epsa \f$ (\a Cellule.alpha fraction occup&eacute;e par du solide dans la cellule, et \a epsa: fraction de cellule coup&eacute;e d&eacute;finit dans parametres.hpp ). Afin ne pas modifier le pas de temps tout en garantissant la condition de CFL, les petites cellules sont fusionn&eacute;es avec leurs voisines. On note \a p une petite cellule et \a g une cellule voisine avec \a p compl&egrave;tement fluide~(\f$ alpha_g = 0 \f$ ). On d&eacute;finit les termes d'&eacute;change suivants :
\f{eqnarray*}{ E_{pg} = \frac{1}{2 - alpha_p} (U_g - U_{p}), \quad  E_{gp} = \frac{1-  alpha_p}{2 - alpha_p} (U_p - U_{g}) \f}
et on pose:
\f{eqnarray*}{
U_p = U_{p} + E_{pg}, \quad  \quad U_g = U_{g} + E_{gp} \f}
*	\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
*	\return void
*/
void Grille::Mixage(){
	
	Cellule cp, cg;
	
	bool test_fini = true;
	
	for(int i=marge;i<Nx+marge;i++){
		for(int j=marge;j<Ny+marge;j++){ 
			for(int k=marge;k<Nz+marge;k++){
				cp = grille[i][j][k];
				bool test=true;
				if( (cp.alpha>epsa ||cp.p <0. || cp.rho<0.) && abs(cp.alpha-1.)>eps){
					
					for(int ii=-1; ii<=1 && test; ii++){
						for(int jj=-1; jj<=1 && test; jj++){
							for(int kk=-1; kk<=1 && test; kk++){
								if (grille[i+ii][j+jj][k+kk].alpha <eps && grille[i+ii][j+jj][k+kk].p>0. && grille[i+ii][j+jj][k+kk].rho>0. && i+ii>=marge && i+ii<Nx+marge && j+jj>=marge && j+jj<Ny+marge && k+kk>=marge && k+kk<Nz+marge)
								{
									test=false;
									cg = grille[i+ii][j+jj][k+kk];
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
									
									grille[i][j][k] = cp;
									grille[i+ii][j+jj][k+kk] = cg;
									
// 									if( std::abs((1.-cp.alpha)*cp.Mrho + (1.-cg.alpha)*cg.Mrho)>eps){
// 										cout<<" rho p initial "<<temp_rhop<<" rho g initial "<<temp_rhog<<endl;
// 										cout<<" Mrho p "<<cp.Mrho<<" Mrho g  "<<cg.Mrho<<endl;
// 										std::cout<<"Probleme mixage"<< ((1-cp.alpha)*cp.Mrho + (1-cg.alpha)*cg.Mrho)<<std::endl; 
// 										std::cout<< "position du centre de la cellule : "<<cp.x << " "<<cp.y << " "<<cp.z << " "<< " rho "<<cp.rho  << " p "<<cp.p <<" alpha " << cp.alpha<<std::endl;
// 										std::cout<< "position du centre de la cellule de mixage: "<<cg.x << " "<<cg.y << " "<<cg.z << " "<< " rho "<<cg.rho  << " p "<<cg.p <<" alpha " << cg.alpha<<std::endl;
// 									}
								} //if cg.alpha==0
								
							}
						}
					}
					
					if(test){
						
						if (grille[i-2][j][k].alpha == 0. && grille[i-2][j][k].p>0. && grille[i-2][j][k].rho>0. && i-2>=marge)
						{
							cg = grille[i-2][j][k];
							
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
							
							grille[i][j][k] = cp;
							grille[i-2][j][k] = cg;
							test = false;
						}
						else if (grille[i+2][j][k].alpha == 0. && grille[i+2][j][k].p>0. && grille[i+2][j][k].rho>0. &&  i+2<Nx+marge)
						{
							cg = grille[i+2][j][k];
							
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
							
							grille[i][j][k] = cp;
							grille[i+2][j][k] = cg;
							test = false;
						}
						
						else if (grille[i][j-2][k].alpha == 0. && grille[i][j-2][k].p>0. && grille[i][j-2][k].rho>0. && j-2>=marge)
						{
							cg = grille[i][j-2][k];
							
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
							
							grille[i][j][k] = cp;
							grille[i][j-2][k] = cg;
							test = false;
							
						}
						else if (grille[i][j+2][k].alpha == 0. && grille[i][j+2][k].p>0. && grille[i][j+2][k].rho>0.&& j+2<Ny+marge)
						{
							cg = grille[i][j+2][k];
							
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
							
							grille[i][j][k] = cp;
							grille[i][j+2][k] = cg;
							test = false;
						}
						else if (grille[i][j][k-2].alpha == 0. && grille[i][j][k-2].p>0. && grille[i][j][k-2].rho>0.&& k-2>=marge)
						{
							cg = grille[i][j][k-2];
							
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
							
							grille[i][j][k] = cp;
							grille[i][j][k-2] = cg;
							test = false;
						}
						else if(grille[i][j][k+2].alpha == 0. && grille[i][j][k+2].p>0. && grille[i][j][k+2].rho>0. && k+2 < Nz+marge)
						{
							cg = grille[i][j][k+2];
							
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
					if(grille[i][j][k].p<0. || grille[i][j][k].rho<0.){
						test_fini = false;
					}
				}// 0.5<c.alpha<1.
			} //fin boucle sur la grille
		}
	}
	
	if(!test_fini){
		cout<<" non test_fini "<<endl;
		Mixage();
	}
}




/*!
* \fn void Grille::Fill_cel(Solide& S)
*  \brief Remplissage des cellules fictives (\a alpha = 1)
*  \details Afin de calculer les flux pr&egrave;s de l'interface solide-fluide, on d&eacute;finit dans les Cellule compl&egrave;tement occup&eacute;es par le Solide (\a alpha = 1) un &eacute;tat fictif qui sera &eacute;gal &agrave; la valeur de l'&eacute;tat de la maille miroir par rapport &agrave; l'interface. \n
Algo: on cherche l'interface la plus proche du centre de la cellule (boucle sur toutes les faces du Solide) et on calcule la projection du centre de la cellule par rapport &agrave; cette interface via la fonction <b> CGAL::projection(Point_3) </b>.
*	\param S  Solide 
*	\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
*	\return void
*/

void Grille::Fill_cel(Solide& S){
	
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
					//test 18 septembre 2013
					//if((std::abs(c.rho)>eps)){cout<<"cellule solide non vide "<< c.rho<<endl;}
					//fin test 18 septembre 2013
				  Point_3 center_cell(c.x, c.y, c.z);
				  int nbx=0, nby=0,nbz=0;
				  Point_3 projete(0.,0.,0.); //Projete sur la face la plus proche
					Vector_3 V_f(0.,0.,0.); //Vitesse de la paroi au point projete
				  double dist_min = 10000000.;
				  for(int iter=0; iter<nb_part; iter++){
						for(int it=0;it<S.solide[iter].triangles.size();it++){
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

/*!
* \fn double volume_prisme(const Triangle_3& T1,const Triangle_3& T2)
*  \brief Calcul de volume sign&eacute; d'un prisme.
* \details Le volume sign&eacute; du prisme ayant comme basses les triangles \f$ T1(A_1,B_1,C_1)  \f$ et \f$ T2(A_2,B_2,C_2) \f$ est donn&eacute; par: \n
\f{eqnarray*}{
	{\Vert A_1 B_1 C_1 A_2 B_2 C_2 \Vert}_P  = \frac{1}{36}  \left( 2 \vec{A_1 B_1} \wedge \vec{A_1 C_1} + 2 \vec{A_2 B_2} \wedge \vec{A_2 C_2} +  \vec{A_1 B_1}\wedge \vec{A_2 C_2}  +  \vec{A_2 B_2}\wedge \vec{A_1 C_1} \right) \cdot
	\f} 
	\f{eqnarray*}{  \left( \vec{A_1 A_2}  + \vec{B_1 B_2} + \vec{C_1 C_2}\right)\f}
	*	\param T1 Triangle_3  base du prisme
	*	\param T2 Triangle_3  base du prisme
	*	\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
	*	\return double
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

/*!
* \fn double volume_tetra(const Tetrahedron& Tet)
*  \brief Calcul de volume sign&eacute; d'un t&eacute;tra&egrave;dre.
* \details Le volume sign&eacute; du t&eacute;tra&egrave;dre \f$ T(A,B,C,D)  \f$ est donn&eacute; par: \n
\f{eqnarray*}{
	{\Vert A B C D \Vert}_{sign} = \frac{1}{6} \vec{A D} \cdot \left(  \vec{A B} \wedge \vec{A C}  \right) 
	\f}
	*	\param Tet Tetrahedron
	*	\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
	*	\return double
	*/
double volume_tetra(const Tetrahedron& Tet){
	
	double volume=0.;
	
	Vector_3 V = cross_product( Vector_3(Tet.operator[](0),Tet.operator[](1)), Vector_3(Tet.operator[](0),Tet.operator[](2)) );
	
	
	volume = CGAL::to_double( (Vector_3(Tet.operator[](0),Tet.operator[](3)) *V) );
	volume /= 6.;
	
	return volume;
}

/*!
* \fn Point_3 tr(Triangle_3 Tn, Triangle_3 Tn1, Point_3 Xn)
*  \brief Transformation barycentrique du point Xn.
* \details Soit Xn un point appartenant au triangle Tn(A_1,B_1,C_1). Le transform&eacute; barycentrique du Xn est le point Xn1(appartenant au triangle Tn1(A_2,B_2,C_2) ) donn&eacute; par: 
\f{eqnarray*}{
	\lambda = \frac{\Vert  \vec{C_1 Xn} \wedge \vec{C_1 B_1}  \Vert }{\Vert  \vec{C_1 A_1} \wedge \vec{C_1 B_1}  \Vert} 
	\f}
	\f{eqnarray*}{
		\mu = \frac{\Vert  \vec{C_1 Xn} \wedge \vec{C_1 A_1}  \Vert }{\Vert  \vec{C_1 B_1} \wedge \vec{C_1 A_1}  \Vert} 
		\f}
		\f{eqnarray*}{
			Xn1 = \lambda A_2 + \mu B_2 + (1-\lambda -\mu)C_2 
			\f}
	*	\param Tn Triangle_3
	*\param Tn1 Triangle_3 
	*\param Xn  Point_3
	*	\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
	*	\return Point_3
	*/
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

/*!
*\fn Triangle_3 tr(Triangle_3 Tn, Triangle_3 Tn1, Triangle_3 T)
*\brief Transformation barycentrique du Triangle T.
*\details Appel &agrave; la fonction tr(Triangle_3, Triangle_3, Point_3) pour chaque somment du triangle.
*\param Tn Triangle_3
*\param Tn1 Triangle_3 
*\param T  Triangle_3 
*\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
*\return Triangle_3 
*/
Triangle_3 tr(Triangle_3 Tn, Triangle_3 Tn1, Triangle_3 T){
	
	Point_3 s = tr( Tn,Tn1, T.operator[](0) );
	Point_3 r = tr( Tn,Tn1, T.operator[](1) );
	Point_3 v = tr( Tn,Tn1, T.operator[](2) );
	
	return Triangle_3(s, r, v);
}
//transformation inverse tr(Tn1,Tn, tr(Tn,Tn1,T))


/*!
* \fn Point_2 tr(Triangle_3 Tn1, Point_3 Xn)
*  \brief Transformation d'un Point_3 en Point_2.
* \details  Soit le triangle Tn1(A,B,C), le transform&eacute; du Point 3d Xn (appartenant &agrave; Tn1) en un point 2d est donn&eacute; par: 

\f{eqnarray*}{
	\lambda = \frac{\Vert  \vec{C Xn} \wedge \vec{C B}  \Vert }{\Vert  \vec{C A} \wedge \vec{C B}  \Vert} 
	\f}
	\f{eqnarray*}{
		\mu = \frac{\Vert  \vec{C A} \wedge \vec{C Xn}  \Vert }{\Vert  \vec{C A} \wedge \vec{C B}  \Vert} 
		\f}
		\f{eqnarray*}{
			X_{2d} = (\mu, (1-\lambda-\mu))
			\f}
*\param Tn1 Triangle_3 
*\param Xn  Point_3
*\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
*\return Point_2
*/
Point_2 tr(Triangle_3 Tn1, Point_3 Xn){
		
  double lambda = 0., mu = 0.;

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

/*!
*\fn Triangle_2 tr(Triangle_3 Tn1, Triangle_3 T)
*\brief Transformation d'un Triangle_3 en Triangle_2 
	*\details Appel &agrave; la fonction tr(Triangle_3, Point_3) pour chaque somment du triangle.
	*\param Tn1 Triangle_3 
	*\param T  Triangle_3 
	*\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
	*\return Triangle_2 
	*/
Triangle_2 tr(Triangle_3 Tn1, Triangle_3 T){
	
	Point_2 s = tr( Tn1, T.operator[](0) );
	Point_2 r = tr( Tn1, T.operator[](1) );
	Point_2 v = tr( Tn1, T.operator[](2) );
	
	return Triangle_2(s, r, v);
}


/*!
* \fn Point_3 tr(Triangle_3 Tn1, Point_2 Xn)
*  \brief Transformation d'un Point_2 en Point_3.
* \details  Soit le triangle Tn1(A,B,C), le transform&eacute; du Point 2d Xn(X,Y) en un point 3d (appartenant &agrave; Tn1) est donn&eacute; par: 

\f{eqnarray*}{
\lambda = 1- X - Y
\f}
\f{eqnarray*}{
\mu = X
\f}
\f{eqnarray*}{
	X_{3d} = \lambda A + \mu B + (1-\lambda -\mu)C 
	\f}
*\param Tn1 Triangle_3 
*\param Xn  Point_2
*\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
*\return Point_3
*/
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
/*!
*\fn Triangle_3 tr(Triangle_3 Tn1, Triangle_2 T)
*\brief  Transformation d'un Triangle_2 en Triangle_3
*\details Appel &agrave; la fonction tr(Triangle_3, Point_2) pour chaque somment du triangle.
*\param Tn1 Triangle_3 
*\param T  Triangle_2 
*\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
*\return Triangle_3 
*/
Triangle_3 tr(Triangle_3 Tn1, Triangle_2 T){
	
	Point_3 s = tr( Tn1, T.operator[](0) );
	Point_3 r = tr( Tn1, T.operator[](1) );
	Point_3 v = tr( Tn1, T.operator[](2) );
	
	return Triangle_3(s, r, v);
}


/* 
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
			
			// Transformation barycentrique du Triangle T 
			Triangle_3 tr2(Triangle_3 Tn, Triangle_3 Tn1, Triangle_3 T){
				
				Point_3 s = tr2( Tn,Tn1, T.operator[](0) );
				Point_3 r = tr2( Tn,Tn1, T.operator[](1) );
				Point_3 v = tr2( Tn,Tn1, T.operator[](2) );
				
				return Triangle_3(s, r, v);
			}*/

/*!
*\fn void Grille::cells_intersection_face(int& in,int& jn,int& kn,int& in1,int& jn1,int& kn1, std::vector<Bbox>& box_cells, std::vector<Cellule>& Cells)
*\brief Liste de cellules fluide intersect&eacute;es par un triangle d'interface entre t et t-dt
*\details Soit
\f$ C_1= grille[in][jn][kn]\f$ la cellule o&ugrave; se trouve le triangle d'interface au temps t et \f$ C_2= grille[in1][jn1][kn1]\f$ la cellule o&ugrave; se trouve le triangle d'interface au temps t-dt. \n
Par CFL, \f$ max(|in-in1 |, |jn-jn1 |, |kn-kn1 |)<=1 \f$ . 
On aura au plus 8 cellules intersect&eacute;es pas &quot;le prisme r&eacute;gl&eacute;&quot;(prisme ayant comme bases le triangle d'interface aux temps t-dt (\a Particule.triangles_prev)  et t (\a Particule.triangles) ):\n
-\f$ grille[in][jn][kn]\f$ \n
-\f$ grille[in1][jn][kn]\f$ \n
-\f$ grille[in][jn1][kn]\f$ \n
-\f$ grille[in][jn1][kn1]\f$ \n
-\f$ grille[in1][jn1][kn]\f$ \n
-\f$ grille[in1][jn][kn1]\f$ \n
-\f$ grille[in][jn1][kn1]\f$ \n
-\f$ grille[in1][jn1][kn1]\f$ \n

*\param (in,jn,kn) index d'une Cellule (cellule o&ugrave; se trouve le triangle d'interface au temps t-dt (\a Particule.triangles_prev))
*\param (in1,jn1,kn1) index d'une Cellule (cellule o&ugrave; se trouve le triangle d'interface au temps t (\a Particule.triangles))
*\param box_cells vecteur de Box 3d
*\param Cells vecteur de Cellule
*\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
*\return void
*/
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
/*!
*\fn void Grille::swap_face(Triangles& T3d_prev, Triangles& T3d_n, const double dt,  Particule & P)
*\brief Calcul de la quantit&eacute; balay&eacute;e par un morceau de parois entre t et t-dt. Calcul du flux &agrave; la parois.
*\details Algorithme:\n
- Construction du vecteur des Box contenant les prismes ayant comme bases T3d_prev et T3d_n. Boucle sur les prismes ainsi obtenus:
- On cherche l'index de la cellule qui contient T3d_prev(le triangle est enti&egrave;rement contenu dans une cellule) et celui de la cellule qui contient T3d_n(le triangle est enti&egrave;rement contenu dans une cellule).
- Si le prisme est contenu dans une seule cellule on calcule le volume du prisme via la fonction volume_prisme(const Triangle_3&,const Triangle_3&) et la quantit&eacute; balay&eacute;e par la face (\a Particule.triangles) est donn&eacute;e par : \f$  volume\_prisme*U^n/volume\_cellule \f$. Sinon,
 - On liste les cellules fluide intersect&eacute;es par le prisme via la fonction \a cells_intersection_face(int& ,int& ,int& ,int& ,int& ,int& , std::vector<Bbox>& s, std::vector<Cellule>& s).
 - On d&eacute;coupe le prisme en  t&eacute;tra&egrave;dres: soit  \f$ T1(A_1,B_1,C_1)\f$  et \f$ T2(A_2,B_2,C_2)\f$  les bases du prisme, on d&eacute;finit les points: \f$ A = \frac{1}{4}(B_1 + B_2 + C_1 +C_2) \f$ , \f$ B = \frac{1}{4}(A_1 + A_2 + C_1 +C_2) \f$ et \f$ C = \frac{1}{4}(A_1 + A_2 + B_1 + B_2 ) \f$. Les t&eacute;tra&egrave;dres d&eacute;coupant \f$ A_1,B_1,C_1 A_2,B_2,C_2 \f$ sont: \f$ A_1 A_2 C B \f$, \f$ B_1 B_2 A C \f$, \f$ C_1 C_2 B A \f$, \f$ A_1 C C_1 B \f$, \f$ B_1 A C_1 C \f$, \f$ A C B C_1 \f$, \f$ A B C C_2 \f$, \f$ A B_2 C_2 C \f$, \f$ A_1 B_1 C_1 C \f$, \f$ A_2 C_2 C B \f$, \f$ A_2 B_2 C C_2. \f$
 - Intersections de ces  t&eacute;tra&egrave;dres avec les cellules fluide intersect&eacute;es par le prisme via la fonction intersect_cube_tetrahedron(Bbox&, Tetrahedron&). La quantité balayée par la face est donnée par la somme des: \f$  volume\_{intersection\_cellule\_tetrahedre}*U^n/volume\_cellule. \f$ \n

Calcul du flux &agrave; la parois: soit \a f un morceau d'interface, le flux &agrave; la parois est donné par :
\f{eqnarray*}{
	\Phi_f  =  \left(0,  p^x \, A_f n^{x}_f, \, p^y \,A_f n^{y}_f, \, p^z \,A_f n^{z}_f, V_f \cdot \left( p^x \, A_f n^{x}_f,p^y \,A_f n^{y}_f,p^z \,A_f n^{z}_f \right)^t \right)^t
	\f} \n
	o&ugrave; \f$ A_f \f$ l'aire de l'interface f,  \f$ n_f \f$ la normale sortante &agrave; l'interface f, \f$ V_f \f$ la vitesse au centre de la parois calculée via la fonction \a vitesse_parois(Point_3& ) et \f$ p^x, p^y, p^z \f$ les pressions efficaces selon les directions x, y et z pendant le pas de temps (\a Cellule.pdtx, \a  Cellule.pdty et \a Cellule.pdtz).

*\param T3d_prev Triangles_3 (triangles d'interface au temps t: \a Particule.Triangles_interface)
*\param T3d_n    Triangles_3 (triangles d'interface au temps t-dt: \a Particule.Triangles_interface_prev)
*\param dt pas de temps
*\param P Particule 
*\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
*\return void
*/
void Grille::swap_face(Triangles& T3d_prev, Triangles& T3d_n, const double dt,  Particule & P, double & volume_test){
	
	//CGAL::Timer user_time, user_time2;
	//double time=0.;
	std::vector<Bbox> box_prismes(T3d_prev.size());
	for (int i=0; i< T3d_prev.size(); i++){
		Bbox box_triangles_prev = T3d_prev[i].bbox();
		Bbox box_triangles_n = T3d_n[i].bbox();
		box_prismes[i]= box_triangles_prev.operator+(box_triangles_n);
	} 
	
	for (int i=0; i< box_prismes.size(); i++){
		//double vol_test=0.; 
		int in=0, jn=0, kn=0, in1=0, jn1=0, kn1=0;
		bool interieur = true;
		Point_3 center_prev= centroid(T3d_prev[i].operator[](0),T3d_prev[i].operator[](1),T3d_prev[i].operator[](2));
		Point_3 center_n= centroid(T3d_n[i].operator[](0),T3d_n[i].operator[](1),T3d_n[i].operator[](2));
		in_cell(center_prev, in, jn, kn, interieur);
		in_cell(center_n, in1, jn1, kn1, interieur);
		
		//Cellule c_cur= grille[in1][jn1][kn1];
		if((std::abs(grille[in1][jn1][kn1].alpha -1.)<eps)  && (interieur==true)){
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
		double volume_cel = c.dx*c.dy*c.dz;  
		if ( (in==in1) && (jn==jn1) && (kn==kn1) && (interieur==true)){
			// le prisme est contenu dans une seule cellule 
			double volume_p=volume_prisme(T3d_prev[i],T3d_n[i]);
			//calcul du volume
			if( (std::abs(volume_p)>eps) && (std::abs(1.-c.alpha)>eps)){
				c.delta_w[0] += volume_p*c.rho0/volume_cel; 
				c.delta_w[1] += volume_p*c.impx0/volume_cel;
				c.delta_w[2] += volume_p*c.impy0/volume_cel; 
				c.delta_w[3] += volume_p*c.impz0/volume_cel; 
				c.delta_w[4] += volume_p*c.rhoE0/volume_cel;
			grille[in1][jn1][kn1] = c;
			}
			volume_test += volume_p;
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
							//user_time.start();
							volume += (intersect_cube_tetrahedron(box_cells[iter], vect_Tet[it]) * sign(volume_tetra(vect_Tet[it])) );
							//time+=user_time.time();
							//user_time.reset();
						}
					} //if intersect Box_Tetra avec Box_Cell
				} // boucle sur tetra			
			}//if inter box_cell inter box_prisme
			
			//if(std::abs(volume)>eps){ //11 cotobre 2013
				c.delta_w[0] += volume*Cells[iter].rho0/volume_cel; 
				c.delta_w[1] += volume*Cells[iter].impx0/volume_cel;
				c.delta_w[2] += volume*Cells[iter].impy0/volume_cel; 
				c.delta_w[3] += volume*Cells[iter].impz0/volume_cel; 
				c.delta_w[4] += volume*Cells[iter].rhoE0/volume_cel;
				grille[in1][jn1][kn1] = c;
				
			//}
			volume_test += volume;
		 } // boucle sur les box_cells
		}//end else 
		if (explicite){//explicit algo
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
//test 11 octobre 2013			
// 				if (abs(c.phi_x)<=eps) {c.phi_x = 0.;} 
// 				if (abs(c.phi_y)<=eps) {c.phi_y = 0.;} 
// 				if (abs(c.phi_z)<=eps) {c.phi_z = 0.;} 
// 				if (abs(c.phi_v)<=eps) {c.phi_v = 0.;} 
//fin test 11 octobre 2013			
				grille[in1][jn1][kn1] = c;
			}
		}//explicit algo
		else {//semi_implicit algo
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
//test 11 octobre 2013				
// 				if (abs(c.phi_x)<=eps) {c.phi_x = 0.;} 
// 				if (abs(c.phi_y)<=eps) {c.phi_y = 0.;} 
// 				if (abs(c.phi_z)<=eps) {c.phi_z = 0.;} 
// 				if (abs(c.phi_v)<=eps) {c.phi_v = 0.;} 
// test 11 octobre 2013
				grille[in1][jn1][kn1] = c;
			}
		} //semi_implicit algo
		
	} //end boucle sur les prismes
}	

/**
\fn void Sous_Maillage_2d(const Triangles_2& Tn, const Triangles_2& Tn1, Triangles_2& tri2)
\brief Construction sous-maillage 2d d'une face 2d du Solide.
\details Un d&eacute;coupage en triangles de la face aux temps t et t-dt tel que chaque triangle soit enti&egrave;rement contenu dans une cellule aux temps t et t-dt (pas n&eacute;cessairement la m&ecirc;me cellule).\n
Algorithme:\n
- On associe &agrave; chaque triangle un Box 2d (une bo&icirc;te contenant le triangle).  \n
- Boucle sur les Box 2d.
- Test d'intersection des Box via la fonction <b> CGAL::do_overlap(Bbox_2, Bbox_2)</b>. Si oui: \n
 - Test d'intersection des triangles contenues dans les Box via la fonction <b> CGAL::do_intersect(Triangle_2,Triangle_2)</b>. Si oui:  \n
  - Calcul d'intersections entre les deux triangles via la fonction <b>CGAL::intersection(Triangle_2,Triangle_2)</b>. \n
  - Si le r&eacute;sultat de l'intersection est un triangle on le rajoute dans le sous-maillage. Si le r&eacute;sultat de l'intersection est un  polygone (une liste de points) on le triangularise en utilisant la classe <b>CGAL::Triangulation</b>  et la fonction <b> CGAL::insert </b> de cette classe. (Si le r&eacute;sultat de l'intersection est un point ou un segment on ne fait rien car le volume balayé est nul.)
\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
\param Tn Triangles_2 
\param Tn1 Triangles_2
\param tri2 Triangles_2
\return void
*/
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
/**
\fn void sous_maillage_faceTn_faceTn1_2d(Triangle_3& Tn, Triangles& tn, Triangle_3& Tn1, Triangles& tn1, Vector_3& N,Triangles& T3d_n,Triangles& T3d_n1)
\brief D&eacute;coupage en triangles de la face aux temps t et t-dt
\details A partir de la position de l'interface au temps t (\a Tn) et au temps t-dt (\a Tn1) on va d&eacute;couper cette face en triangles enti&egrave;rement contenues dans une cellule aux temps t et t-dt (pas n&eacute;cessairement la m&ecirc;me cellule). \n
Algorithme:\n
- Transformation barycentrique de \a tn via la fonction \a tr(Triangle_3, Triangle_3, Triangle_3) et transformation des triangles r&eacute;sultants en triangles 2d via la fonction \a tr(Triangle_3, Triangle_3).
- Transformation de tn1 en triangles 2d via la fonction \a tr(Triangle_3, Triangle_3).
- Construction du sous-maillage 2d de la face via la fonction \a Sous_Maillage_2d(const Triangles_2&, const Triangles_2&, Triangles_2&).
- Transformation du sous-maillage 2d dans un sous-maillage 3d de la face via la fonction \a tr(Triangle_3, Triangle_2).

\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage!</b> 
\param Tn Triangle_3 : interface au temps t-dt (\a Particule.triangles_prev)
\param tn vecteur de Triangle_3 : triangulation de la face Tn (\a Particule.Triangles_interface_prev)
\param Tn1 Triangle_3 interface au temps t (\a Particule.triangles)
\param tn1 vecteur de Triangle_3 : triangulation de la face Tn1 (\a Particule.Triangles_interface)
\param N   Vector_3 : normale sortante au Tn1 (\a Particule.normales)
\param T3d_n vecteur de Triangle_3 : Sous-maillage triangulaire de la face \a Particule.Triangles_interface_prev au temps t-dt 
\param T3d_n1 vecteur de Triangle_3 : Sous-maillage triangulaire de la face \a Particule.Triangles_interface au temps t 
\return void
*/
void sous_maillage_faceTn_faceTn1_2d(Triangle_3& Tn, Triangles& tn, Triangle_3& Tn1, Triangles& tn1, Vector_3& N,Triangles& T3d_n,Triangles& T3d_n1){
	
	//CGAL::Timer user_time;
	//double time=0.;	
	//user_time.start();
	
	Triangles tn_n1(tn.size());
	for(int i=0; i<tn.size(); i++){
		tn_n1[i] = tr(Tn, Tn1, tn[i]); // transf barycentrique de tn 
	}
	//user_time.reset();
	
	Triangles_2 Tn_n1_2(tn_n1.size());
	for(int i=0; i<tn_n1.size(); i++){
		Tn_n1_2[i] = tr(Tn1, tn_n1[i]);
	}
	
	Triangles_2 Tn1_2(tn1.size());
	for(int i=0; i<tn1.size(); i++){
		Tn1_2[i] =tr(Tn1, tn1[i]);
	}
	Triangles_2 tri2;
	Sous_Maillage_2d(Tn1_2, Tn_n1_2, tri2);
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

/**
\fn void Grille::Swap_2d(const double dt, Solide& S)
\brief Calcul de la quantit&eacute; balay&eacute;e par le Solide entre t et t-dt.
\details Algorithme:\n
- Sous-découpage des faces du Solide (\a Particule.triangles et \a Particule.triangles_prev) en triangles contenues enti&egrave;rement dans une cellule au temps t et t-dt (pas nécessairement la m&ecirc;me cellule) via la fonction sous_maillage_faceTn_faceTn1_2d(Triangle_3&, Triangles&, Triangle_3&, Triangles&, Vector_3& ,Triangles& ,Triangles&).\n
- Calcul de la quantité balayée par les faces et du flux &agrave; la parois via la fonction swap_face(Triangles&, Triangles&, const double ,  Particule &).
\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b> 
\param S Solide
\param dt pas de temps
\return void
*/

void Grille::Swap_2d(const double dt, Solide& S){
	
	//CGAL::Timer user_time, user_time2;
	//double time_1=0., time_2=0.;
	double volume_test=0.;
	for(int i=0;i<S.solide.size();i++){
		for (int j=0; j<S.solide[i].triangles.size(); j++){
			if (S.solide[i].fluide[j]){
			Triangles T3d_n,T3d_n1;
			//user_time.start();
			sous_maillage_faceTn_faceTn1_2d(S.solide[i].triangles_prev[j], S.solide[i].Triangles_interface_prev[j] ,
										 S.solide[i].triangles[j], S.solide[i].Triangles_interface[j],
										 S.solide[i].normales[j], T3d_n,T3d_n1);
				
			//time_1+=CGAL::to_double(user_time.time());
		  //user_time.reset();
		  //user_time2.start();
			swap_face(T3d_n,T3d_n1,dt, S.solide[i],volume_test );
			//time_2+=CGAL::to_double(user_time2.time());
		  //user_time2.reset();
			}
		}
	}
	cout<<"volume balayee = "<< volume_test<<endl;
}

/**
\fn CDT Sous_Maillage_3d(Triangles_2& Tn1, Triangles_2& Tn_n1, CDT &cdt)
\brief Construction sous-maillage sous contrainte pour une face du Solide.
\details Un d&eacute;coupage en triangles de la face aux temps t et t-dt tel que chaque triangle soit enti&egrave;rement contenu dans une cellule aux temps t et t-dt (pas n&eacute;cessairement la m&ecirc;me cellule).
\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! Il  n'est pas n&eacute;cessaire de la re-coder car c'est une ancienne m&eacute;thode(plus co&ucirc;teuse que la nouvelle version)!!! </b> 
\return CDT
*/
CDT Sous_Maillage_3d(Triangles_2& Tn1, Triangles_2& Tn_n1, CDT &cdt){
	
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
/**
\fn void sous_maillage_faceTn_faceTn1_3d(Triangle_3& Tn, Triangles& tn, Triangle_3& Tn1, Triangles& tn1, Vector_3& N,Triangles& T3d_n,Triangles& T3d_n1)
\brief D&eacute;coupage en triangles de la face aux temps t et t-dt.
\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! Il  n'est pas n&eacute;cessaire de la re-coder car c'est une ancienne m&eacute;thode(plus co&ucirc;teuse que la nouvelle version)!!! </b> 
\return void
*/
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
	Sous_Maillage_3d(Tn1_2, Tn_n1_2, cdt);
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
/**
\fn void Grille::Swap_3d(const double dt, Solide& S)
\brief Calcul de la quantit&eacute; balay&eacute;e
\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! Il  n'est pas n&eacute;cessaire de la re-coder car c'est une ancienne m&eacute;thode(plus co&ucirc;teuse que la nouvelle version)!!! </b> 
\return void
*/
void Grille::Swap_3d(const double dt, Solide& S){
	//CGAL::Timer user_time, user_time2;
	//double time_1=0., time_2=0.;
	double volume_test=0.;
	for(int i=0;i<S.solide.size();i++){
		for (int j=0; j<S.solide[i].triangles.size(); j++){
			if (S.solide[i].fluide[j]){
				Triangles T3d_n,T3d_n1;
				//user_time.start();
				sous_maillage_faceTn_faceTn1_3d(S.solide[i].triangles_prev[j], S.solide[i].Triangles_interface_prev[j] ,
																				S.solide[i].triangles[j], S.solide[i].Triangles_interface[j],
																				S.solide[i].normales[j], T3d_n,T3d_n1);
				 //time_1+=CGAL::to_double(user_time.time());
				//user_time.reset();
				// user_time2.start();
				swap_face(T3d_n,T3d_n1,dt,S.solide[i], volume_test );
				//time_2+=CGAL::to_double(user_time2.time());
				// user_time2.reset();
			}
		}
	}
}
/*!
* \fn void Grille:: Mixage_cible()
*  \brief M&eacute;lange conservatif de petites cellules coup&eacute;es.
*  \details On d&eacute;finit une petite cellule tel que \f$ alpha > epsa \f$ (\a Cellule.alpha fraction occup&eacute;e par du solide dans la cellule, et \a epsa: fraction de cellule coup&eacute;e d&eacute;finit dans parametres.hpp ). Afin ne pas modifier le pas de temps tout en garantissant la condition de CFL, les petites cellules sont fusionn&eacute;es avec leurs voisines.
	*	\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
	*	\return void
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
				Cellule cp = grille[i][j][k];
				int ii=i, jj=j, kk=k;
				if((cp.alpha>epsa || cp.p<0. || cp.rho<0.) && abs(cp.alpha-1.)>eps){
					Cellule cg = cible(grille[i][j][k], i, j,k, ii,jj,kk);
					
					cg.cible_alpha += (1.-cp.alpha);
					cg.cible_rho  += (1.-cp.alpha)*cp.rho;
					cg.cible_impx += (1.-cp.alpha)*cp.impx;
					cg.cible_impy += (1.-cp.alpha)*cp.impy;
					cg.cible_impz += (1.-cp.alpha)*cp.impz;
					cg.cible_rhoE += (1.-cp.alpha)*cp.rhoE;
					
					cp.cible_i= ii;
					cp.cible_j = jj;
					cp.cible_k = kk;
					
					grille[i][j][k] = cp;
					grille[ii][jj][kk] = cg;

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
				Cellule cp = grille[i][j][k];

				if(std::abs(cp.cible_alpha)>0.){
					cp.rho = ((1.-cp.alpha)*cp.rho + cp.cible_rho)/((1.-cp.alpha) + cp.cible_alpha);
					cp.impx = ((1.-cp.alpha)*cp.impx + cp.cible_impx)/((1.-cp.alpha) + cp.cible_alpha);
					cp.impy = ((1.-cp.alpha)*cp.impy + cp.cible_impy)/((1.-cp.alpha) + cp.cible_alpha);
					cp.impz = ((1.-cp.alpha)*cp.impz + cp.cible_impz)/((1.-cp.alpha) + cp.cible_alpha);
					cp.rhoE = ((1.-cp.alpha)*cp.rhoE + cp.cible_rhoE)/((1.-cp.alpha) + cp.cible_alpha);
					cp.u = cp.impx/cp.rho;
					cp.v = cp.impy/cp.rho;
					cp.w = cp.impz/cp.rho;
					cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
					grille[i][j][k] = cp;
				}
			}
		}
	}
	bool test_fini = true;
	for(int i=marge;i<Nx+marge;i++){
		for(int j=marge;j<Ny+marge;j++){ 
			for(int k=marge;k<Nz+marge;k++){
				Cellule cp = grille[i][j][k];
				Cellule cible=grille[cp.cible_i][cp.cible_j][cp.cible_k];
					cp.rho = cible.rho;
					cp.impx = cible.impx;
					cp.impy = cible.impy;
					cp.impz = cible.impz;
					cp.rhoE = cible.rhoE;
					cp.u = cp.impx/cp.rho;
					cp.v = cp.impy/cp.rho;
					cp.w = cp.impz/cp.rho;
					cp.p = (gam-1.)*(cp.rhoE-cp.rho*cp.u*cp.u/2.-cp.rho*cp.v*cp.v/2. - cp.rho*cp.w*cp.w/2.);
					grille[i][j][k] = cp;
					
					if(grille[i][j][k].p<0. || grille[i][j][k].rho<0.){
						test_fini = false;
					}
			}
		}
	}
	
	if(!test_fini){
		cout<<" non test_fini "<<endl;
		Mixage_cible();
	}
}
