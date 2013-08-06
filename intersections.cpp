/*!
 *  \file intersections.cpp
 *  \brief Intersection de la grille fluide avec le solide.
 \details Impl&eacute;mentation de la fonctions intersect_cube_tetrahedron(Bbox& cube, Tetrahedron& Tet) et Grille::Parois(Solide& S,double dt). 
   \warning  <b> Proc&eacute;dures sp&eacute;cifiques au couplage! </b>
 */

#include "intersections.hpp"
#include "fluide.hpp"
#include "solide.hpp"
#include "solide.cpp"

using std:: cout;
using std:: endl;


/*!
* \fn void Grille::Parois(Solide& S,double dt) 
*\brief Intersection de la grille fluide avec le solide.
\details Intersection de la Grille fluide avec le Solide et calcul des diff&eacute;rents quantit&eacute;s d'int&eacute;r&ecirc;t: occupation du Solide dans la Cellule : \a Cellule.alpha, occupation des faces de la cellule par le solide: \a Cellule.kappai, \a Cellule.kappaj et \a Cellule.kappak. Remplissage des vecteurs : \n
- \a Particule.Points_interface: points d'intersections de la cellule avec les faces triangulaires du solide; \n
- \a Particule.Triangles_interface: d&eacute;coupage des faces des particules en morceaux triangulaires d'interface contenus dans une seule cellule de la grille fluide; \n
- \a Particule.Position_Triangles_interface: index de la Cellule fluide contenant \a Particule.Triangles_interface. \n

Algorithme: \n
- Construction du vecteur \a box_grille contenant les cellules cubiques de la grille fluide sous la forme de Box 3d (\a Bbox). Une Box repr&eacute;sente une bo&icirc;te rectangulaire. Cette fa&ccedil;on de voir les cellules fluide permet de faire appel aux fonctions membres de la classe <b> CGAL::Bbox_3  </b>. 
- Construction du vecteur \a solide des Bbox associes aux \a Particule. 
- Boucle sur \a box_grille.
- Boucle sur \a solide.
- Test d'intersection entre \a box_grille et \a solide  via la fonction \b CGAL::do_intersect(Bbox, Bbox). Sinon, il n'y a pas d'intersection et on passe directement &agrave; la cellule suivante. Si oui:
- On test si \a box_grille est contenu dans \a solide via la fonction box_inside_convex_polygon(const Particule&, const Bbox&). Si oui,  l'intersection est \a box_grille, sinon:
- Boucle sur les faces triangulaires du Solide.
- Recherche des sommets de faces triangulaires du Solide contenues dans \a box_grille via la fonction inside_box(const Bbox&, const Point_3&). 
- Test d'intersections entre \a box_grille est les faces triangulaires du Solide via la fonction \b CGAL::do_intersect(Bbox,Triangle_3). Si non, il n'y a rien a faire, si oui:
- Triangulation des faces du \a box_grille via la fonction triang_cellule(const Bbox&, Triangles&) .
- Boucle sur les faces triangulaires du \a box_grille.
- Recherche des sommets de faces du \a box_grille contenues dans Solide via la fonction inside_convex_polygon(const Particule&, const Point_3&). 
- Test d'intersection entre les faces du \a box_grille et les faces triangulaires du solide (\a Particule.triangle) via la fonction \b CGAL::do_intersect(Triangle_3, Triangle_3). Si non, il n'y a rien a faire, si oui:
- Intersection entre les faces triangulaires du \a box_grille et les faces triangulaires du Solide via la fonction
\b  CGAL::intersection(Triangle_3, Triangle_3).\n




\remark Le vecteur \a Particule.Points_interface, contenant les points d'intersections entre les cellules fluides et les faces triangulaires des particules (\a Particule.triangles) est rempli fur et &agrave; mesure dans l'algorithme de recherche des intersections. \n
On exploite les r&eacute;sultats des intersections afin de calculer les diff&eacute;rents quantit&eacute;s d'int&eacute;r&ecirc;t: \a Cellule.alpha, \a Cellule.kappai, \a Cellule.kappaj et \a Cellule.kappak. \n
Les fonctions utilis&eacute;es dans ce calcul sont: 
<b> CGAL::Triangulation(vector<Point_3>) </b> et <b> CGAL::convex_hull_3(vector<Point_3>, Polyhedron_3) </b>, <b> CGAL::tetrahedron.volume() </b> et <b> CGAL::Triangle_3.squared_area()</b>. \n
La derni&egrave;re &eacute;tape consiste dans le remplissage du vecteur \a Particule.Triangles_interface qui contient
le d&eacute;coupage des faces des particules en morceaux triangulaires d'interface contenus dans une seule
cellule de la grille fluide et, \a Particule.Position_Triangles_interface contenant l'index de la Cellule fluide contenant \a Particule.Triangles_interface. 
Le d&eacute;coupage des faces des particules en morceaux triangulaires d'interface se fait &agrave; partir du vecteur
\a Particule.Points_interface qui contient pour chaque face triangulaire des particules la liste de points d'intersections
de celles-ci avec la grille fluide. On fait simplement une triangulation contenant ces points via la fonction
<b> CGAL::Triangulation(vector<Point_3>)</b>.
*\param S Solide
*\param dt pas de temps
*\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
*\return void 
*/

void Grille::Parois(Solide& S,double dt) {
	
	const double eps_relat = numeric_limits<double>::epsilon( );
	
	std::vector<Bbox> box_grille;
	const int nx_m=Nx+2*marge;
	const int ny_m=Ny+2*marge;
	const int nz_m=Nz+2*marge;
	const int Ns = (nx_m+1)*(ny_m+1)*(nz_m+1);
	const double volume_cel = deltax*deltay*deltaz;
	
	vector<vector<double> > Sommets(Ns, vector<double>(3,0.));
	int l=0;
	for(int i=0;i<nx_m+1;i++){
		for(int j=0;j<ny_m+1;j++){
			for(int k=0;k<nz_m+1;k++){
				Sommets[l][0] = (i-marge)*deltax;
				Sommets[l][1] = (j-marge)*deltay;
				Sommets[l][2] = (k-marge)*deltaz;
				l++;        
			}
		}
	}
	
	
	double x_min=0.,y_min=0.,z_min=0.,x_max=0.,y_max=0.,z_max=0.;
	
	for(int i=0; i<nx_m; i++){
		for(int j=0; j<ny_m; j++){ 
			for(int k=0; k<nz_m; k++){ 
				x_min = Sommets[k+j*(nz_m+1)+i*(nz_m+1)*(ny_m+1)][0];
				y_min = Sommets[k+j*(nz_m+1)+i*(nz_m+1)*(ny_m+1)][1];
				z_min = Sommets[k+j*(nz_m+1)+i*(nz_m+1)*(ny_m+1)][2];
				
				x_max = Sommets[(k+1)+(j+1)*(nz_m+1)+(i+1)*(nz_m+1)*(ny_m+1)][0];
				y_max = Sommets[(k+1)+(j+1)*(nz_m+1)+(i+1)*(nz_m+1)*(ny_m+1)][1];
				z_max = Sommets[(k+1)+(j+1)*(nz_m+1)+(i+1)*(nz_m+1)*(ny_m+1)][2];
				
				box_grille.push_back(Bbox(x_min,y_min,z_min,x_max,y_max,z_max));
				
			}
		}
	}
	
	
	
	int nb_particules = S.size();
	int taille = box_grille.size();
	cout<<"the grille size is : "<<taille<<endl;
	cout<<"nombres de particules : "<<nb_particules<<endl;
	
	
	const double eps_box = 0.1;
	std::vector<Bbox> solide(nb_particules);
	for(int it=0; it<nb_particules; it++){
		solide[it] = Bbox(S.solide[it].min_x-eps_box,S.solide[it].min_y-eps_box, S.solide[it].min_z-eps_box, S.solide[it].max_x+eps_box, S.solide[it].max_y+eps_box, S.solide[it].max_z+eps_box);
	}
	
	
	double volume_s=0.;
	
	Cellule cel;
	int i=0;
	CGAL::Timer user_time, user_time2;
	user_time.start();
	double time=0.;
	for (int a=0; a< nx_m; a++){
		for (int b=0; b< ny_m; b++){
			for (int c=0; c< nz_m; c++){
				
				std::vector< std::vector< std::vector<Point_3> > > Points_interface;
				Points_interface.resize(nb_particules);
				for(int count=0; count<nb_particules; count++){
					Points_interface[count].resize(S.solide[count].triangles.size(), std::vector<Point_3>(0));
				}
				cel = grille[a][b][c]; 
				cel.alpha = 0.; cel.kappai = 0.; cel.kappaj = 0.; cel.kappak = 0.;
				cel.phi_x = 0.; cel.phi_y = 0.; cel.phi_z = 0.; cel.phi_v = 0.; 
				cel.delta_w[0]= 0.; cel.delta_w[1]=0.; cel.delta_w[2]=0.; cel.delta_w[3]=0.; cel.delta_w[4] = 0.;
				std::vector<Point_3> Points_poly; 
				double alpha = 0.0;
				std::vector<double>  kappa(6,0.0);
				bool intersection = false;
				bool exterieur = true;
				bool box_in_solide = false;
				bool point_in_solide = false;
				bool point_in_cell = false;
				Triangles trianglesB;
				for(int iter_s=0; iter_s<nb_particules; iter_s++){ //boucle sur les particules 
					if (CGAL::do_intersect(box_grille[i],solide[iter_s]) ) {	
						intersection = true;
						box_in_solide = box_inside_convex_polygon(S.solide[iter_s],box_grille[i]);
						if(box_in_solide){
							exterieur = false;
							cel.alpha = 1.;
							cel.kappai = 1.;
							cel.kappaj = 1.;
							cel.kappak = 1.;
							//test 3 janvier 2013
							if(a>0) {grille[a-1][b][c].kappai = 1.;}
							if(b>0) {grille[a][b-1][c].kappaj = 1.;}
							if(c>0) {grille[a][b][c-1].kappak = 1.;}
							//fin test 3 janvier 2013
							volume_s += volume_cel;
							box_in_solide = false;
						}
						else 
						{    
							triang_cellule(box_grille[i] , trianglesB); 
							
							for ( int j = 0; j < S.solide[iter_s].triangles.size(); j++){ 
								
								//test if point is in cell_box
								point_in_cell = inside_box(box_grille[i], S.solide[iter_s].triangles[j].operator[](0));
								if(point_in_cell) {Points_poly.push_back(S.solide[iter_s].triangles[j].operator[](0)); point_in_cell=false;
								Points_interface[iter_s][j].push_back(S.solide[iter_s].triangles[j].operator[](0));
								}
								
								point_in_cell = inside_box(box_grille[i], S.solide[iter_s].triangles[j].operator[](1));
								if(point_in_cell) {Points_poly.push_back(S.solide[iter_s].triangles[j].operator[](1)); point_in_cell=false;
								Points_interface[iter_s][j].push_back(S.solide[iter_s].triangles[j].operator[](1));
								}
								
								point_in_cell = inside_box(box_grille[i], S.solide[iter_s].triangles[j].operator[](2));
								if(point_in_cell) {Points_poly.push_back(S.solide[iter_s].triangles[j].operator[](2)); point_in_cell=false;
								Points_interface[iter_s][j].push_back(S.solide[iter_s].triangles[j].operator[](2));
								}
								
								
								if (CGAL::do_intersect(box_grille[i],S.solide[iter_s].triangles[j]) ) {
									for ( int k = 0; k < trianglesB.size(); k++){
										//test if point is in solide
										point_in_solide = inside_convex_polygon(S.solide[iter_s],trianglesB[k].operator[](0));
										if(point_in_solide) {Points_poly.push_back(trianglesB[k].operator[](0)); point_in_solide=false;
										}
										
										point_in_solide = inside_convex_polygon(S.solide[iter_s],trianglesB[k].operator[](1));
										if(point_in_solide) {Points_poly.push_back(trianglesB[k].operator[](1)); point_in_solide=false;
										}
										
										point_in_solide = inside_convex_polygon(S.solide[iter_s],trianglesB[k].operator[](2));
										if(point_in_solide) {Points_poly.push_back(trianglesB[k].operator[](2)); point_in_solide=false;
										}
										
										if (CGAL::do_intersect(S.solide[iter_s].triangles[j],trianglesB[k]) ) {
											
											Triangle_3 t;
											Point_3 P;
											Segment_3 seg;
											std::vector<Point_3> vPoints; 
											CGAL::Object result = CGAL::intersection(S.solide[iter_s].triangles[j], trianglesB[k]);
											
											if(CGAL::assign(P,result)){
												Points_poly.push_back(P);
												Points_interface[iter_s][j].push_back(P);
											}
											else if(CGAL::assign(seg,result)){
												Points_poly.push_back(seg.operator[](0));
												Points_poly.push_back(seg.operator[](1));
												Points_interface[iter_s][j].push_back(seg.operator[](0));
												Points_interface[iter_s][j].push_back(seg.operator[](1));
											}
											else if(CGAL::assign(t,result)){
												Points_poly.push_back(t.operator[](0));
												Points_poly.push_back(t.operator[](1));
												Points_poly.push_back(t.operator[](2));
												Points_interface[iter_s][j].push_back(t.operator[](0));
												Points_interface[iter_s][j].push_back(t.operator[](1));
												Points_interface[iter_s][j].push_back(t.operator[](2));
												
											}
											else if(CGAL::assign(vPoints,result)){ 
												for(int l= 0; l<vPoints.size(); l++)
												{
													Points_poly.push_back(vPoints[l]);
													Points_interface[iter_s][j].push_back(vPoints[l]);
												}
												
											}
											else {cout<<"Intersection type: ? in intersect particule with grille"<< trianglesB[k]<<endl;}
										}
										//end intr t1 et t2
									} //end boucle sur trianglesB
									
								} // if intersection cellule[i] avec triangle[j]	 
								
							} // end boucle sur triangles
						} //else 
						
					} // if inter grille[i] et solide[iter_s]  
				} //fin boucle sur les particules
				user_time2.start();
				//traitement calcul de alpha et kappa pour la cellule=grille[i]!!!!!!
				if(intersection && exterieur){
					Triangulation T(Points_poly.begin(), Points_poly.end());
					
					if (T.dimension() == 3){
						Polyhedron_3 poly;
						CGAL::convex_hull_3(T.points_begin(), T.points_end(), poly);
						//int l=0, n=0;
						Finite_cells_iterator cit;
						for (cit = T.finite_cells_begin(); cit!= T.finite_cells_end(); cit++){
							alpha+= CGAL::to_double(T.tetrahedron( cit).volume());
						}
						
						Facet_iterator fiter;
						for (fiter = poly.facets_begin(); fiter!= poly.facets_end(); fiter++){
							
							Triangle_3 K((*fiter).halfedge()->vertex()->point(),(*fiter).halfedge()->next()->vertex()->point(),
													 (*fiter).halfedge()->opposite()->vertex()->point());
													 
													 if (abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](0).operator[](2))<=eps_relat && abs(trianglesB[0].operator[](0).operator[](2) - K.operator[](1).operator[](2))<=eps_relat && abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](2).operator[](2))<=eps_relat )
													 { 
														 kappa[0] +=sqrt(CGAL::to_double(K.squared_area()));
													 }
													 
													 else if (abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](0).operator[](2))<=eps_relat && abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](1).operator[](2))<=eps_relat && abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](2).operator[](2))<=eps_relat )
													 { 
														 kappa[1] +=sqrt(CGAL::to_double(K.squared_area()));
													 }
													 
													 else if (abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](0).operator[](0))<=eps_relat && abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](1).operator[](0))<=eps_relat && abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](2).operator[](0))<=eps_relat)
													 { 
														 kappa[2] +=sqrt(CGAL::to_double(K.squared_area()));
													 }
													 
													 else if (abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](0).operator[](0))<=eps_relat && abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](1).operator[](0))<=eps_relat && abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](2).operator[](0))<=eps_relat)
													 { 
														 kappa[3] +=sqrt(CGAL::to_double(K.squared_area()));
													 }
													 
													 else if (abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](0).operator[](1))<=eps_relat && abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](1).operator[](1))<=eps_relat && abs(trianglesB[8].operator[](0).operator[](1) - K.operator[](2).operator[](1))<=eps_relat)
													 { 
														 kappa[4] +=sqrt(CGAL::to_double(K.squared_area()));
													 }
													 
													 else if (abs(trianglesB[10].operator[](0).operator[](1) -  K.operator[](0).operator[](1))<=eps_relat && abs(trianglesB[10].operator[](0).operator[](1) -  K.operator[](1).operator[](1))<=eps_relat && abs(trianglesB[10].operator[](0).operator[](1) - K.operator[](2).operator[](1))<=eps_relat)
													 { 
														 kappa[5] +=sqrt(CGAL::to_double(K.squared_area()));
													 }
													 
													 else{
														 
													 } //calcul des aires parietales
													 
						}
					}
					
					if (T.dimension() == 2){
						//int l=0, n=0;
						Finite_faces_iterator it;
						for (it = T.finite_facets_begin(); it != T.finite_facets_end(); it++){
							Triangle_3 K= T.triangle(*it);
							
							if (abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](0).operator[](2))<=eps_relat && abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](1).operator[](2))<=eps_relat && abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](2).operator[](2))<=eps_relat)
							{ 
								kappa[0] +=sqrt(CGAL::to_double(K.squared_area()));
								
							}
							
							else if (abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](0).operator[](2))<=eps_relat && abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](1).operator[](2))<=eps_relat && abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](2).operator[](2))<=eps_relat)
							{ 
								kappa[1] +=sqrt(CGAL::to_double(K.squared_area()));
								
							}
							
							else if (abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](0).operator[](0))<=eps_relat && abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](1).operator[](0))<=eps_relat && abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](2).operator[](0))<=eps_relat)
							{ 
								kappa[2] +=sqrt(CGAL::to_double(K.squared_area()));
								
							}
							
							else if (abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](0).operator[](0))<=eps_relat && abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](1).operator[](0))<=eps_relat && abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](2).operator[](0))<=eps_relat)
							{ 
								kappa[3] +=sqrt(CGAL::to_double(K.squared_area()));
								
							}
							
							else if (abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](0).operator[](1))<=eps_relat && abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](1).operator[](1))<=eps_relat && abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](2).operator[](1))<=eps_relat)
							{ 
								kappa[4] +=sqrt(CGAL::to_double(K.squared_area()));
								
							}
							
							else 
							{ 
								kappa[5] +=sqrt(CGAL::to_double(K.squared_area()));
								
							}
						}
						
					}
					
					cel.alpha = alpha/volume_cel;
					cel.kappai = kappa[3]/(deltay * deltaz);
					cel.kappaj = kappa[4]/(deltax * deltaz);
					cel.kappak = kappa[1]/(deltax * deltay);

					//test 31 janvier 2013
					if(a>0){grille[a-1][b][c].kappai = kappa[2]/(deltay * deltaz);}
					if(b>0){grille[a][b-1][c].kappaj = kappa[5]/(deltax * deltaz);}
					if(c>0){grille[a][b][c-1].kappak = kappa[0]/(deltax * deltay);}
					//fin test 31 janvier 2013
					
					volume_s +=alpha;
					
				}
				grille[a][b][c] = cel;
				//test 3 janvier 2013
				if(std::abs(grille[a][b][c].alpha -1.) <1.e-10) {
					grille[a][b][c].alpha = 1.;
					grille[a][b][c].kappai = 1.;
					grille[a][b][c].kappaj = 1.;
					grille[a][b][c].kappak = 1.;
					if(a>0){grille[a-1][b][c].kappai = 1.;}
					if(b>0){grille[a][b-1][c].kappaj = 1.;}
					if(c>0){grille[a][b][c-1].kappak = 1.;}
				}
				//fin test 31 janvier 2013
				i++;
				time+= user_time2.time();
				user_time2.reset();
				
				//triangularisation de l'interface face par face
				Finite_faces_iterator iter;
				for(int count=0; count<nb_particules;count++){
					for(int it=0; it<S.solide[count].triangles.size(); it++){
						Triangulation T(Points_interface[count][it].begin(), Points_interface[count][it].end());
						assert(T.is_valid());
						if(T.dimension()==2){
							for (iter = T.finite_facets_begin(); iter != T.finite_facets_end(); iter++){
								Triangle_3 Tri= T.triangle(*iter);
								if(std::sqrt(CGAL::to_double(Tri.squared_area())) >eps){
									Vector_3 vect0(Tri.operator[](0),Tri.operator[](1));
									Vector_3 vect1(Tri.operator[](0),Tri.operator[](2));
									Vector_3 normale = CGAL::cross_product(vect0,vect1);
									if (normale*S.solide[count].normales[it] > 0.){
										S.solide[count].Triangles_interface[it].push_back(Tri);
									}
									else {
										S.solide[count].Triangles_interface[it].push_back(Triangle_3(Tri.operator[](0),Tri.operator[](2),Tri.operator[](1)));
									}
									std::vector<int> poz(3); poz[0]= a; poz[1] = b; poz[2] = c;
									S.solide[count].Position_Triangles_interface[it].push_back(poz);
								}
							}
						}
					} //fin boucle sur Triangles
				} //fin boucle sur les particules
			} //fin boucle sur grille
		}
	}	
	//cout << "temps Parois : " << user_time.time() - time << " seconds." << endl;
	user_time.reset();
	cout<<"volume solide := "<<volume_s<<endl;
}





/*!
* \fn double intersect_cube_tetrahedron(Bbox& cube, Tetrahedron& Tet) 
*\brief Intersection d'une box avec un t&eacute;tra&egrave;dre. 
\details Intersection de la bo&icirc;te rectangulaire (Bbox) \b cube avec le t&eacute;tra&egrave;dre \b Tet. Renvoie le volume de 
l'intersection. Fonction appelle lors du calcul de la quantit&eacute; balay&eacute;e. \n
Algorithme: \n
- Triangulation des faces du \a cube via la fonction  void triang_cellule(const Bbox&, Triangles& ). \n
- Recherche des sommets du t&eacute;tra&egrave;dre contenus dans \a cube via la fonction inside_box(const Bbox&, const Point_3& ). \n
- Boucle sur les faces du t&eacute;tra&egrave;dre.
- Test d'intersection entre \a cube et les faces du t&eacute;tra&egrave;dre via la fonction \b  CGAL::do_intersect(Bbox, Triangle_3). Si oui:
- Boucle sur les faces triangulaires du \a cube.
- Test d'intersection entre les faces du \a cube et les faces du t&eacute;tra&egrave;dre via la fonction  \b CGAL::do_intersect(Triangle_3, Triangle_3). Si oui: 
- Intersections des faces du \a t&eacute;tra&egrave;dre avec les faces du \a cube via la fonction  
\b CGAL::do_intersect(Triangle_3, Triangle_3).
- Calcul du volume du poly&egrave;dre r&eacute;sultant de l'intersection du \a cube et \a Tet. Pour le cacul du 
volume on construit avec les points d'intersections une triangulation(t&eacute;tra&egrave;dres) via la 
fonction \b  CGAL::Triangulation(vector<Point_3>) et le volume de l'intersection est la somme des volumes de ces t&eacute;tra&egrave;dres. 
Le volume d'un t&eacute;tra&egrave;dre est calcul&eacute; via la fonction \b CGAL::tetrahedron.volume().

*\param cube Box 3d 
*\param Tet T&eacute;tra&egrave;dre
*\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
*\return double
*/
double intersect_cube_tetrahedron(Bbox& cube, Tetrahedron& Tet){
	
	double volume=0.;
	
	if(abs(Tet.volume())>eps){
		
		Triangles triangCube;
		triang_cellule(cube , triangCube); 
		std::vector<Triangle_3> triangTet(4);
		std::vector<Point_3> Points_intersect;
		
		for(int i=0; i<4; i++){
			if (inside_box(cube, Tet.operator[](i)) ){
				Points_intersect.push_back(Tet.operator[](i));
			}
		}
		
		triangTet[0]= Triangle_3(Tet.operator[](0), Tet.operator[](1), Tet.operator[](2));
		triangTet[1]= Triangle_3(Tet.operator[](0), Tet.operator[](2), Tet.operator[](3));
		triangTet[2]= Triangle_3(Tet.operator[](0), Tet.operator[](1), Tet.operator[](3));
		triangTet[3]= Triangle_3(Tet.operator[](1), Tet.operator[](2), Tet.operator[](3));
		
		for(int i= 0; i<4; i++){
			if (CGAL::do_intersect(cube,triangTet[i]) ) {
				for(int j= 0; j<triangCube.size(); j++)
					if (CGAL::do_intersect(triangCube[j],triangTet[i]) ) {
						Triangle_3 t;
						Point_3 P;
						Segment_3 seg;
						std::vector<Point_3> vPoints; 
						
						CGAL::Object result = CGAL::intersection(triangCube[j],triangTet[i]);
						
						if(CGAL::assign(P,result)){
							Points_intersect.push_back(P);
						}
						else if(CGAL::assign(seg,result)){
							Points_intersect.push_back(seg.operator[](0));
							Points_intersect.push_back(seg.operator[](1));
						}
						else if(CGAL::assign(t,result)){
							Points_intersect.push_back(t.operator[](0));
							Points_intersect.push_back(t.operator[](1));
							Points_intersect.push_back(t.operator[](2));	
						}
						else if(CGAL::assign(vPoints,result)){ 
							for(int l= 0; l<vPoints.size(); l++)
							{
								Points_intersect.push_back(vPoints[l]);
							}
							
						}
						else {cout<<"Intersection type: ? dans <<intersect_cube_tetrahedron>> "<< endl;}
					}
			}
		}
		
		if( Points_intersect.size() >=4 ){
			Triangulation T(Points_intersect.begin(), Points_intersect.end());
			if(T.is_valid()){
				Finite_cells_iterator cit;
				for (cit = T.finite_cells_begin(); cit!= T.finite_cells_end(); cit++){
					volume+= CGAL::to_double(T.tetrahedron( cit).volume());
				}
			}
		}
	}
	return std::abs(volume);
}	




/*!
 * \fn void Grille::Parois(Solide& S,double dt) 
 *\brief Intersection de la grille fluide avec le solide.
 \details Intersection de la Grille fluide avec le Solide et calcul des diff&eacute;rents quantit&eacute;s d'int&eacute;r&ecirc;t: occupation du Solide dans la Cellule : \a Cellule.alpha, occupation des faces de la cellule par le solide: \a Cellule.kappai, \a Cellule.kappaj et \a Cellule.kappak. Remplissage des vecteurs : \n
- \a Particule.Points_interface: points d'intersections de la cellule avec les faces triangulaires du solide; \n
- \a Particule.Triangles_interface: d&eacute;coupage des faces des particules en morceaux triangulaires d'interface contenus dans une seule cellule de la grille fluide; \n
- \a Particule.Position_Triangles_interface: index de la Cellule fluide contenant \a Particule.Triangles_interface. \n
 
 Algorithme: \n
 - Construction du vecteur \a box_grille contenant les cellules cubiques de la grille fluide sous la forme de Box 3d (\a Bbox). Une Box repr&eacute;sente une bo&icirc;te rectangulaire. Cette fa&ccedil;on de voir les cellules fluide permet de faire appel aux fonctions membres de la classe <b> CGAL::Bbox_3  </b>. 
 - Construction du vecteur \a solide des Bbox associes aux \a Particule. 
 - Boucle sur \a box_grille.
 - Boucle sur \a solide.
 - Test d'intersection entre \a box_grille et \a solide  via la fonction \b CGAL::do_intersect(Bbox, Bbox). Sinon, il n'y a pas d'intersection et on passe directement &agrave; la cellule suivante. Si oui:
  - On test si \a box_grille est contenu dans \a solide via la fonction box_inside_convex_polygon(const Particule&, const Bbox&). Si oui,  l'intersection est \a box_grille, sinon:
   - Boucle sur les faces triangulaires du Solide.
   - Recherche des sommets de faces triangulaires du Solide contenues dans \a box_grille via la fonction inside_box(const Bbox&, const Point_3&). 
   - Test d'intersections entre \a box_grille est les faces triangulaires du Solide via la fonction \b CGAL::do_intersect(Bbox,Triangle_3). Si non, il n'y a rien a faire, si oui:
     - Triangulation des faces du \a box_grille via la fonction triang_cellule(const Bbox&, Triangles&) .
     - Boucle sur les faces triangulaires du \a box_grille.
     - Recherche des sommets de faces du \a box_grille contenues dans Solide via la fonction inside_convex_polygon(const Particule&, const Point_3&). 
     - Test d'intersection entre les faces du \a box_grille et les faces triangulaires du solide (\a Particule.triangle) via la fonction \b CGAL::do_intersect(Triangle_3, Triangle_3). Si non, il n'y a rien a faire, si oui:
     - Intersection entre les faces triangulaires du \a box_grille et les faces triangulaires du Solide via la fonction
     \b  CGAL::intersection(Triangle_3, Triangle_3).\n
      
      
      
      
\remark Le vecteur \a Particule.Points_interface, contenant les points d'intersections entre les cellules fluides et les faces triangulaires des particules (\a Particule.triangles) est rempli fur et &agrave; mesure dans l'algorithme de recherche des intersections. \n
On exploite les r&eacute;sultats des intersections afin de calculer les diff&eacute;rents quantit&eacute;s d'int&eacute;r&ecirc;t: \a Cellule.alpha, \a Cellule.kappai, \a Cellule.kappaj et \a Cellule.kappak. \n
Les fonctions utilis&eacute;es dans ce calcul sont: 
<b> CGAL::Triangulation(vector<Point_3>) </b> et <b> CGAL::convex_hull_3(vector<Point_3>, Polyhedron_3) </b>, <b> CGAL::tetrahedron.volume() </b> et <b> CGAL::Triangle_3.squared_area()</b>. \n
La derni&egrave;re &eacute;tape consiste dans le remplissage du vecteur \a Particule.Triangles_interface qui contient
le d&eacute;coupage des faces des particules en morceaux triangulaires d'interface contenus dans une seule
cellule de la grille fluide et, \a Particule.Position_Triangles_interface contenant l'index de la Cellule fluide contenant \a Particule.Triangles_interface. 
Le d&eacute;coupage des faces des particules en morceaux triangulaires d'interface se fait &agrave; partir du vecteur
\a Particule.Points_interface qui contient pour chaque face triangulaire des particules la liste de points d'intersections
de celles-ci avec la grille fluide. On fait simplement une triangulation contenant ces points via la fonction
<b> CGAL::Triangulation(vector<Point_3>)</b>.
 *\param S Solide
 *\param dt pas de temps
 *\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
 *\return void 
 */

void Grille::Parois_particles(Solide& S,double dt) {
	
	const double eps_relat = numeric_limits<double>::epsilon( );
	
	std::vector<Bbox> box_grille;
	const int nx_m=Nx+2*marge;
	const int ny_m=Ny+2*marge;
	const int nz_m=Nz+2*marge;
	const int Ns = (nx_m+1)*(ny_m+1)*(nz_m+1);
	const double volume_cel = deltax*deltay*deltaz;
	
	vector<vector<double> > Sommets(Ns, vector<double>(3,0.));
	int l=0;
	for(int i=0;i<nx_m+1;i++){
		for(int j=0;j<ny_m+1;j++){
			for(int k=0;k<nz_m+1;k++){
				Sommets[l][0] = (i-marge)*deltax;
				Sommets[l][1] = (j-marge)*deltay;
				Sommets[l][2] = (k-marge)*deltaz;
				l++;        
			}
		}
	}
	
	
	double x_min=0.,y_min=0.,z_min=0.,x_max=0.,y_max=0.,z_max=0.;
	
	for(int i=0; i<nx_m; i++){
		for(int j=0; j<ny_m; j++){ 
			for(int k=0; k<nz_m; k++){ 
				x_min = Sommets[k+j*(nz_m+1)+i*(nz_m+1)*(ny_m+1)][0];
				y_min = Sommets[k+j*(nz_m+1)+i*(nz_m+1)*(ny_m+1)][1];
				z_min = Sommets[k+j*(nz_m+1)+i*(nz_m+1)*(ny_m+1)][2];
				
				x_max = Sommets[(k+1)+(j+1)*(nz_m+1)+(i+1)*(nz_m+1)*(ny_m+1)][0];
				y_max = Sommets[(k+1)+(j+1)*(nz_m+1)+(i+1)*(nz_m+1)*(ny_m+1)][1];
				z_max = Sommets[(k+1)+(j+1)*(nz_m+1)+(i+1)*(nz_m+1)*(ny_m+1)][2];
				
				box_grille.push_back(Bbox(x_min,y_min,z_min,x_max,y_max,z_max));
				
			}
		}
	}
	
	
	
	int nb_particules = S.size();
	int taille = box_grille.size();
	cout<<"the grille size is : "<<taille<<endl;
	cout<<"nombres de particules : "<<nb_particules<<endl;
	
	
	const double eps_box = 0.1;
	std::vector<Bbox> solide(nb_particules);
	for(int it=0; it<nb_particules; it++){
		solide[it] = Bbox(S.solide[it].min_x-eps_box,S.solide[it].min_y-eps_box, S.solide[it].min_z-eps_box, S.solide[it].max_x+eps_box, S.solide[it].max_y+eps_box, S.solide[it].max_z+eps_box);
	}
	

	double volume_s=0.;
	
	Cellule cel;
	int i=0;
	CGAL::Timer user_time, user_time2;
	user_time.start();
	double time=0.;
	for (int a=0; a< nx_m; a++){
		for (int b=0; b< ny_m; b++){
			for (int c=0; c< nz_m; c++){
				
				std::vector< std::vector< std::vector<Point_3> > > Points_interface;
				Points_interface.resize(nb_particules);
				for(int count=0; count<nb_particules; count++){
				Points_interface[count].resize(S.solide[count].triangles.size(), std::vector<Point_3>(0));
				}
				cel = grille[a][b][c]; 
				cel.alpha = 0.; cel.kappai = 0.; cel.kappaj = 0.; cel.kappak = 0.;
				cel.phi_x = 0.; cel.phi_y = 0.; cel.phi_z = 0.; cel.phi_v = 0.; 
				cel.delta_w[0]= 0.; cel.delta_w[1]=0.; cel.delta_w[2]=0.; cel.delta_w[3]=0.; cel.delta_w[4] = 0.;
				Triangles trianglesB;
				bool exterieur = true;
				for(int iter_s=0; iter_s<nb_particules && exterieur; iter_s++){ //boucle sur les particules 
					std::vector<Point_3> Points_poly; 
					double alpha = 0.0;
					std::vector<double>  kappa(6,0.0);
					bool intersection = false;
					bool box_in_solide = false;
					bool point_in_solide = false;
					bool point_in_cell = false;
					if (CGAL::do_intersect(box_grille[i],solide[iter_s]) ) {	
						intersection = true;
						box_in_solide = box_inside_convex_polygon(S.solide[iter_s],box_grille[i]);
						if(box_in_solide){
							exterieur = false;
							cel.alpha = 1.;
							cel.kappai = 1.;
							cel.kappaj = 1.;
							cel.kappak = 1.;
							//test 3 janvier 2013
							if(a>0) {grille[a-1][b][c].kappai = 1.;}
							if(b>0) {grille[a][b-1][c].kappaj = 1.;}
							if(c>0) {grille[a][b][c-1].kappak = 1.;}
							//fin test 3 janvier 2013
							volume_s += volume_cel;
							box_in_solide = false;
						}
						else 
						{    
							triang_cellule(box_grille[i] , trianglesB); 
							
							for ( int j = 0; j < S.solide[iter_s].triangles.size(); j++){ 
								
								//test if point is in cell_box
								point_in_cell = inside_box(box_grille[i], S.solide[iter_s].triangles[j].operator[](0));
								if(point_in_cell) {Points_poly.push_back(S.solide[iter_s].triangles[j].operator[](0)); point_in_cell=false;
								Points_interface[iter_s][j].push_back(S.solide[iter_s].triangles[j].operator[](0));
								}
								
								point_in_cell = inside_box(box_grille[i], S.solide[iter_s].triangles[j].operator[](1));
								if(point_in_cell) {Points_poly.push_back(S.solide[iter_s].triangles[j].operator[](1)); point_in_cell=false;
								Points_interface[iter_s][j].push_back(S.solide[iter_s].triangles[j].operator[](1));
								}
								
								point_in_cell = inside_box(box_grille[i], S.solide[iter_s].triangles[j].operator[](2));
								if(point_in_cell) {Points_poly.push_back(S.solide[iter_s].triangles[j].operator[](2)); point_in_cell=false;
								Points_interface[iter_s][j].push_back(S.solide[iter_s].triangles[j].operator[](2));
								}
								
								
								if (CGAL::do_intersect(box_grille[i],S.solide[iter_s].triangles[j]) ) {
									for ( int k = 0; k < trianglesB.size(); k++){
										//test if point is in solide
										point_in_solide = inside_convex_polygon(S.solide[iter_s],trianglesB[k].operator[](0));
										if(point_in_solide) {Points_poly.push_back(trianglesB[k].operator[](0)); point_in_solide=false;
										}
										
										point_in_solide = inside_convex_polygon(S.solide[iter_s],trianglesB[k].operator[](1));
										if(point_in_solide) {Points_poly.push_back(trianglesB[k].operator[](1)); point_in_solide=false;
										}
										
										point_in_solide = inside_convex_polygon(S.solide[iter_s],trianglesB[k].operator[](2));
										if(point_in_solide) {Points_poly.push_back(trianglesB[k].operator[](2)); point_in_solide=false;
										}
										
										if (CGAL::do_intersect(S.solide[iter_s].triangles[j],trianglesB[k]) ) {
											
											Triangle_3 t;
											Point_3 P;
											Segment_3 seg;
											std::vector<Point_3> vPoints; 
											CGAL::Object result = CGAL::intersection(S.solide[iter_s].triangles[j], trianglesB[k]);
											
											if(CGAL::assign(P,result)){
												Points_poly.push_back(P);
												Points_interface[iter_s][j].push_back(P);
											}
											else if(CGAL::assign(seg,result)){
												Points_poly.push_back(seg.operator[](0));
												Points_poly.push_back(seg.operator[](1));
												Points_interface[iter_s][j].push_back(seg.operator[](0));
												Points_interface[iter_s][j].push_back(seg.operator[](1));
											}
											else if(CGAL::assign(t,result)){
												Points_poly.push_back(t.operator[](0));
												Points_poly.push_back(t.operator[](1));
												Points_poly.push_back(t.operator[](2));
												Points_interface[iter_s][j].push_back(t.operator[](0));
												Points_interface[iter_s][j].push_back(t.operator[](1));
												Points_interface[iter_s][j].push_back(t.operator[](2));
												
											}
											else if(CGAL::assign(vPoints,result)){ 
												for(int l= 0; l<vPoints.size(); l++)
												{
													Points_poly.push_back(vPoints[l]);
													Points_interface[iter_s][j].push_back(vPoints[l]);
												}
												
											}
											else {cout<<"Intersection type: ? in intersect particule with grille"<< trianglesB[k]<<endl;}
										}
										//end intr t1 et t2
									} //end boucle sur trianglesB
									
								} // if intersection cellule[i] avec triangle[j]	 
								
							} // end boucle sur triangles
						} //else 
						
					} // if inter grille[i] et solide[iter_s] 
					
					//traitement calcul de alpha et kappa pour la cellule=grille[i]!!!!!!
					user_time2.start();
					if(intersection && exterieur){
						Triangulation T(Points_poly.begin(), Points_poly.end());
						
						if (T.dimension() == 3){
							Polyhedron_3 poly;
							CGAL::convex_hull_3(T.points_begin(), T.points_end(), poly);
							Finite_cells_iterator cit;
							for (cit = T.finite_cells_begin(); cit!= T.finite_cells_end(); cit++){
								alpha+= CGAL::to_double(T.tetrahedron( cit).volume());
							}
							
							Facet_iterator fiter;
							for (fiter = poly.facets_begin(); fiter!= poly.facets_end(); fiter++){
								
								Triangle_3 K((*fiter).halfedge()->vertex()->point(),(*fiter).halfedge()->next()->vertex()->point(),
														 (*fiter).halfedge()->opposite()->vertex()->point());
														 
														 if (abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](0).operator[](2))<=eps_relat && abs(trianglesB[0].operator[](0).operator[](2) - K.operator[](1).operator[](2))<=eps_relat && abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](2).operator[](2))<=eps_relat )
														 { 
															 kappa[0] +=sqrt(CGAL::to_double(K.squared_area()));
														 }
														 
														 else if (abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](0).operator[](2))<=eps_relat && abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](1).operator[](2))<=eps_relat && abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](2).operator[](2))<=eps_relat )
														 { 
															 kappa[1] +=sqrt(CGAL::to_double(K.squared_area()));
														 }
														 
														 else if (abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](0).operator[](0))<=eps_relat && abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](1).operator[](0))<=eps_relat && abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](2).operator[](0))<=eps_relat)
														 { 
															 kappa[2] +=sqrt(CGAL::to_double(K.squared_area()));
														 }
														 
														 else if (abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](0).operator[](0))<=eps_relat && abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](1).operator[](0))<=eps_relat && abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](2).operator[](0))<=eps_relat)
														 { 
															 kappa[3] +=sqrt(CGAL::to_double(K.squared_area()));
														 }
														 
														 else if (abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](0).operator[](1))<=eps_relat && abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](1).operator[](1))<=eps_relat && abs(trianglesB[8].operator[](0).operator[](1) - K.operator[](2).operator[](1))<=eps_relat)
														 { 
															 kappa[4] +=sqrt(CGAL::to_double(K.squared_area()));
														 }
														 
														 else if (abs(trianglesB[10].operator[](0).operator[](1) -  K.operator[](0).operator[](1))<=eps_relat && abs(trianglesB[10].operator[](0).operator[](1) -  K.operator[](1).operator[](1))<=eps_relat && abs(trianglesB[10].operator[](0).operator[](1) - K.operator[](2).operator[](1))<=eps_relat)
														 { 
															 kappa[5] +=sqrt(CGAL::to_double(K.squared_area()));
														 }
														 
														 else{
															 
														 } //calcul des aires parietales
														 
							}
						}
						
						if (T.dimension() == 2){
							Finite_faces_iterator it;
							for (it = T.finite_facets_begin(); it != T.finite_facets_end(); it++){
								Triangle_3 K= T.triangle(*it);
								
								if (abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](0).operator[](2))<=eps_relat && abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](1).operator[](2))<=eps_relat && abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](2).operator[](2))<=eps_relat)
								{ 
									kappa[0] +=sqrt(CGAL::to_double(K.squared_area()));
									
								}
								
								else if (abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](0).operator[](2))<=eps_relat && abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](1).operator[](2))<=eps_relat && abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](2).operator[](2))<=eps_relat)
								{ 
									kappa[1] +=sqrt(CGAL::to_double(K.squared_area()));
									
								}
								
								else if (abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](0).operator[](0))<=eps_relat && abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](1).operator[](0))<=eps_relat && abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](2).operator[](0))<=eps_relat)
								{ 
									kappa[2] +=sqrt(CGAL::to_double(K.squared_area()));
									
								}
								
								else if (abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](0).operator[](0))<=eps_relat && abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](1).operator[](0))<=eps_relat && abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](2).operator[](0))<=eps_relat)
								{ 
									kappa[3] +=sqrt(CGAL::to_double(K.squared_area()));
									
								}
								
								else if (abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](0).operator[](1))<=eps_relat && abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](1).operator[](1))<=eps_relat && abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](2).operator[](1))<=eps_relat)
								{ 
									kappa[4] +=sqrt(CGAL::to_double(K.squared_area()));
									
								}
								
								else 
								{ 
									kappa[5] +=sqrt(CGAL::to_double(K.squared_area()));
									
								}
							}
							
						}
						
						cel.alpha  += alpha/volume_cel;
						cel.kappai += kappa[3]/(deltay * deltaz);
						cel.kappaj += kappa[4]/(deltax * deltaz);
						cel.kappak += kappa[1]/(deltax * deltay);
						//test 30 juillet 2013
						if(cel.kappai >=1.) {cel.kappai=1.;}
						if(cel.kappaj >=1.) {cel.kappaj=1.;}
						if(cel.kappak >=1.) {cel.kappak=1.;}
						if(cel.alpha >=1.) {cel.alpha=1.;}
						//fin test 30 juillet 2013
						
						volume_s +=alpha;
						time+= user_time2.time();
						user_time2.reset();
						
					}// fin calcul alpha et kappa
					
				} //fin boucle sur les particules
				
         grille[a][b][c] = cel;
				//test 31 janvier 2013
				if(std::abs(grille[a][b][c].alpha -1.) <1.e-10) {
					grille[a][b][c].alpha = 1.;
					grille[a][b][c].kappai = 1.;
					grille[a][b][c].kappaj = 1.;
					grille[a][b][c].kappak = 1.;
					if(a>0){grille[a-1][b][c].kappai = 1.;}
					if(b>0){grille[a][b-1][c].kappaj = 1.;}
					if(c>0){grille[a][b][c-1].kappak = 1.;}
				}
				//fin test 31 janvier 2013
				i++;
				
				//triangularisation de l'interface face par face
	      Finite_faces_iterator iter;
				for(int count=0; count<nb_particules;count++){//boucle sur les Particules
					for(int it=0; it<S.solide[count].triangles.size(); it++){ //boucle sur les Particules.triangles
						Triangulation T(Points_interface[count][it].begin(), Points_interface[count][it].end());
							assert(T.is_valid());
							if(T.dimension()==2){
								for (iter = T.finite_facets_begin(); iter != T.finite_facets_end(); iter++){
									  Triangle_3 Tri= T.triangle(*iter);
										if(std::sqrt(CGAL::to_double(Tri.squared_area())) >eps){
										Vector_3 vect0(Tri.operator[](0),Tri.operator[](1));
										Vector_3 vect1(Tri.operator[](0),Tri.operator[](2));
										Vector_3 normale = CGAL::cross_product(vect0,vect1);
										if (normale*S.solide[count].normales[it] > 0.){
											S.solide[count].Triangles_interface[it].push_back(Tri);
										}
										else {
											S.solide[count].Triangles_interface[it].push_back(Triangle_3(Tri.operator[](0),Tri.operator[](2),Tri.operator[](1)));
										}
										std::vector<int> poz(3); poz[0]= a; poz[1] = b; poz[2] = c;
										S.solide[count].Position_Triangles_interface[it].push_back(poz);
									}
								}
							}
					} //fin boucle sur Triangles
				} //fin boucle sur les particules
				
			} //fin boucle sur grille
		}
	}	
	//cout << "temps Parois : " << user_time.time() - time << " seconds." << endl;
	user_time.reset();
	cout<<"volume solide := "<<volume_s<<endl;
}

