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
							//box_in_solide = false;
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
											
										  //Triangle_3 t;
										  //Point_3 P;
										  //Segment_3 seg;
										  //std::vector<Point_3> vPoints; 
										  //CGAL::Object result = CGAL::intersection(S.solide[iter_s].triangles[j], trianglesB[k]);
										  std::vector<Point_3> result = intersection_bis(S.solide[iter_s].triangles[j], trianglesB[k]);
										  
										  
										  /*if(CGAL::assign(P,result)){
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
										  else {cout<<"Intersection type: ? in intersect particule with grille"<< trianglesB[k]<<endl;}*/
										  for(int l= 0; l<result.size(); l++)
										  {
										    Points_poly.push_back(result[l]);
										    Points_interface[iter_s][j].push_back(result[l]);
										  }
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
				if(intersection && exterieur && Points_poly.size()>0){
				  //Triangulation T(Points_poly.begin(), Points_poly.end());
				  std::vector<Point_3> Points_poly2 = redondances(Points_poly.begin(), Points_poly.end());
					
				  //if (T.dimension() == 3){
						Polyhedron_3 poly;
						//CGAL::convex_hull_3(T.points_begin(), T.points_end(), poly);
						CGAL::convex_hull_3(Points_poly2.begin(), Points_poly2.end(), poly);
						//int l=0, n=0;
						/*Finite_cells_iterator cit;
						for (cit = poly.finite_cells_begin(); cit!= poly.finite_cells_end(); cit++){
							alpha+= CGAL::to_double(T.tetrahedron( cit).volume());
							}*/
						Point_3 P = (*(poly.facets_begin())).halfedge()->vertex()->point();
						Facet_iterator fiter;
						for (fiter = poly.facets_begin(); fiter!= poly.facets_end(); fiter++){
						  Tetrahedron T(P,(*fiter).halfedge()->vertex()->point(),(*fiter).halfedge()->next()->vertex()->point(), (*fiter).halfedge()->opposite()->vertex()->point());
						  alpha+= CGAL::to_double(T.volume());
						}
						
						//Facet_iterator fiter;
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
						/*}
					
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
							}*/
						
						//}
					
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
	//cout<<"volume solide := "<<volume_s<<endl;
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
		
/*
		Point_3 p1(cube.xmin(),cube.ymin(),cube.zmin());
		if(inside_tetra(Tet,p1)) {Points_intersect.push_back(p1);}
		
		Point_3 p2(cube.xmax(),cube.ymax(),cube.zmax());
		if(inside_tetra(Tet,p2)) {Points_intersect.push_back(p2);}
		
		Point_3 p3(cube.xmax(),cube.ymin(),cube.zmin());
		if(inside_tetra(Tet,p3)) {Points_intersect.push_back(p3);}
		
		Point_3 p4(cube.xmax(),cube.ymin(),cube.zmax());
		if(inside_tetra(Tet,p4)) {Points_intersect.push_back(p4);}
		
		Point_3 p5(cube.xmax(),cube.ymax(),cube.zmin());
		if(inside_tetra(Tet,p5)) {Points_intersect.push_back(p5);}
		
		Point_3 p6(cube.xmin(),cube.ymax(),cube.zmin());
		if(inside_tetra(Tet,p6)) {Points_intersect.push_back(p6);}
		
		Point_3 p7(cube.xmin(),cube.ymax(),cube.zmax());
		if(inside_tetra(Tet,p7)) {Points_intersect.push_back(p7);}
		
		Point_3 p8(cube.xmin(),cube.ymin(),cube.zmax());
		if(inside_tetra(Tet,p8)) {Points_intersect.push_back(p8);}*/

		
		triangTet[0]= Triangle_3(Tet.operator[](0), Tet.operator[](1), Tet.operator[](2));
		triangTet[1]= Triangle_3(Tet.operator[](0), Tet.operator[](2), Tet.operator[](3));
		triangTet[2]= Triangle_3(Tet.operator[](0), Tet.operator[](1), Tet.operator[](3));
		triangTet[3]= Triangle_3(Tet.operator[](1), Tet.operator[](2), Tet.operator[](3));
		
		for(int i= 0; i<4; i++){
			if (CGAL::do_intersect(cube,triangTet[i]) ) {
				for(int j= 0; j<triangCube.size(); j++)
					if (CGAL::do_intersect(triangCube[j],triangTet[i]) ) {
					  //Triangle_3 t;
					  //Point_3 P;
					  //Segment_3 seg;
					  //std::vector<Point_3> vPoints; 
					  
					  //CGAL::Object result = CGAL::intersection(triangCube[j],triangTet[i]);
					  std::vector<Point_3> result = intersection_bis(triangCube[j],triangTet[i]);
					  
					  /*if(CGAL::assign(P,result)){
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
					  else {cout<<"Intersection type: ? dans <<intersect_cube_tetrahedron>> "<< endl;}*/
					  for(int l= 0; l<result.size(); l++)
					  {
					    Points_intersect.push_back(result[l]);
					  }
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
* \fn void Grille::Parois_particles(Solide& S,double dt) 
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
   - Test d'intersections entre \a box_grille et les faces triangulaires du Solide via la fonction \b CGAL::do_intersect(Bbox,Triangle_3). Si non, il n'y a rien a faire, si oui:
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
  CGAL::Timer total_time,bbox_time,do_intersect_time,triangularisation_time,test_time,test_inside_time,alpha_time,intersect_time,convex_hull_time,volume_time,kappa_time1,kappa_time2,triangulation_time;
  double temps_total=0.,temps_bbox=0.,temps_do_intersect=0.,temps_triangularisation=0.,temps_test=0.,temps_test_inside=0.,temps_alpha=0.,temps_intersect=0.,nb_intersect=0.,temps_convex_hull=0.,nb_convex_hull=0.,temps_volume=0.,temps_kappa1=0.,temps_kappa2=0.,nb_kappa1=0.,nb_kappa2=0.,temps_triangulation=0.;
  
  total_time.start();
  bbox_time.start();
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
	int nb_triangles=0.;
	for(int iter=0; iter<nb_particules; iter++){
		nb_triangles += S.solide[iter].triangles.size();
	}
	cout<<"nombres de traiangles : "<<nb_triangles<<endl;
	
	const double eps_box = 0.1;
	std::vector<Bbox> solide(nb_particules);
	for(int it=0; it<nb_particules; it++){
		solide[it] = Bbox(S.solide[it].min_x-eps_box,S.solide[it].min_y-eps_box, S.solide[it].min_z-eps_box, S.solide[it].max_x+eps_box, S.solide[it].max_y+eps_box, S.solide[it].max_z+eps_box);
	}
	
	temps_bbox += CGAL::to_double(bbox_time.time());
	

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
				do_intersect_time.start();
				for(int iter_s=0; iter_s<nb_particules && exterieur; iter_s++){ //boucle sur les particules 
				  test_time.start();
				  bool test = CGAL::do_intersect(box_grille[i],solide[iter_s]);
				  temps_test += CGAL::to_double(test_time.time());
				  test_time.reset();
				  if (test){//CGAL::do_intersect(box_grille[i],solide[iter_s]) ) {
				    std::vector<Point_3> Points_poly; 
				    double alpha = 0.0;
				    std::vector<double>  kappa(6,0.0);
				    bool intersection = false;
				    bool box_in_solide = false;
				    bool point_in_solide = false;
				    bool point_in_cell = false;
				    intersection = true;
				    test_inside_time.start();
				    box_in_solide = box_inside_convex_polygon(S.solide[iter_s],box_grille[i]);
				    temps_test_inside += CGAL::to_double(test_inside_time.time());
				    test_inside_time.reset();
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
					  //Test pour savoir si les coins de la cellule sont dans le solide
					  point_in_solide = inside_convex_polygon(S.solide[iter_s],Point_3(box_grille[i].xmin(),box_grille[i].ymin(),box_grille[i].zmin()));
					  if(point_in_solide) {Points_poly.push_back(Point_3(box_grille[i].xmin(),box_grille[i].ymin(),box_grille[i].zmin())); point_in_solide=false;
					  }

					  point_in_solide = inside_convex_polygon(S.solide[iter_s],Point_3(box_grille[i].xmin(),box_grille[i].ymin(),box_grille[i].zmax()));
					  if(point_in_solide) {Points_poly.push_back(Point_3(box_grille[i].xmin(),box_grille[i].ymin(),box_grille[i].zmax())); point_in_solide=false;
					  }

					  point_in_solide = inside_convex_polygon(S.solide[iter_s],Point_3(box_grille[i].xmin(),box_grille[i].ymax(),box_grille[i].zmin()));
					  if(point_in_solide) {Points_poly.push_back(Point_3(box_grille[i].xmin(),box_grille[i].ymax(),box_grille[i].zmin())); point_in_solide=false;
					  }

					  point_in_solide = inside_convex_polygon(S.solide[iter_s],Point_3(box_grille[i].xmin(),box_grille[i].ymax(),box_grille[i].zmax()));
					  if(point_in_solide) {Points_poly.push_back(Point_3(box_grille[i].xmin(),box_grille[i].ymax(),box_grille[i].zmax())); point_in_solide=false;
					  }

					  point_in_solide = inside_convex_polygon(S.solide[iter_s],Point_3(box_grille[i].xmax(),box_grille[i].ymin(),box_grille[i].zmin()));
					  if(point_in_solide) {Points_poly.push_back(Point_3(box_grille[i].xmax(),box_grille[i].ymin(),box_grille[i].zmin())); point_in_solide=false;
					  }

					  point_in_solide = inside_convex_polygon(S.solide[iter_s],Point_3(box_grille[i].xmax(),box_grille[i].ymin(),box_grille[i].zmax()));
					  if(point_in_solide) {Points_poly.push_back(Point_3(box_grille[i].xmax(),box_grille[i].ymin(),box_grille[i].zmax())); point_in_solide=false;
					  }

					  point_in_solide = inside_convex_polygon(S.solide[iter_s],Point_3(box_grille[i].xmax(),box_grille[i].ymax(),box_grille[i].zmin()));
					  if(point_in_solide) {Points_poly.push_back(Point_3(box_grille[i].xmax(),box_grille[i].ymax(),box_grille[i].zmin())); point_in_solide=false;
					  }

					  point_in_solide = inside_convex_polygon(S.solide[iter_s],Point_3(box_grille[i].xmax(),box_grille[i].ymax(),box_grille[i].zmax()));
					  if(point_in_solide) {Points_poly.push_back(Point_3(box_grille[i].xmax(),box_grille[i].ymax(),box_grille[i].zmax())); point_in_solide=false;
					  }

					  //Intersections des aretes de la cellule avec le triangle
					  for(int kx=0;kx<2;kx++){
					    for(int ky=0;ky<2;ky++){
					      for(int kz=0;kz<2;kz++){
						double x1 = box_grille[i].xmin()+kx*(box_grille[i].xmax()-box_grille[i].xmin());
						double y1 = box_grille[i].ymin()+ky*(box_grille[i].ymax()-box_grille[i].ymin());
						double z1 = box_grille[i].zmin()+kz*(box_grille[i].zmax()-box_grille[i].zmin());
						if(kx==0){
						  double x2 = box_grille[i].xmax();
						  Segment_3 seg(Point_3(x1,y1,z1),Point_3(x2,y1,z1));
						  if (CGAL::do_intersect(seg,S.solide[iter_s].triangles[j]) ) {
						    intersect_time.start();
						    nb_intersect+=1.;
						    std::vector<Point_3> result = intersection_bis(seg,S.solide[iter_s].triangles[j]);
						    temps_intersect += CGAL::to_double(intersect_time.time());
						    intersect_time.reset();
						    
						    for(int l= 0; l<result.size(); l++)
						    {
						      Points_poly.push_back(result[l]);
						      Points_interface[iter_s][j].push_back(result[l]);
						    }
						  }
						  //end intersection seg et t1
						}
						if(ky==0){
						  double y2 = box_grille[i].ymax();
						  Segment_3 seg(Point_3(x1,y1,z1),Point_3(x1,y2,z1));
						  if (CGAL::do_intersect(seg,S.solide[iter_s].triangles[j]) ) {
						    intersect_time.start();
						    nb_intersect+=1.;
						    std::vector<Point_3> result = intersection_bis(seg,S.solide[iter_s].triangles[j]);
						    temps_intersect += CGAL::to_double(intersect_time.time());
						    intersect_time.reset();
						    
						    for(int l= 0; l<result.size(); l++)
						    {
						      Points_poly.push_back(result[l]);
						      Points_interface[iter_s][j].push_back(result[l]);
						    }
						  }
						  //end intersection seg et t1
						}
						if(kz==0){
						  double z2 = box_grille[i].xmax();
						  Segment_3 seg(Point_3(x1,y1,z1),Point_3(x1,y1,z2));
						  if (CGAL::do_intersect(seg,S.solide[iter_s].triangles[j]) ) {
						    intersect_time.start();
						    nb_intersect+=1.;
						    std::vector<Point_3> result = intersection_bis(seg,S.solide[iter_s].triangles[j]);
						    temps_intersect += CGAL::to_double(intersect_time.time());
						    intersect_time.reset();
						    
						    for(int l= 0; l<result.size(); l++)
						    {
						      Points_poly.push_back(result[l]);
						      Points_interface[iter_s][j].push_back(result[l]);
						    }
						  }
						  //end intersection seg et t1
						}
					      }
					    }
					  }
					  

					  //Intersections des aretes de la face solide avec les triangles de la cellule
					  for ( int k = 0; k < trianglesB.size(); k++){	
					    if (CGAL::do_intersect(S.solide[iter_s].triangles[j],trianglesB[k]) ) {
					      for(int l=0;l<3;l++){
						int lp = (l+1)%3;
						Segment_3 seg(S.solide[iter_s].triangles[j].operator[](l),S.solide[iter_s].triangles[j].operator[](lp));
						
						intersect_time.start();
						nb_intersect+=1.;
						std::vector<Point_3> result = intersection_bis(seg,trianglesB[k]);
						temps_intersect += CGAL::to_double(intersect_time.time());
						intersect_time.reset();
						
						for(int l= 0; l<result.size(); l++)
						{
						  Points_poly.push_back(result[l]);
						  Points_interface[iter_s][j].push_back(result[l]);
						}
					      }
					      //end intersection seg et t2
					    }
					  } //end boucle sur trianglesB
									
					} // if intersection cellule[i] avec triangle[j]	 
								
				      } // end boucle sur triangles
				    } //else 
						
					
					
				    //traitement calcul de alpha et kappa pour la cellule=grille[i]!!!!!!
				    user_time2.start();
				    alpha_time.start();
				    if(intersection && exterieur && Points_poly.size()>0){
				      
				      triangulation_time.start();
				      //Triangulation Tr(Points_poly.begin(), Points_poly.end());
				      std::vector<Point_3> Points_poly2 = redondances(Points_poly.begin(), Points_poly.end());
				      /*cout << "Points_poly.size()=" << Points_poly.size() << " Points_poly2.size()=" << Points_poly2.size() << endl;
				      //getchar();
				      cout << "Points_poly:" << endl;
				      for(std::vector<Point_3>::iterator it=Points_poly.begin();it!=Points_poly.end();it++){
					cout << (*it).x() << " " << (*it).y() << " " << (*it).z() << endl;
					getchar();
				      }*/
				      /*cout << "Triangulation:" << endl;
				      for(Point_iterator it=Tr.points_begin();it!=Tr.points_end();it++){
					cout << (*it).x() << " " << (*it).y() << " " << (*it).z() << endl;
					getchar();
					}*/
				      
				      
				      temps_triangulation += CGAL::to_double(triangulation_time.time());
				      triangulation_time.reset();
				      
				      //if (T.dimension() == 3){
				      Polyhedron_3 poly;
				      convex_hull_time.start();
				      //cout << "before convex_hull size=" << Points_poly.size() << endl;
				      CGAL::convex_hull_3(Points_poly2.begin(), Points_poly2.end(), poly);
				      //CGAL::convex_hull_3(Tr.points_begin(), Tr.points_end(), poly);
//getchar();
				      //cout << "after convex_hull facets=" << poly.size_of_facets() << " vertices=" << poly.size_of_vertices() << endl;
					nb_convex_hull += 1.;
					temps_convex_hull += CGAL::to_double(convex_hull_time.time());
					convex_hull_time.reset();
					volume_time.start();
					//Finite_cells_iterator cit;
					//for (cit = poly.finite_cells_begin(); cit!= poly.finite_cells_end(); cit++){
					//alpha+= CGAL::to_double(poly.tetrahedron( cit).volume());
					//}
					Point_3 P = (*(poly.facets_begin())).halfedge()->vertex()->point();
					Facet_iterator fiter;
					for (fiter = poly.facets_begin(); fiter!= poly.facets_end(); fiter++){
					  Tetrahedron T(P,(*fiter).halfedge()->vertex()->point(),(*fiter).halfedge()->next()->vertex()->point(), (*fiter).halfedge()->opposite()->vertex()->point());
					  alpha+= CGAL::to_double(T.volume());
					}
					
					
					temps_volume += CGAL::to_double(volume_time.time());
					volume_time.reset();
					
					kappa_time1.start();
					//Facet_iterator fiter;
					for (fiter = poly.facets_begin(); fiter!= poly.facets_end(); fiter++){
					  nb_kappa1 +=1.;
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
					temps_kappa1 += CGAL::to_double(kappa_time1.time());
					kappa_time1.reset();
					//}
						
				    /*if (T.dimension() == 2){
					kappa_time2.start();
					Finite_faces_iterator it;
					for (it = T.finite_facets_begin(); it != T.finite_facets_end(); it++){
					  nb_kappa2 +=1.;
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
					temps_kappa2 += CGAL::to_double(kappa_time2.time());
					kappa_time2.reset();
					
					}*/
						
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
				    temps_alpha += CGAL::to_double(alpha_time.time());
				    alpha_time.reset();
				    
				  }//fin du test do_intersect
					
					
				} //fin boucle sur les particules
				temps_do_intersect += CGAL::to_double(do_intersect_time.time());
				do_intersect_time.reset();
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
				triangularisation_time.reset();
				triangularisation_time.start();
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
				temps_triangularisation += CGAL::to_double(triangularisation_time.time());
				triangularisation_time.reset();
			} //fin boucle sur grille
		}
	}	
	//cout << "temps Parois : " << user_time.time() - time << " seconds." << endl;
	user_time.reset();
	cout<<"volume solide parois := "<<volume_s<<endl;
	temps_total = CGAL::to_double(total_time.time());
	
	cout << "######### COUTS INTERSECTIONS ##########" << endl;
	cout << "Bbox=" << 100*temps_bbox/temps_total << "%" << endl;
	cout << "do_intersect=" << 100*temps_do_intersect/temps_total << "%" << endl;
	cout << "   test intersect=" << 100*temps_test/temps_total << "%" << endl;
	cout << "   test_inside=" << 100*temps_test_inside/temps_total << "%" << endl;
	cout << "   intersect=" << 100*temps_intersect/temps_total << "%          t_moy=" << temps_intersect/nb_intersect << " nb_intersect=" << nb_intersect << endl;
	cout << "   alpha=" << 100*temps_alpha/temps_total << "%" << endl;
	cout << "      triangulation=" << 100*temps_triangulation/temps_total << "%          t_moy=" << temps_triangulation/nb_convex_hull << " nb_triangulation=" << nb_convex_hull << endl;
	cout << "      convex_hull=" << 100*temps_convex_hull/temps_total << "%          t_moy=" << temps_convex_hull/nb_convex_hull << " nb_convex_hull=" << nb_convex_hull << endl;
	cout << "      volume=" << 100*temps_volume/temps_total << "%" << endl;
	cout << "      kappa 3d=" << 100*temps_kappa1/temps_total << "%          t_moy=" << temps_kappa1/nb_kappa1 << " nb_kappa1=" << nb_kappa1 << endl;
	cout << "      kappa 2d=" << 100*temps_kappa2/temps_total << "%          t_moy=" << temps_kappa2/nb_kappa2 << " nb_kappa2=" << nb_kappa2 << endl;
	cout << "triangularisation=" << 100*temps_triangularisation/temps_total << "%" << endl;
	cout << "Reste=" << 100-100*(temps_bbox+temps_do_intersect+temps_triangularisation)/temps_total << "%" << endl;
	cout << "########################################" << endl;
	
}

//test septembre 2013
//cout<<"volume intersection solide decomposition en tetra avec la grille fluide "<<endl;
//Fluide.Parois_tetra(S,dt);
// void Grille::Parois_tetra(Solide& S,double dt) {
// 	
// 	const double eps_relat = numeric_limits<double>::epsilon( );
// 	
// 	std::vector<Bbox> box_grille;
// 	const int nx_m=Nx+2*marge;
// 	const int ny_m=Ny+2*marge;
// 	const int nz_m=Nz+2*marge;
// 	const int Ns = (nx_m+1)*(ny_m+1)*(nz_m+1);
// 	const double volume_cel = deltax*deltay*deltaz;
// 	
// 	vector<vector<double> > Sommets(Ns, vector<double>(3,0.));
// 	int l=0;
// 	for(int i=0;i<nx_m+1;i++){
// 		for(int j=0;j<ny_m+1;j++){
// 			for(int k=0;k<nz_m+1;k++){
// 				Sommets[l][0] = (i-marge)*deltax;
// 				Sommets[l][1] = (j-marge)*deltay;
// 				Sommets[l][2] = (k-marge)*deltaz;
// 				l++;        
// 			}
// 		}
// 	}
// 	
// 	
// 	double x_min=0.,y_min=0.,z_min=0.,x_max=0.,y_max=0.,z_max=0.;
// 	
// 	for(int i=0; i<nx_m; i++){
// 		for(int j=0; j<ny_m; j++){ 
// 			for(int k=0; k<nz_m; k++){ 
// 				x_min = Sommets[k+j*(nz_m+1)+i*(nz_m+1)*(ny_m+1)][0];
// 				y_min = Sommets[k+j*(nz_m+1)+i*(nz_m+1)*(ny_m+1)][1];
// 				z_min = Sommets[k+j*(nz_m+1)+i*(nz_m+1)*(ny_m+1)][2];
// 				
// 				x_max = Sommets[(k+1)+(j+1)*(nz_m+1)+(i+1)*(nz_m+1)*(ny_m+1)][0];
// 				y_max = Sommets[(k+1)+(j+1)*(nz_m+1)+(i+1)*(nz_m+1)*(ny_m+1)][1];
// 				z_max = Sommets[(k+1)+(j+1)*(nz_m+1)+(i+1)*(nz_m+1)*(ny_m+1)][2];
// 				
// 				box_grille.push_back(Bbox(x_min,y_min,z_min,x_max,y_max,z_max));
// 				
// 			}
// 		}
// 	}
// 	
// 	
// 	
// 	int nb_particules = S.size();
// 	int taille = box_grille.size();
// 	
// 	
// 	const double eps_box = 0.1;
// 	std::vector<Bbox> solide(nb_particules);
// 	for(int it=0; it<nb_particules; it++){
// 		solide[it] = Bbox(S.solide[it].min_x-eps_box,S.solide[it].min_y-eps_box, S.solide[it].min_z-eps_box, S.solide[it].max_x+eps_box, S.solide[it].max_y+eps_box, S.solide[it].max_z+eps_box);
// 	}
// 	
// 	
// 	double volume_s=0.;
// 	Cellule cel;
// 	int i=0;
// 
// 	for (int a=0; a< nx_m; a++){
// 		for (int b=0; b< ny_m; b++){
// 			for (int c=0; c< nz_m; c++){
// 
// 				cel = grille[a][b][c]; 
// 				for(int iter_s=0; iter_s<nb_particules; iter_s++){
// 					
// 					if (CGAL::do_intersect(box_grille[i],solide[iter_s]) ) {
// 					
// 								std::vector<Point_3> Points_particule; 
// 								Point_3 center;
// 								for(int l= 0; l< S.solide[iter_s].triangles.size(); l++)
// 								{
// 									Points_particule.push_back(S.solide[iter_s].triangles[l].operator[](0));
// 									Points_particule.push_back(S.solide[iter_s].triangles[l].operator[](1));
// 									Points_particule.push_back(S.solide[iter_s].triangles[l].operator[](2));
// 								}	
// 								center = centroid(Points_particule.begin(),Points_particule.end());
// 								
// 								for(int l= 0; l<S.solide[iter_s].triangles.size(); l++)
// 								{
// 									Tetrahedron tetra(S.solide[iter_s].triangles[l].operator[](0), S.solide[iter_s].triangles[l].operator[](1), S.solide[iter_s].triangles[l].operator[](2),center);
// 									
// 									
// 									Bbox box_tetra= tetra.bbox();
// 									if (CGAL::do_intersect(box_grille[i],box_tetra)) {
// 										
// 										if(box_inside_tetra(tetra,box_grille[i])){
// 											volume_s += std::abs(CGAL::to_double(tetra.volume()));
// 										}
// 										else 
// 										{
// 										if (inside_box(box_grille[i],tetra.operator[](0)) &&  inside_box(box_grille[i], tetra.operator[](1))
// 											&& inside_box(box_grille[i],tetra.operator[](2)) && inside_box(box_grille[i], tetra.operator[](3))){
// 											volume_s += std::abs(CGAL::to_double(tetra.volume())); 
// 										}
// 										else {
// 											volume_s += intersect_cube_tetrahedron(box_grille[i], tetra);
// 										}
// 									}
// 								}
// 								
// 						}//boucle sur les tetra
// 					} // if inter grille[i] et solide[iter_s] 
// 				} //fin boucle sur les particules			
// 				i++;
//      }
// 		}
// 	}
// 	cout<<"volume solide := "<<volume_s<<endl;
// }
// 
//fin test septembre 2013




void Grille::parois_cellule_vide(Solide& S) {
	
// 	for(int i=0;i<Nx+2*marge;i++){
// 		for(int j=0;j<Ny+2*marge;j++){
// 			for(int k=0;k<Nz+2*marge;k++){
// 				Cellule c = grille[i][j][k]; 
// 				
//  				if(c.x>0.3 && c.x<0.6 &&  c.y>0.18 && c.y<0.47 && c.z>0.4 && c.z<0.6){
// 					c.vide=true;
// 					c.rho=0.;
// 					c.p=0.;
// 					c.impx=0.; c.impy=0.; c.impz=0.; c.rhoE=0.; 
// 					c.u=0.; c.v=0.; c.w=0.;
// 					grille[i][j][k]=c;
// 				}
// 			}
// 		}
// 	}
	
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
	int nb_triangles=0.;
	for(int iter=0; iter<nb_particules; iter++){
		nb_triangles += S.solide[iter].triangles.size();
	}
	
	const double eps_box = 0.1;
	std::vector<Bbox> solide(nb_particules);
	for(int it=0; it<nb_particules; it++){
		solide[it] = Bbox(S.solide[it].min_x-eps_box,S.solide[it].min_y-eps_box, S.solide[it].min_z-eps_box, S.solide[it].max_x+eps_box, S.solide[it].max_y+eps_box, S.solide[it].max_z+eps_box);
	}
	
		
	Cellule cel;
	int i=0;
	for (int a=0; a< nx_m; a++){
		for (int b=0; b< ny_m; b++){
			for (int c=0; c< nz_m; c++){
				cel = grille[a][b][c]; 
				Triangles trianglesB;
				for(int iter_s=0; iter_s<nb_particules; iter_s++){ //boucle sur les particules 
					if (CGAL::do_intersect(box_grille[i],solide[iter_s]) ) {
							triang_cellule(box_grille[i] , trianglesB); 
							for ( int j = 0; j < S.solide[iter_s].triangles.size(); j++){ 
								if (CGAL::do_intersect(box_grille[i],S.solide[iter_s].triangles[j]) ) {
									if(S.solide[iter_s].vide[j] && abs(grille[a][b][c].alpha0 -1.)<eps){
										//cout<<" vide "<< grille[a][b][c].y << endl;  getchar();
										grille[a][b][c].vide = true;
										grille[a][b][c].rho = 0.;
										grille[a][b][c].p = 0.;
										grille[a][b][c].u = grille[a][b][c].v = grille[a][b][c].w = 0.;
										//grille[a][b][c].rhoE = 0.;
										//grille[a][b][c].impx =  grille[a][b][c].impy = grille[a][b][c].impz = 0.;
									}
								} // if intersection cellule[i] avec triangle[j]	 
							} // end boucle sur triangles
					} // if inter grille[i] et solide[iter_s] 
				} //fin boucle sur les particules
				i++;
			} //fin boucle sur grille
		}
	}
	for(int iter_s=0; iter_s<nb_particules; iter_s++){
		for ( int j = 0; j < S.solide[iter_s].faces.size(); j++){ 
			if(S.solide[iter_s].faces[j].voisin == -2)  {S.solide[iter_s].faces[j].voisin = -1;}
		}
	}
	
}

