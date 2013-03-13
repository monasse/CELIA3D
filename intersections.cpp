#include "intersections.hpp"
#include "fluide.hpp"
#include "solide.hpp"
#include "solide.cpp"

using std:: cout;
using std:: endl;

void Grille::parois(Solide& S,double dt) {
	
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
		//S.solide[it].Affiche();
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
// 				grille[a][b][c].alpha = 0.; grille[a][b][c].kappai = 0.; grille[a][b][c].kappaj = 0.; grille[a][b][c].kappak = 0.;
// 				grille[a][b][c].phi_x = 0.; grille[a][b][c].phi_y = 0.; grille[a][b][c].phi_z = 0.; grille[a][b][c].phi_v = 0.;
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
					
					std::vector<double>  v_lambda;
					std::vector<Vector_3> v_n_lambda;
					std::vector<Point_3> X_f; //centre de la face
					///////////////////////////////////////////////////////////////////////////////	  
					if (T.dimension() == 3){
						Polyhedron_3 poly;
						CGAL::convex_hull_3(T.points_begin(), T.points_end(), poly);
						int l=0, n=0;
						double norm2=0.;
						Vector_3 norm;
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
														 
														 v_lambda.push_back(sqrt(CGAL::to_double(K.squared_area())) );
														 norm= orthogonal_vector(K.operator[](0),K.operator[](1),K.operator[](2));
														 norm2= sqrt(CGAL::to_double(norm*norm));
														 v_n_lambda.push_back(norm/norm2);
														 X_f.push_back(centroid(K.operator[](0),K.operator[](1),K.operator[](2)));
													 } //calcul des aires parietales
													 
						}
					}
					
					if (T.dimension() == 2){
						int l=0, n=0;
						double norm2=0.;
						Vector_3 norm;
						Finite_faces_iterator it;
						for (it = T.finite_facets_begin(); it != T.finite_facets_end(); it++){
							Triangle_3 K= T.triangle(*it);
							
							if (abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](0).operator[](2))<=eps_relat && abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](1).operator[](2))<=eps_relat && abs(trianglesB[0].operator[](0).operator[](2) -  K.operator[](2).operator[](2))<=eps_relat)
							{ 
								kappa[0] +=sqrt(CGAL::to_double(K.squared_area()));
								v_lambda.push_back(sqrt(CGAL::to_double(K.squared_area())));
								norm= orthogonal_vector(K.operator[](0),K.operator[](1),K.operator[](2));
								norm2= sqrt(CGAL::to_double(norm*norm));
								v_n_lambda.push_back(norm/norm2);
								X_f.push_back(centroid(K.operator[](0),K.operator[](1),K.operator[](2)));
							}
							
							else if (abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](0).operator[](2))<=eps_relat && abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](1).operator[](2))<=eps_relat && abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](2).operator[](2))<=eps_relat)
							{ 
								kappa[1] +=sqrt(CGAL::to_double(K.squared_area()));
								v_lambda.push_back(sqrt(CGAL::to_double(K.squared_area())));
								norm= orthogonal_vector(K.operator[](0),K.operator[](1),K.operator[](2));
								norm2= sqrt(CGAL::to_double(norm*norm));
								v_n_lambda.push_back(norm/norm2);
								X_f.push_back(centroid(K.operator[](0),K.operator[](1),K.operator[](2)));
							}
							
							else if (abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](0).operator[](0))<=eps_relat && abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](1).operator[](0))<=eps_relat && abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](2).operator[](0))<=eps_relat)
							{ 
								kappa[2] +=sqrt(CGAL::to_double(K.squared_area()));
								v_lambda.push_back(sqrt(CGAL::to_double(K.squared_area())));
								norm= orthogonal_vector(K.operator[](0),K.operator[](1),K.operator[](2));
								norm2= sqrt(CGAL::to_double(norm*norm));
								v_n_lambda.push_back(norm/norm2);
								X_f.push_back(centroid(K.operator[](0),K.operator[](1),K.operator[](2)));
							}
							
							else if (abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](0).operator[](0))<=eps_relat && abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](1).operator[](0))<=eps_relat && abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](2).operator[](0))<=eps_relat)
							{ 
								kappa[3] +=sqrt(CGAL::to_double(K.squared_area()));
								v_lambda.push_back(sqrt(CGAL::to_double(K.squared_area())));
								norm= orthogonal_vector(K.operator[](0),K.operator[](1),K.operator[](2));
								norm2= sqrt(CGAL::to_double(norm*norm));
								v_n_lambda.push_back(norm/norm2);
								X_f.push_back(centroid(K.operator[](0),K.operator[](1),K.operator[](2)));
							}
							
							else if (abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](0).operator[](1))<=eps_relat && abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](1).operator[](1))<=eps_relat && abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](2).operator[](1))<=eps_relat)
							{ 
								kappa[4] +=sqrt(CGAL::to_double(K.squared_area()));
								v_lambda.push_back(sqrt(CGAL::to_double(K.squared_area())));
								norm= orthogonal_vector(K.operator[](0),K.operator[](1),K.operator[](2));
								norm2= sqrt(CGAL::to_double(norm*norm));
								v_n_lambda.push_back(norm/norm2);
								X_f.push_back(centroid(K.operator[](0),K.operator[](1),K.operator[](2)));
							}
							
							else 
							{ 
								kappa[5] +=sqrt(CGAL::to_double(K.squared_area()));
								v_lambda.push_back(sqrt(CGAL::to_double(K.squared_area())));
								norm= orthogonal_vector(K.operator[](0),K.operator[](1),K.operator[](2));
								norm2= sqrt(CGAL::to_double(norm*norm));
								v_n_lambda.push_back(norm/norm2);
								X_f.push_back(centroid(K.operator[](0),K.operator[](1),K.operator[](2)));
							}
						}
						
					}
					
// 					for (int it= 0; it< v_lambda.size(); it++)
// 					{
// 						cel.phi_x += cel.pdtx * v_lambda[it] *( CGAL::to_double(v_n_lambda[it].x()))/volume_cel;
// 						cel.phi_y += cel.pdty *v_lambda[it] *( CGAL::to_double(v_n_lambda[it].y()))/volume_cel;
// 						cel.phi_z += cel.pdtz *v_lambda[it] *( CGAL::to_double(v_n_lambda[it].z()))/volume_cel;
// 						bool in = false; Vector_3 V_f;
// 						for (int iter_s=0; iter_s<nb_particules && in==false; iter_s++){
// 							in = inside_convex_polygon(S.solide[iter_s] ,X_f[it]);
// 							V_f = S.solide[iter_s].vitesse_parois(X_f[it]);
// 						}
//             cel.phi_v += v_lambda[it] * (CGAL::to_double(cel.pdtx*v_n_lambda[it].x()*V_f.x()  + cel.pdty*v_n_lambda[it].y()*V_f.y()+
//                                           cel.pdtz*v_n_lambda[it].z()*V_f.z()))/volume_cel;
//              }
// 					
// 					if (abs(cel.phi_x)<=eps_relat) {cel.phi_x = 0.;} 
// 					if (abs(cel.phi_y)<=eps_relat) {cel.phi_y = 0.;} 
// 					if (abs(cel.phi_z)<=eps_relat) {cel.phi_z = 0.;} 
// 					if (abs(cel.phi_v)<=eps_relat) {cel.phi_v = 0.;} 
					
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
				if(std::abs(grille[a][b][c].alpha -1.) <eps) {
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
