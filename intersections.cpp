#include "intersections.hpp"
#include "fluide.hpp"
#include "solide.hpp"
#include "solide.cpp"

using std:: cout;
using std:: endl;

void Grille::parois(Solide& S) {
	
	const double eps_relat = numeric_limits<double>::epsilon( );
	
	//cout<<"Eps relative :"<<1- 100/100<<endl;
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
		//S[it].Affiche();
	}
	
	
	
	
	
	Triangles interface;
	
	int count=0;
	double volume_s=0.;
	
	Cellule cel;
	int i=0;
	CGAL::Timer user_time;
	user_time.start();
	for (int a=0; a< nx_m; a++){
		for (int b=0; b< ny_m; b++){
			for (int c=0; c< nz_m; c++){
				cel = grille[a][b][c]; 
				cel.phi_x = 0.;
				cel.phi_y = 0.;
				cel.phi_z = 0.;
				std::vector<Point_3> Points_poly; 
				double alpha = 0.0;
				std::vector<double>  kappa(6,0.0);
				bool intersection = false;
				bool exterieur = true;
				bool box_in_solide = false;
				bool point_in_solide = false;
				bool point_in_cell = false;
				Triangles trianglesB;
				for(int iter_s=0; iter_s<nb_particules; iter_s++){
					if (CGAL::do_intersect(box_grille[i],solide[iter_s]) ) {	
						intersection = true;
						box_in_solide = box_inside_convex_polygon(S.solide[iter_s],box_grille[i]);
						if(box_in_solide){
							exterieur = false;
							cel.alpha = 1.;
							cel.kappai = 1.;
							cel.kappaj = 1.;
							cel.kappak = 1.;
							volume_s += volume_cel;
							box_in_solide = false;
						}
						else 
						{    
							triang_cellule(box_grille[i] , trianglesB); 
							
							for ( int j = 0; j < S.solide[iter_s].triangles.size(); j++){ 
								
								//test if point is in cell_box
								point_in_cell = inside_box(box_grille[i], S.solide[iter_s].triangles[j].operator[](0));
								if(point_in_cell) {Points_poly.push_back(S.solide[iter_s].triangles[j].operator[](0)); point_in_cell=false;}
								
								point_in_cell = inside_box(box_grille[i], S.solide[iter_s].triangles[j].operator[](1));
								if(point_in_cell) {Points_poly.push_back(S.solide[iter_s].triangles[j].operator[](1)); point_in_cell=false;}
								
								point_in_cell = inside_box(box_grille[i], S.solide[iter_s].triangles[j].operator[](2));
								if(point_in_cell) {Points_poly.push_back(S.solide[iter_s].triangles[j].operator[](2)); point_in_cell=false;}
								
								
								if (CGAL::do_intersect(box_grille[i],S.solide[iter_s].triangles[j]) ) {
									for ( int k = 0; k < trianglesB.size(); k++){
										//test if point is in solide
										point_in_solide = inside_convex_polygon(S.solide[iter_s],trianglesB[k].operator[](0));
										if(point_in_solide) {Points_poly.push_back(trianglesB[k].operator[](0)); point_in_solide=false;}
										
										point_in_solide = inside_convex_polygon(S.solide[iter_s],trianglesB[k].operator[](1));
										if(point_in_solide) {Points_poly.push_back(trianglesB[k].operator[](1)); point_in_solide=false;}
										
										point_in_solide = inside_convex_polygon(S.solide[iter_s],trianglesB[k].operator[](2));
										if(point_in_solide) {Points_poly.push_back(trianglesB[k].operator[](2)); point_in_solide=false;}
										
										if (CGAL::do_intersect(S.solide[iter_s].triangles[j],trianglesB[k]) ) {
											
											Triangle_3 t;
											Point_3 P;
											Segment_3 seg;
											std::vector<Point_3> vPoints; 
											CGAL::Object result = CGAL::intersection(S.solide[iter_s].triangles[j], trianglesB[k]);
											
											if(CGAL::assign(P,result)){
												Points_poly.push_back(P);
											}
											else if(CGAL::assign(seg,result)){
												Points_poly.push_back(seg.operator[](0));
												Points_poly.push_back(seg.operator[](1));
											}
											else if(CGAL::assign(t,result)){
												Points_poly.push_back(t.operator[](0));
												Points_poly.push_back(t.operator[](1));
												Points_poly.push_back(t.operator[](2));
												
												
											}
											else if(CGAL::assign(vPoints,result)){ 
												for(int l= 0; l<vPoints.size(); l++)
												{
													Points_poly.push_back(vPoints[l]);
												}
												
											}
											else {cout<<"Intersection type: ? "<< trianglesB[k]<<endl; 
											count++;
											}
										}
										//end intr t1 et t2
									} //end boucle sur trianglesB
									
								} // if intersection cellule[i] avec triangle[j]	 
								
							} // end boucle sur triangles
						} //else 
						
					} // if inter grille[i] et solide[iter_s]  
				} //fin boucle sur les particules
				
				//traitement calcul de alpha et kappa pour la cellule=grille[i]!!!!!!
				if(intersection && exterieur){
					Triangulation T(Points_poly.begin(), Points_poly.end());
					
					std::vector<double>  v_lambda;
					std::vector<Vector_3> v_n_lambda;
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
														 interface.push_back(K);
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
								interface.push_back(K);
							}
							
							else if (abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](0).operator[](2))<=eps_relat && abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](1).operator[](2))<=eps_relat && abs(trianglesB[2].operator[](0).operator[](2) -  K.operator[](2).operator[](2))<=eps_relat)
							{ 
								kappa[1] +=sqrt(CGAL::to_double(K.squared_area()));
								v_lambda.push_back(sqrt(CGAL::to_double(K.squared_area())));
								norm= orthogonal_vector(K.operator[](0),K.operator[](1),K.operator[](2));
								norm2= sqrt(CGAL::to_double(norm*norm));
								v_n_lambda.push_back(norm/norm2);
								interface.push_back(K);
								
							}
							
							else if (abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](0).operator[](0))<=eps_relat && abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](1).operator[](0))<=eps_relat && abs(trianglesB[4].operator[](0).operator[](0) -  K.operator[](2).operator[](0))<=eps_relat)
							{ 
								kappa[2] +=sqrt(CGAL::to_double(K.squared_area()));
								v_lambda.push_back(sqrt(CGAL::to_double(K.squared_area())));
								norm= orthogonal_vector(K.operator[](0),K.operator[](1),K.operator[](2));
								norm2= sqrt(CGAL::to_double(norm*norm));
								v_n_lambda.push_back(norm/norm2);
								interface.push_back(K);
								
							}
							
							else if (abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](0).operator[](0))<=eps_relat && abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](1).operator[](0))<=eps_relat && abs(trianglesB[6].operator[](0).operator[](0) -  K.operator[](2).operator[](0))<=eps_relat)
							{ 
								kappa[3] +=sqrt(CGAL::to_double(K.squared_area()));
								v_lambda.push_back(sqrt(CGAL::to_double(K.squared_area())));
								norm= orthogonal_vector(K.operator[](0),K.operator[](1),K.operator[](2));
								norm2= sqrt(CGAL::to_double(norm*norm));
								v_n_lambda.push_back(norm/norm2);
								interface.push_back(K);
								
							}
							
							else if (abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](0).operator[](1))<=eps_relat && abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](1).operator[](1))<=eps_relat && abs(trianglesB[8].operator[](0).operator[](1) -  K.operator[](2).operator[](1))<=eps_relat)
							{ 
								kappa[4] +=sqrt(CGAL::to_double(K.squared_area()));
								v_lambda.push_back(sqrt(CGAL::to_double(K.squared_area())));
								norm= orthogonal_vector(K.operator[](0),K.operator[](1),K.operator[](2));
								norm2= sqrt(CGAL::to_double(norm*norm));
								v_n_lambda.push_back(norm/norm2);
								interface.push_back(K);
								
							}
							
							else 
							{ 
								kappa[5] +=sqrt(CGAL::to_double(K.squared_area()));
								v_lambda.push_back(sqrt(CGAL::to_double(K.squared_area())));
								norm= orthogonal_vector(K.operator[](0),K.operator[](1),K.operator[](2));
								norm2= sqrt(CGAL::to_double(norm*norm));
								v_n_lambda.push_back(norm/norm2);
								interface.push_back(K);
								
								
							}
						}
						
					}
					
					
					for (int it= 0; it< v_lambda.size(); it++)
					{
						cel.phi_x += v_lambda[it] *( CGAL::to_double(v_n_lambda[it].x()))/volume_cel;
						cel.phi_y += v_lambda[it] *( CGAL::to_double(v_n_lambda[it].y()))/volume_cel;
						cel.phi_z += v_lambda[it] *( CGAL::to_double(v_n_lambda[it].z()))/volume_cel;
					}
					
					if (abs(cel.phi_x)<=eps_relat) {cel.phi_x = 0.;} 
					if (abs(cel.phi_y)<=eps_relat) {cel.phi_y = 0.;} 
					if (abs(cel.phi_z)<=eps_relat) {cel.phi_z = 0.;} 
					
					cel.alpha = alpha/volume_cel;
					cel.kappai = kappa[3]/(deltay * deltaz);
					cel.kappaj = kappa[4]/(deltax * deltaz);
					cel.kappak = kappa[1]/(deltax * deltay);
					volume_s +=alpha;
				}
				grille[a][b][c] = cel;
				i++; 
			} //fin boucle sur grille
		}
	}
	
	
	cout << "Intersection time is: " << user_time.time() << " seconds." << endl;
	cout<<"count triangles degenerate is : "<<count<<endl;
	user_time.reset();
	
	cout<<"volume solide := "<<volume_s<<endl;
	
}
