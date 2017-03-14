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
  \brief Intersection of the fluid grid with the solid. 
  \warning  <b> Specific coupling procedures ! </b>
*/

#include "intersections.hpp"
#include "fluide.hpp"
#include "solide.hpp"
#include "solide.cpp"

using std:: cout;
using std:: endl;






/*!\brief Intersection of a bounding box with a tetrahedron. 
  \details Intersection of the cubic bounding box \b cube with tetrahedron \b Tet. Returns the volume of the intersection. This function is called in the swept quantity calculation. \n
  Algorithm: \n
  - Triangulate the faces of \b cube using function triang_cellule(const Bbox&, Triangles& ). \n
  - Search for the vertices of the tetrahedron contained in \a cube using function inside_box(const Bbox&, const Point_3& ). \n
  - Loop on the tetrahedron faces.
  - Test the intersection between \a cube and the tetrahedron faces using function \b  CGAL::do_intersect(Bbox, Triangle_3). Si oui:
  - Loop on the triangular faces of \a cube.
  - Test the intersection between faces of \a cube and faces of the tetrahedron using function  \b CGAL::do_intersect(Triangle_3, Triangle_3). If so: 
  - Compute the intersection of the faces of the tetrahedron with the faces of \a cube using function
  \b CGAL::do_intersect(Triangle_3, Triangle_3).
  - Compute the volume of the polyhedron which results from the intersection of \a cube and \a Tet. For the volume computation, construct with the intersection points a tetrahedron tesselation using function \b  CGAL::Triangulation(vector<Point_3>). The volume is the sum of the volumes of the tetrahedra, computed using function \b CGAL::tetrahedron.volume().

  \param cube Box 3d 
  \param Tet Tetrahedron
  \warning <b> Specific coupling procedure ! </b>
  \return double
*/
double intersect_cube_tetrahedron(const Bbox& cube, const Tetrahedron& Tet, double& temps_intersections, double& temps_triangulation){
	
  double volume=0.,volume2=0.;
  CGAL::Timer intersections_time,triangulation_time;
  intersections_time.start();
  

  if(abs(Tet.volume())>eps){
    
    Triangles triangCube;
    triang_cellule(cube , triangCube); 
    Triangle_3 triangTet[4];
    std::vector<Point_3> Points_intersect;//, Points_intersect2;
	  
    //Vertices of the tetrahedron
    for(int i=0; i<4; i++){
      if (CGAL::do_overlap(cube, Tet.operator[](i).bbox()) ){
	Points_intersect.push_back(Tet.operator[](i));
      }
    }
		
    triangTet[0]= Triangle_3(Tet.operator[](0), Tet.operator[](1), Tet.operator[](2));
    triangTet[1]= Triangle_3(Tet.operator[](0), Tet.operator[](2), Tet.operator[](3));
    triangTet[2]= Triangle_3(Tet.operator[](0), Tet.operator[](1), Tet.operator[](3));
    triangTet[3]= Triangle_3(Tet.operator[](1), Tet.operator[](2), Tet.operator[](3));
	
    //Cube vertices
    for(int kx=0;kx<2;kx++){
      for(int ky=0;ky<2;ky++){
	for(int kz=0;kz<2;kz++){
	  Point_3 s1(cube.xmin()+kx*(cube.xmax()-cube.xmin()),cube.ymin()+ky*(cube.ymax()-cube.ymin()),cube.zmin()+kz*(cube.zmax()-cube.zmin()));
		
	  if(CGAL::do_overlap(s1.bbox(),Tet.bbox())){
	    if(inside_tetra(Tet,s1)){
	      Points_intersect.push_back(s1);
	    }
	  }
	}
      }
    }
	  
    //Cube edges
    for(int kx=0;kx<2;kx++){
      for(int ky=0;ky<2;ky++){
	for(int kz=0;kz<2;kz++){
	  Point_3 s1(cube.xmin()+kx*(cube.xmax()-cube.xmin()),cube.ymin()+ky*(cube.ymax()-cube.ymin()),cube.zmin()+kz*(cube.zmax()-cube.zmin()));
	  if(kx==0){
	    Point_3 s2(cube.xmax(),s1.y(),s1.z());
	    Segment_3 seg(s1,s2);
	    if(CGAL::do_overlap(seg.bbox(),Tet.bbox())){
	      for(int f=0;f<4;f++){
		if(CGAL::do_intersect(seg,triangTet[f])){
		  Point_3 P;
		  Segment_3 a;
		  CGAL::Object result = CGAL::intersection(seg,triangTet[f]);
		  if(assign(P,result)){
		    Points_intersect.push_back(P);
		  }
		  else if(assign(a,result)){
		    Points_intersect.push_back(a.vertex(0));
		    Points_intersect.push_back(a.vertex(1));
		  } else {
		    cout << "Problem edge/triangle intersection" << endl;
		    getchar();
		  }
		}
	      }
	    }
	  }
	  if(ky==0){
	    Point_3 s2(s1.x(),cube.ymax(),s1.z());
	    Segment_3 seg(s1,s2);
	    if(CGAL::do_overlap(seg.bbox(),Tet.bbox())){
	      for(int f=0;f<4;f++){
		if(CGAL::do_intersect(seg,triangTet[f])){
		  Point_3 P;
		  Segment_3 a;
		  CGAL::Object result = CGAL::intersection(seg,triangTet[f]);
		  if(assign(P,result)){
		    Points_intersect.push_back(P);
		  }
		  else if(assign(a,result)){
		    Points_intersect.push_back(a.vertex(0));
		    Points_intersect.push_back(a.vertex(1));
		  } else {
		    cout << "Problem edge/triangle intersection" << endl;
		    getchar();
		  }
		}
	      }
	    }
	  }
	  if(kz==0){
	    Point_3 s2(s1.x(),s1.y(),cube.zmax());
	    Segment_3 seg(s1,s2);
	    if(CGAL::do_overlap(seg.bbox(),Tet.bbox())){
	      for(int f=0;f<4;f++){
		if(CGAL::do_intersect(seg,triangTet[f])){
		  Point_3 P;
		  Segment_3 a;
		  CGAL::Object result = CGAL::intersection(seg,triangTet[f]);
		  if(assign(P,result)){
		    Points_intersect.push_back(P);
		  }
		  else if(assign(a,result)){
		    Points_intersect.push_back(a.vertex(0));
		    Points_intersect.push_back(a.vertex(1));
		  } else {
		    cout << "Problem edge/triangle intersection" << endl;
		    getchar();
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    //Intersection of the tetrahedron edges with the cube faces
    Segment_3 areteTet[6];
    areteTet[0] = Segment_3(Tet.vertex(0),Tet.vertex(1));
    areteTet[1] = Segment_3(Tet.vertex(0),Tet.vertex(2));
    areteTet[2] = Segment_3(Tet.vertex(0),Tet.vertex(3));
    areteTet[3] = Segment_3(Tet.vertex(1),Tet.vertex(2));
    areteTet[4] = Segment_3(Tet.vertex(1),Tet.vertex(3));
    areteTet[5] = Segment_3(Tet.vertex(2),Tet.vertex(3));
	  
    Triangles trianglesB;
    triang_cellule(cube,trianglesB);
    for(int a=0;a<6;a++){
      for(int t=0;t<trianglesB.size();t++){
	if(CGAL::do_overlap(areteTet[a].bbox(),trianglesB[t].bbox())){
	  if(CGAL::do_intersect(areteTet[a],trianglesB[t])){
	    Point_3 P;
	    Segment_3 seg;
	    CGAL::Object result = CGAL::intersection(areteTet[a],trianglesB[t]);
	    if(assign(P,result)){
	      Points_intersect.push_back(P);
	    }
	    else if(assign(seg,result)){
	      Points_intersect.push_back(seg.vertex(0));
	      Points_intersect.push_back(seg.vertex(1));
	    } else {
	      cout << "Problem edge/triangle intersection" << endl;
	      getchar();
	    }
	  }
	}
      }
    }
    temps_intersections += intersections_time.time();
	  
    triangulation_time.start();
	  	  
    //Triangulation of the intersection
    if( Points_intersect.size() >=4 ){
      Exact_to_Inexact to_inexact;
      std::vector<InexactPoint_3> Points_intersect3;
      for(std::vector<Point_3>::iterator it = Points_intersect.begin();it!=Points_intersect.end();it++){
	Points_intersect3.push_back(to_inexact(*it));
      }
      if(Points_intersect3.size()>=4){
	      
	InexactTriangulation T(Points_intersect3.begin(), Points_intersect3.end());
	if(T.is_valid()){
	  InexactFinite_cells_iterator cit;
	  for (cit = T.finite_cells_begin(); cit!= T.finite_cells_end(); cit++){
	    volume+= std::abs(CGAL::to_double(T.tetrahedron( cit).volume()));
	  }
	}
      }
    }
    temps_triangulation += triangulation_time.time();
	  
    return std::abs(volume);
  } else {
    return 0.;
  }
  
	
}	

/*!\brief Alternative version of the intersection of a cube with a tetrahedron
 */
double intersect_cube_tetrahedron_bis(const Bbox& cube, const Tetrahedron& Tet, double& temps_intersections, double& temps_triangulation){
  double volume = 0.;
  if(abs(Tet.volume())>eps){
    //Intersection test
    bool intersect=false;
    //Test on the tetrahedron vertices
    for(int i=0; i<4 && !intersect; i++){
      if (CGAL::do_overlap(cube, Tet.vertex(i).bbox()) ){
	intersect = true;
      }
    }
    //Test on the cube vertices
    for(int kx=0;kx<2 && !intersect;kx++){
      for(int ky=0;ky<2 && !intersect;ky++){
	for(int kz=0;kz<2 && !intersect;kz++){
	  Point_3 s1(cube.xmin()+kx*(cube.xmax()-cube.xmin()),cube.ymin()+ky*(cube.ymax()-cube.ymin()),cube.zmin()+kz*(cube.zmax()-cube.zmin()));
	  
	  if(CGAL::do_overlap(s1.bbox(),Tet.bbox())){
	    if(inside_tetra(Tet,s1)){
	      intersect = true;
	    }
	  }
	}
      }
    }
    //Test on the cube edges
    if(!intersect){
      Triangle_3 triangTet[4];
      triangTet[0]= Triangle_3(Tet.vertex(0), Tet.vertex(1), Tet.vertex(2));
      triangTet[1]= Triangle_3(Tet.vertex(0), Tet.vertex(2), Tet.vertex(3));
      triangTet[2]= Triangle_3(Tet.vertex(0), Tet.vertex(1), Tet.vertex(3));
      triangTet[3]= Triangle_3(Tet.vertex(1), Tet.vertex(2), Tet.vertex(3));
      for(int kx=0;kx<2 && !intersect;kx++){
	for(int ky=0;ky<2 && !intersect;ky++){
	  for(int kz=0;kz<2 && !intersect;kz++){
	    Point_3 s1(cube.xmin()+kx*(cube.xmax()-cube.xmin()),cube.ymin()+ky*(cube.ymax()-cube.ymin()),cube.zmin()+kz*(cube.zmax()-cube.zmin()));
	    if(kx==0){
	      Point_3 s2(cube.xmax(),s1.y(),s1.z());
	      Segment_3 seg(s1,s2);
	      if(CGAL::do_overlap(seg.bbox(),Tet.bbox())){
		for(int f=0;f<4;f++){
		  if(CGAL::do_intersect(seg,triangTet[f])){
		    intersect = true;
		  }
		}
	      }
	    }
	    if(ky==0){
	      Point_3 s2(s1.x(),cube.ymax(),s1.z());
	      Segment_3 seg(s1,s2);
	      if(CGAL::do_overlap(seg.bbox(),Tet.bbox())){
		for(int f=0;f<4;f++){
		  if(CGAL::do_intersect(seg,triangTet[f])){
		    intersect = true;
		  }
		}
	      }
	    }
	    if(kz==0){
	      Point_3 s2(s1.x(),s1.y(),cube.zmax());
	      Segment_3 seg(s1,s2);
	      if(CGAL::do_overlap(seg.bbox(),Tet.bbox())){
		for(int f=0;f<4;f++){
		  if(CGAL::do_intersect(seg,triangTet[f])){
		    intersect = true;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    //Test the intersection of the tetrahedron edges with the cube faces
    if(!intersect){
      Segment_3 areteTet[6];
      areteTet[0] = Segment_3(Tet.vertex(0),Tet.vertex(1));
      areteTet[1] = Segment_3(Tet.vertex(0),Tet.vertex(2));
      areteTet[2] = Segment_3(Tet.vertex(0),Tet.vertex(3));
      areteTet[3] = Segment_3(Tet.vertex(1),Tet.vertex(2));
      areteTet[4] = Segment_3(Tet.vertex(1),Tet.vertex(3));
      areteTet[5] = Segment_3(Tet.vertex(2),Tet.vertex(3));
      
      Triangles trianglesB;
      triang_cellule(cube,trianglesB);
      for(int a=0;a<6 && !intersect;a++){
	for(int t=0;t<trianglesB.size() && !intersect;t++){
	  if(CGAL::do_overlap(areteTet[a].bbox(),trianglesB[t].bbox())){
	    if(CGAL::do_intersect(areteTet[a],trianglesB[t])){
	      intersect = true;
	    }
	  }
	}
      }
    }
    
    //In case there is an intersection
    if(intersect){
      //Translation with regards to the tetrahedron center
      Point_3 centre = CGAL::centroid(Tet.vertex(0), Tet.vertex(1), Tet.vertex(2),Tet.vertex(3));
      if(!CGAL::do_overlap(centre.bbox(),cube)){
	//TO DO !
      }
      Vector_3 vec(centre,Point_3(0,0,0));
      //Projective transformation of planes
      Plane_3 plan[10];
      //Tetrahedron planes
      plan[0]= Plane_3(Tet.vertex(0)+vec, Tet.vertex(1)+vec, Tet.vertex(2)+vec);
      plan[1]= Plane_3(Tet.vertex(0)+vec, Tet.vertex(2)+vec, Tet.vertex(3)+vec);
      plan[2]= Plane_3(Tet.vertex(0)+vec, Tet.vertex(1)+vec, Tet.vertex(3)+vec);
      plan[3]= Plane_3(Tet.vertex(1)+vec, Tet.vertex(2)+vec, Tet.vertex(3)+vec);
      //Cube planes
      plan[4] = Plane_3(1,0,0,-cube.xmin()-vec.x());
      plan[5] = Plane_3(1,0,0,-cube.xmax()-vec.x());
      plan[6] = Plane_3(0,1,0,-cube.ymin()-vec.y());
      plan[7] = Plane_3(0,1,0,-cube.ymax()-vec.y());
      plan[8] = Plane_3(0,0,1,-cube.zmin()-vec.z());
      plan[9] = Plane_3(0,0,1,-cube.zmax()-vec.z());
      //Projective transformation
      std::vector<Point_3> dual;
      for(int i=0;i<10;i++){
	if(plan[i].d()==0){
	  cout << "Problem dual construction d()=" << plan[i].d() << " i=" << i << endl;
	  getchar();
	}
	dual.push_back(Point_3(plan[i].a(),plan[i].b(),plan[i].c(),plan[i].d()));
      }	
      //Construction of the dual convex hull
      if(!coplanar(dual.begin(),dual.end())){
	Polyhedron_3 poly_dual;
	CGAL::convex_hull_3(dual.begin(), dual.end(), poly_dual);
	
	//Computation of the volume
	//Browse the primal faces (dual vertices)
	Plane_3 plane_0(poly_dual.vertices_begin()->vertex_begin()->vertex()->point(),poly_dual.vertices_begin()->vertex_begin()->next()->vertex()->point(),poly_dual.vertices_begin()->vertex_begin()->opposite()->vertex()->point());
	//cout << "plane_0=" << plane_0 << endl;
	Point_3 P0(plane_0.a(),plane_0.b(),plane_0.c(),plane_0.d());
	for(Vertex_iterator v=poly_dual.vertices_begin();v!=poly_dual.vertices_end();v++){
	  Point_3 P_dual = v->point();
	  //Browse the primal vertices (faces around the dual vertex) and store them in sommets
	  std::vector<Point_3> sommets;
	  Halfedge_around_vertex_circulator he=(*v).vertex_begin(), done(he);
	  do 
	  {
	    Plane_3 plane_dual(he->vertex()->point(),he->opposite()->vertex()->point(),he->next()->vertex()->point());
	    if(plane_dual.d()==0){
	      cout << "Problem recovery of primal d()=" << plane_dual.d() << endl;
	      getchar();
	    }
	    sommets.push_back(Point_3(plane_dual.a(),plane_dual.b(),plane_dual.c(),plane_dual.d()));
	    Point_3 Pt(plane_dual.a(),plane_dual.b(),plane_dual.c(),plane_dual.d());
	  } while(++he != done);
	  for(int i=0;i<sommets.size();i++){
	    cout << "Sommet_bis=" << sommets[i]-vec << endl;
	  }
	  for(int i=1;i<sommets.size()-1;i++){
	    Tetrahedron T(P0,sommets[0],sommets[i],sommets[i+1]);
	    volume += std::abs(CGAL::to_double(T.volume()));
	  }
	}
      }
    }
    
  }
  return volume;
  
}





/*!\brief Intersection of the fluid grid with solid.
  \details Intersection of the fluid grid with the solid and computation of the quantities of interest: solid occupancy ratio in the cell (\a Cellule.alpha), solid occupancy ratio on the cell faces (\a Cellule.kappai, \a Cellule.kappaj and \a Cellule.kappak). Definition of the interface objects: \n
  - \a Particule.Points_interface: intersection points of the cell with the triangular faces of the solid; \n
  - \a Particule.Triangles_interface: partition of the solid faces into interface triangles contained in one single cell of the fluid grid; \n
  - \a Particule.Position_Triangles_interface: index of the fluid grid cell containing \a Particule.Triangles_interface. \n
 
  Algorithm: \n
  - Construct vector \a box_grille containing the cubic cells of the fluid grid in the form of 3d bounding boxes (\a Bbox). 
  - Construct the vector \a solide of the bounding boxes associated with the \a Particule. 
  - Loop on \a box_grille.
  - Loop on \a solide.
  - Test the intersection between \a box_grille and \a solide using function \b CGAL::do_overlap(Bbox, Bbox). If so:
  - Test whether \a box_grille is fully contained in \a solide via la fonction box_inside_convex_polygon(const Particule&, const Bbox&). If so, the intersection is \a box_grille. Otherwise:
  - Browse the solid vertices contained in \a box_grille
  - Intersect the edges of \a box_grille with the solid faces
  - Loop on the triangular faces of the Solide.
  - Browse the vertices of the triangular faces of Solide contained in \a box_grille using function inside_box(const Bbox&, const Point_3&). 
  - Test the intersections between \a box_grille and the triangular faces of Solide using function \b CGAL::do_intersect(Bbox,Triangle_3). If so:
  - Triangulate the faces of \a box_grille  using function triang_cellule(const Bbox&, Triangles&) .
  - Loop on the triangular faces of \a box_grille.
  - Browse the vertices of the faces of \a box_grille contained in Solide using function inside_convex_polygon(const Particule&, const Point_3&). 
  - Test the intersection between the faces of \a box_grille and the triangular faces of Solide (\a Particule.triangle) using function \b CGAL::do_intersect(Triangle_3, Triangle_3). If so:
  - Compute the intersection between the triangular faces of \a box_grille and the triangular faces of Solide using function \b  CGAL::intersection(Triangle_3, Triangle_3).\n
      
      
      
      
  \remark The vector \a Particule.Points_interface containing the intersection points between fluid cells and triangular faces (\a Particule.triangles) is progressively filled during the intersection algorithm. \n
  The intersection result is used  to compute the interest quantities \a Cellule.alpha, \a Cellule.kappai, \a Cellule.kappaj and \a Cellule.kappak. \n
  The function used in the intersection computation are
  <b> CGAL::Triangulation(vector<Point_3>) </b>, <b> CGAL::convex_hull_3(vector<Point_3>, Polyhedron_3) </b>, <b> CGAL::tetrahedron.volume() </b> and <b> CGAL::Triangle_3.squared_area()</b>. \n
  *\param S Solide
  *\param dt Time-step
  *\warning <b> Specific coupling procedure ! </b>
  *\return void 
  */

void Grille::Parois_particles(Solide& S,double dt) {
  CGAL::Timer total_time,bbox_time,do_intersect_time,triangularisation_time,test_time,test_inside_time,alpha_time,intersect_time,convex_hull_time,volume_time,kappa_time1,kappa_time2,triangulation_time,vertices_time,coin_time,sommet_time,sommet_interface_time,aretes_cellule_time,aretes_solide_time,test_sommet_interface_time,push_back_time,triangulation_time2;
  total_time.start();bbox_time.start();do_intersect_time.start();triangularisation_time.start();test_time.start();test_inside_time.start();alpha_time.start();intersect_time.start();convex_hull_time.start();volume_time.start();kappa_time1,kappa_time2,triangulation_time.start();vertices_time.start();coin_time.start();sommet_time.start();sommet_interface_time.start();aretes_cellule_time.start();aretes_solide_time.start();test_sommet_interface_time.start();push_back_time.start();triangulation_time2.start();
  double temps_total=0.,temps_bbox=0.,temps_do_intersect=0.,temps_triangularisation=0.,temps_test=0.,temps_test_inside=0.,temps_alpha=0.,temps_intersect=0.,nb_intersect=0.,temps_convex_hull=0.,nb_convex_hull=0.,temps_volume=0.,temps_kappa1=0.,temps_kappa2=0.,nb_kappa1=0.,nb_kappa2=0.,temps_triangulation=0.,temps_vertices=0.,temps_coin=0.,temps_sommet=0.,temps_sommet_interface=0.,temps_aretes_cellule=0.,temps_aretes_solide=0.,temps_test_sommet_interface=0.,temps_push_back=0.,nb_sommet_interface=0.,nb_sommet=0.,temps_triangulation2=0.;
  
  total_time.reset();
  bbox_time.reset();
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
  cout<<"Number of particles: "<<nb_particules<<endl;
  int nb_triangles=0.;
  for(int iter=0; iter<nb_particules; iter++){
    nb_triangles += S.solide[iter].triangles.size();
  }
  cout<<"Number of triangles: "<<nb_triangles<<endl;

  vertices_time.reset();
  const double eps_box = 0.1;
  for(int it=0; it<nb_particules; it++){
    //Update the particle vertices
    Particule& P = S.solide[it];
    P.vertices.clear();
    for(std::vector<Face>::iterator f=P.faces.begin();f!=P.faces.end();f++){
      for(std::vector<Vertex>::iterator v=(*f).vertex.begin();v!=(*f).vertex.end();v++){
	bool test_new = true;
	Point_3 p1 = (*v).pos;
	for(std::vector<Point_3>::iterator pit=P.vertices.begin();pit!=P.vertices.end() && test_new;pit++){
	  test_new = ((*pit)!=p1);
	}
	if(test_new){
	  P.vertices.push_back(p1);
	}
      }
    }
		
		
  }
  temps_vertices += CGAL::to_double(vertices_time.time());
	
  temps_bbox += CGAL::to_double(bbox_time.time());
	

  double volume_s=0.;
	
  Cellule cel;
  int i=0;
  CGAL::Timer user_time, user_time2;
  user_time.reset();
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
	do_intersect_time.reset();
	for(int iter_s=0; iter_s<nb_particules && exterieur; iter_s++){ 
	  test_time.reset();
	  bool test = CGAL::do_overlap(box_grille[i],S.solide[iter_s].bbox);
	  temps_test += CGAL::to_double(test_time.time());
	  if (test){
	    std::vector<Point_3> Points_poly; 
	    double alpha = 0.0;
	    std::vector<double>  kappa(6,0.0);
	    bool intersection = false;
	    bool box_in_solide = false;
	    bool point_in_solide = false;
	    bool point_in_cell = false;
	    intersection = true;
	    test_inside_time.reset();
	    box_in_solide = box_inside_convex_polygon(S.solide[iter_s],box_grille[i]);
	    temps_test_inside += CGAL::to_double(test_inside_time.time());
	    if(box_in_solide){
	      exterieur = false;
	      cel.alpha = 1.;
	      cel.kappai = 1.;
	      cel.kappaj = 1.;
	      cel.kappak = 1.;
	      if(a>0) {grille[a-1][b][c].kappai = 1.;}
	      if(b>0) {grille[a][b-1][c].kappaj = 1.;}
	      if(c>0) {grille[a][b][c-1].kappak = 1.;}
	      volume_s += volume_cel;
	      box_in_solide = false;
	    }
	    else 
	    {    
	      triang_cellule(box_grille[i] , trianglesB); 
					
	      //Test whether the cell vertices are inside the solid
	      coin_time.reset();
	      for(int kx=0;kx<2;kx++){
		for(int ky=0;ky<2;ky++){
		  for(int kz=0;kz<2;kz++){
		    double x = box_grille[i].xmin()+kx*(box_grille[i].xmax()-box_grille[i].xmin());
		    double y = box_grille[i].ymin()+ky*(box_grille[i].ymax()-box_grille[i].ymin());
		    double z = box_grille[i].zmin()+kz*(box_grille[i].zmax()-box_grille[i].zmin());
		    point_in_solide = inside_convex_polygon(S.solide[iter_s],Point_3(x,y,z));
		    if(point_in_solide) {Points_poly.push_back(Point_3(x,y,z)); point_in_solide=false;
		    }
		  }
		}
	      }
	      temps_coin += CGAL::to_double(coin_time.time());
				      
	      //Test whether the solid vertices are in the cell
	      sommet_time.reset();
	      for(std::vector<Point_3>::iterator sommet=S.solide[iter_s].vertices.begin();sommet!=S.solide[iter_s].vertices.end();sommet++){
		nb_sommet +=1.;
		if(CGAL::do_overlap(box_grille[i],(*sommet).bbox())){Points_poly.push_back((*sommet));
		}
	      }
	      temps_sommet += sommet_time.time();
				      
				      
				      
	      for ( int j = 0; j < S.solide[iter_s].triangles.size(); j++){ 
					
		if (CGAL::do_overlap(box_grille[i],S.solide[iter_s].triangles[j].bbox()) ) {
					  
		  sommet_interface_time.reset();
		  nb_sommet_interface +=1.;
					  
		  if(CGAL::do_overlap(box_grille[i], S.solide[iter_s].triangles[j].operator[](0).bbox())){
		    Points_interface[iter_s][j].push_back(S.solide[iter_s].triangles[j].operator[](0));
		  }
					  
					  
		  if(CGAL::do_overlap(box_grille[i], S.solide[iter_s].triangles[j].operator[](1).bbox())){
		    Points_interface[iter_s][j].push_back(S.solide[iter_s].triangles[j].operator[](1));
		  }
		
		  if(CGAL::do_overlap(box_grille[i], S.solide[iter_s].triangles[j].operator[](2).bbox())){
		    Points_interface[iter_s][j].push_back(S.solide[iter_s].triangles[j].operator[](2));
		  }
		  temps_sommet_interface += sommet_interface_time.time();
		  			  
					  			  
		  //Intersections of the cell edges with the triangle
		  aretes_cellule_time.reset();
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
			    intersect_time.reset();
			    nb_intersect+=1.;
			    std::vector<Point_3> result = intersection_bis(seg,S.solide[iter_s].triangles[j]);
			    temps_intersect += CGAL::to_double(intersect_time.time());
			    for(int l= 0; l<result.size(); l++)
			    {
			      Points_poly.push_back(result[l]);
			      Points_interface[iter_s][j].push_back(result[l]);
			    }
			  }
			}
			if(ky==0){
			  double y2 = box_grille[i].ymax();
			  Segment_3 seg(Point_3(x1,y1,z1),Point_3(x1,y2,z1));
			  if (CGAL::do_intersect(seg,S.solide[iter_s].triangles[j]) ) {
			    intersect_time.reset();
			    nb_intersect+=1.;
			    std::vector<Point_3> result = intersection_bis(seg,S.solide[iter_s].triangles[j]);
			    temps_intersect += CGAL::to_double(intersect_time.time());
						    
			    for(int l= 0; l<result.size(); l++)
			    {
			      Points_poly.push_back(result[l]);
			      Points_interface[iter_s][j].push_back(result[l]);
			    }
			  }
			}
			if(kz==0){
			  double z2 = box_grille[i].zmax();
			  Segment_3 seg(Point_3(x1,y1,z1),Point_3(x1,y1,z2));
			  if (CGAL::do_intersect(seg,S.solide[iter_s].triangles[j]) ) {
			    intersect_time.reset();
			    nb_intersect+=1.;
			    std::vector<Point_3> result = intersection_bis(seg,S.solide[iter_s].triangles[j]);
			    temps_intersect += CGAL::to_double(intersect_time.time());
						    
			    for(int l= 0; l<result.size(); l++)
			    {
			      Points_poly.push_back(result[l]);
			      Points_interface[iter_s][j].push_back(result[l]);
			    }
			  }
			}
		      }
		    }
		  }
		  temps_aretes_cellule += aretes_cellule_time.time();
					  
					  
		  //Intersection of the solid face edges with the cell triangle faces
		  aretes_solide_time.reset();
		  for ( int k = 0; k < trianglesB.size(); k++){	
		    if (CGAL::do_overlap(S.solide[iter_s].triangles[j].bbox(),trianglesB[k].bbox()) ) {
		      for(int l=0;l<3;l++){
			int lp = (l+1)%3;
			Segment_3 seg(S.solide[iter_s].triangles[j].operator[](l),S.solide[iter_s].triangles[j].operator[](lp));
			if (CGAL::do_intersect(seg,trianglesB[k]) ) {
			  intersect_time.reset();
			  nb_intersect+=1.;
			  std::vector<Point_3> result = intersection_bis(seg,trianglesB[k]);
			  temps_intersect += CGAL::to_double(intersect_time.time());
						  
			  for(int lt= 0; lt<result.size(); lt++)
			  {
			    if(S.solide[iter_s].triangles[j].operator[](l)<=S.solide[iter_s].triangles[j].operator[](lp)){
			      Points_poly.push_back(result[lt]);
			    }
			    Points_interface[iter_s][j].push_back(result[lt]);
			  }
			}
		      }
		    }
		  } 
		  temps_aretes_solide += aretes_solide_time.time();
		  			  
					  
		} 	 
					
	      } 
	    }  
						
					
					
	    //Computation of alpha and kappa for the cell grille[i]
	    user_time2.reset();
	    alpha_time.reset();
	    if(intersection && exterieur && Points_poly.size()>3){
	      Exact_to_Inexact to_inexact;
	      std::vector<InexactPoint_3> Points_poly2;
	      for(std::vector<Point_3>::iterator it = Points_poly.begin();it!=Points_poly.end();it++){
		Points_poly2.push_back(to_inexact(*it));
	      }
	      InexactPolyhedron_3 poly;
	      convex_hull_time.reset();
	      CGAL::convex_hull_3(Points_poly2.begin(), Points_poly2.end(), poly);
		
	      nb_convex_hull += 1.;
	      temps_convex_hull += CGAL::to_double(convex_hull_time.time());
	      volume_time.reset();
	      InexactPoint_3 P = (*(poly.facets_begin())).halfedge()->vertex()->point();
	      InexactFacet_iterator fiter;
	      for (fiter = poly.facets_begin(); fiter!= poly.facets_end(); fiter++){
		InexactTetrahedron T(P,(*fiter).halfedge()->vertex()->point(),(*fiter).halfedge()->next()->vertex()->point(), (*fiter).halfedge()->opposite()->vertex()->point());
		alpha+= CGAL::to_double(T.volume());
	      }
					
					
	      temps_volume += CGAL::to_double(volume_time.time());
					
	      kappa_time1.reset();
	      for (fiter = poly.facets_begin(); fiter!= poly.facets_end(); fiter++){
		nb_kappa1 +=1.;
		InexactTriangle_3 K((*fiter).halfedge()->vertex()->point(),(*fiter).halfedge()->next()->vertex()->point(),
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
															 
		} 
														 
	      }
	      temps_kappa1 += CGAL::to_double(kappa_time1.time());
	      					
	      cel.alpha  += alpha/volume_cel;
	      cel.kappai += kappa[3]/(deltay * deltaz);
	      cel.kappaj += kappa[4]/(deltax * deltaz);
	      cel.kappak += kappa[1]/(deltax * deltay);
	      if(cel.kappai >=1.) {cel.kappai=1.;}
	      if(cel.kappaj >=1.) {cel.kappaj=1.;}
	      if(cel.kappak >=1.) {cel.kappak=1.;}
	      if(cel.alpha >=1.) {cel.alpha=1.;}
	      					
	      volume_s +=alpha;
	      time+= user_time2.time();
	      user_time2.reset();
						
	    }
	    temps_alpha += CGAL::to_double(alpha_time.time());
	    			    
	  }				
					
	} 
	temps_do_intersect += CGAL::to_double(do_intersect_time.time());
	grille[a][b][c] = cel;
	if(std::abs(grille[a][b][c].alpha -1.) <1.e-10) {
	  grille[a][b][c].alpha = 1.;
	  grille[a][b][c].kappai = 1.;
	  grille[a][b][c].kappaj = 1.;
	  grille[a][b][c].kappak = 1.;
	  if(a>0){grille[a-1][b][c].kappai = 1.;}
	  if(b>0){grille[a][b-1][c].kappaj = 1.;}
	  if(c>0){grille[a][b][c-1].kappak = 1.;}
	}
	i++;
				
	//Triangulation of the interface face by face
	triangularisation_time.reset();
	ExactFinite_faces_iterator iter;
	for(int count=0; count<nb_particules;count++){
	  for(int it=0; it<S.solide[count].triangles.size(); it++){ 
			
				    
	    if(Points_interface[count][it].size()>2){ 
	      triangulation_time2.reset();
	      ExactTriangulation T(Points_interface[count][it].begin(), Points_interface[count][it].end());
	      temps_triangulation2 += triangulation_time2.time();
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
	    }
	  } 
	} 
	temps_triangularisation += CGAL::to_double(triangularisation_time.time());
	
      } 
    }
  }	
  user_time.reset();
  cout<<"volume solide parois := "<<volume_s<<endl;
  temps_total = CGAL::to_double(total_time.time());
	
  cout << "######### COUTS INTERSECTIONS ##########" << endl;
  cout << "Bbox=" << 100*temps_bbox/temps_total << "%" << endl;
  cout << "tri vertices=" << 100*temps_vertices/temps_total << "%" << endl;
  cout << "do_intersect=" << 100*temps_do_intersect/temps_total << "%" << endl;
  cout << "   test intersect=" << 100*temps_test/temps_total << "%" << endl;
  cout << "   test_inside=" << 100*temps_test_inside/temps_total << "%" << endl;
  cout << "   coin=" << 100*temps_coin/temps_total << "%" << endl;
  cout << "   sommet=" << 100*temps_sommet/temps_total << "%          t_moy=" << temps_sommet/nb_sommet << " nb_sommet=" << nb_sommet << endl;
  cout << "   sommet_interface=" << 100*temps_sommet_interface/temps_total << "%          t_moy=" << temps_sommet_interface/nb_sommet_interface << " nb_sommet_interface=" << nb_sommet_interface << endl;
  cout << "      test_sommet_interface=" << 100*temps_test_sommet_interface/temps_total << "%" << endl;
  cout << "      push_back=" << 100*temps_push_back/temps_total << "%" << endl;
  cout << "   aretes_cellule=" << 100*temps_aretes_cellule/temps_total << "%" << endl;
  cout << "   aretes_solide=" << 100*temps_aretes_solide/temps_total << "%" << endl;
  cout << "      intersect=" << 100*temps_intersect/temps_total << "%          t_moy=" << temps_intersect/nb_intersect << " nb_intersect=" << nb_intersect << endl;
  cout << "   alpha=" << 100*temps_alpha/temps_total << "%" << endl;
  cout << "      triangulation=" << 100*temps_triangulation/temps_total << "%          t_moy=" << temps_triangulation/nb_convex_hull << " nb_triangulation=" << nb_convex_hull << endl;
  cout << "      convex_hull=" << 100*temps_convex_hull/temps_total << "%          t_moy=" << temps_convex_hull/nb_convex_hull << " nb_convex_hull=" << nb_convex_hull << endl;
  cout << "      volume=" << 100*temps_volume/temps_total << "%" << endl;
  cout << "      kappa 3d=" << 100*temps_kappa1/temps_total << "%          t_moy=" << temps_kappa1/nb_kappa1 << " nb_kappa1=" << nb_kappa1 << endl;
  cout << "      kappa 2d=" << 100*temps_kappa2/temps_total << "%          t_moy=" << temps_kappa2/nb_kappa2 << " nb_kappa2=" << nb_kappa2 << endl;
  cout << "triangularisation=" << 100*temps_triangularisation/temps_total << "%" << endl;
  cout << "   triangulation=" << 100*temps_triangulation2/temps_total << "%" << endl;
  cout << "Reste=" << 100-100*(temps_bbox+temps_vertices+temps_do_intersect+temps_triangularisation)/temps_total << "%" << endl;
  cout << "########################################" << endl;
	
}




