#ifndef INTERSECTIONS
#define INTERSECTIONS

#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <cassert>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Timer.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/centroid.h>
#include <CGAL/number_utils.h>
//Constrained Triangulation
#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

//

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_3                  Point_3;
typedef Kernel::Vector_3                 Vector_3;
typedef Kernel::Line_3                   Line_3;
typedef Kernel::Plane_3                  Plane_3;
typedef CGAL::Triangle_3<Kernel>         Triangle_3;
typedef CGAL::Plane_3<Kernel>            Plane_3;
typedef std::vector<Triangle_3>          Triangles;
typedef CGAL::Tetrahedron_3<Kernel>     Tetrahedron;
typedef std::vector<Point_3>             Points;
typedef Kernel::Segment_3                Segment_3;
typedef CGAL::Bbox_3                     Bbox;
typedef CGAL::Polyhedron_3<Kernel>       Polyhedron_3;
typedef CGAL::Triangulation_3<Kernel>    Triangulation;
typedef CGAL::Aff_transformation_3<Kernel>  Aff_transformation_3;

typedef Triangulation::Finite_facets_iterator Finite_faces_iterator;
typedef Triangulation::Finite_cells_iterator Finite_cells_iterator;
typedef Polyhedron_3::Facet_iterator   Facet_iterator;
typedef Triangles::iterator            Triangle3_iterator;


//Triangulation 2d 
typedef CGAL::Triangulation_2<Kernel>    Triangulation_2;

////Constrained Triangulation
typedef CGAL::Triangulation_vertex_base_2<Kernel>                     Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel>           Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                   TDS;
typedef CGAL::Exact_predicates_tag                                    Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag> CDT;
typedef CDT::Vertex  	                                             Vertex_cdt;  	 
typedef CDT::Face 	                                               Face_cdt; 

typedef CDT::Face_handle 	                                           Face_cdt_handle; 
typedef CDT::Vertex_handle  	                                       Vertex_cdt_handle; 

typedef Kernel::Point_2                  Point_2;
typedef CGAL::Triangle_2<Kernel>         Triangle_2;
typedef std::vector<Triangle_2>          Triangles_2;
typedef std::vector<Point_2>             Points_2;
typedef Kernel::Segment_2                Segment_2;
typedef CGAL::Bbox_2                     Bbox_2;
typedef Triangles_2::iterator            Triangle2_iterator;
typedef Triangles::iterator            Triangle3_iterator;

void triang_cellule(const Bbox& cel, Triangles& trianglesB){
  
           //face1
	    Point_3 s1B(cel.xmin(), cel.ymin(), cel.zmin());
	    Point_3 r1B(cel.xmax(), cel.ymin(), cel.zmin());
	    Point_3 t1B(cel.xmax(), cel.ymax(), cel.zmin());
	    Point_3 v1B(cel.xmin(), cel.ymax(), cel.zmin());
	    
	    Triangle_3 Tri1B(s1B,r1B,v1B);
	    Triangle_3 Tri2B(t1B,r1B,v1B);

	    
	    trianglesB.push_back(Tri1B);
	    trianglesB.push_back(Tri2B);

	    
	    //face2
	    Point_3 s2B(cel.xmin(), cel.ymin(), cel.zmax());
	    Point_3 r2B(cel.xmax(), cel.ymin(), cel.zmax());
	    Point_3 t2B(cel.xmax(), cel.ymax(), cel.zmax());
	    Point_3 v2B(cel.xmin(), cel.ymax(), cel.zmax());
	    
	    Triangle_3 Tri5B(s2B,r2B,v2B);
	    Triangle_3 Tri6B(t2B,r2B,v2B);
	 
	    
	    trianglesB.push_back(Tri5B);
	    trianglesB.push_back(Tri6B);

	    
	    
	    //face3
	    
	    
	    Triangle_3 Tri9B(s2B,s1B,v2B);
	    Triangle_3 Tri10B(v1B,s1B,v2B);

	    
	    trianglesB.push_back(Tri9B);
	    trianglesB.push_back(Tri10B);
;
	    
	    //face4
	    
	    
	    Triangle_3 Tri13B(r2B,r1B,t2B);
	    Triangle_3 Tri14B(t1B,r1B,t2B);
	
	    trianglesB.push_back(Tri13B);
	    trianglesB.push_back(Tri14B);

	    
	    //face5
	    
	    
	    Triangle_3 Tri17B(v2B,v1B,t2B);
	    Triangle_3 Tri18B(t1B,v1B,t2B);

	    
	    trianglesB.push_back(Tri17B);
	    trianglesB.push_back(Tri18B);

	    
	    //face6
	    
	    
	    Triangle_3 Tri21B(s2B,s1B,r2B);
	    Triangle_3 Tri22B(r1B,s1B,r2B);

	    
	    trianglesB.push_back(Tri21B);
	    trianglesB.push_back(Tri22B);

  
}

#endif
