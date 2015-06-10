/*!
 *  \file intersections.hpp
 *  \brief Inclusion des fichiers et d&eacute;finitions des types de la libraire \b CGAL.
    \details Les inclusions des fichier et les d&eacute;finitions des types de la libraire CGAL se font uniquement dans ce fichier. Dans ce fichier on d&eacute;finit &eacute;galement une fonction de triangulation des faces d'une cellule cubique fluide.
    \warning  <b> Fichier sp&eacute;cifique au couplage! </b>
 */

#ifndef INTERSECTIONS
#define INTERSECTIONS

#include <iostream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <cassert>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
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

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel   IK;
typedef CGAL::Cartesian_converter<Kernel,IK> Exact_to_Inexact;
typedef CGAL::Cartesian_converter<IK,Kernel> Inexact_to_Exact;
typedef Kernel::Point_3                  Point_3;
typedef IK::Point_3                      InexactPoint_3;
typedef Kernel::Vector_3                 Vector_3;
typedef Kernel::Line_3                   Line_3;
typedef Kernel::Plane_3                  Plane_3;
typedef CGAL::Triangle_3<Kernel>         Triangle_3;
typedef CGAL::Triangle_3<IK>             InexactTriangle_3;
typedef CGAL::Plane_3<Kernel>            Plane_3;
typedef std::vector<Triangle_3>          Triangles;
typedef CGAL::Tetrahedron_3<Kernel>      Tetrahedron;
typedef CGAL::Tetrahedron_3<IK>          InexactTetrahedron;
typedef std::vector<Point_3>             Points;
typedef Kernel::Segment_3                Segment_3;
typedef CGAL::Bbox_3                     Bbox;
typedef CGAL::Polyhedron_3<IK>           InexactPolyhedron_3;
typedef CGAL::Polyhedron_3<Kernel>       Polyhedron_3;
typedef CGAL::Triangulation_3<Kernel>    ExactTriangulation;
typedef CGAL::Triangulation_3<IK>        InexactTriangulation;
typedef CGAL::Aff_transformation_3<Kernel>  Aff_transformation_3;

typedef ExactTriangulation::Finite_facets_iterator ExactFinite_faces_iterator;
typedef ExactTriangulation::Finite_cells_iterator ExactFinite_cells_iterator;
typedef InexactTriangulation::Finite_facets_iterator InexactFinite_faces_iterator;
typedef InexactTriangulation::Finite_cells_iterator InexactFinite_cells_iterator;
typedef Polyhedron_3::Facet   Facet;
typedef Polyhedron_3::Facet_iterator   Facet_iterator;
typedef Polyhedron_3::Vertex_iterator   Vertex_iterator;
typedef Polyhedron_3::Plane_iterator   Plane_iterator;
typedef Polyhedron_3::Halfedge_around_vertex_circulator   Halfedge_around_vertex_circulator;
typedef InexactPolyhedron_3::Facet   InexactFacet;
typedef InexactPolyhedron_3::Facet_iterator   InexactFacet_iterator;
typedef InexactPolyhedron_3::Vertex_iterator   InexactVertex_iterator;
typedef InexactPolyhedron_3::Plane_iterator   InexactPlane_iterator;
typedef InexactPolyhedron_3::Halfedge_around_vertex_circulator   InexactHalfedge_around_vertex_circulator;
typedef Triangles::iterator            Triangle3_iterator;


//Triangulation 2d 
typedef CGAL::Triangulation_2<Kernel>    Triangulation_2;

////Constrained Triangulation
typedef CGAL::Triangulation_vertex_base_2<Kernel>                      Vb;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel>            Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                    TDS;
typedef CGAL::Exact_predicates_tag                                     Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag>  CDT;
typedef CDT::Vertex  	                                           Vertex_cdt;  	 
typedef CDT::Face 	                                               Face_cdt; 

typedef CDT::Face_handle 	             Face_cdt_handle; 
typedef CDT::Vertex_handle  	         Vertex_cdt_handle; 

typedef Kernel::Point_2                  Point_2;

typedef CGAL::Triangle_2<Kernel>         Triangle_2;
typedef std::vector<Triangle_2>          Triangles_2;
typedef std::vector<Point_2>             Points_2;
typedef Kernel::Segment_2                Segment_2;
typedef CGAL::Bbox_2                     Bbox_2;
typedef Triangles_2::iterator            Triangle2_iterator;
typedef Triangles::iterator              Triangle3_iterator;


/**
 \fn void triang_cellule(const Bbox& cel, Triangles& trianglesB)
 \brief Triangulation des faces d'une cellule cubique fluide.
 \details Bbox est une Box 3d; une Box repr&eacute;sente une bo&icirc;te rectangulaire. Cette fa&ccedil;on de voir les cellules fluide permet de faire appel aux fonctions membres de la classe <b> CGAL::Bbox_3  </b>. 
 \param cel Box_3d (cellule cubique fluide)
 \param trianglesB liste des triangles d&eacute;crivant les faces de Box
 \warning Proc&eacute;dure sp&eacute;cifique au couplage!
 \return void
 */
void triang_cellule(const Bbox& cel, Triangles& trianglesB){
  
  trianglesB.clear();
  
	    Point_3 s1B(cel.xmin(), cel.ymin(), cel.zmin());
	    Point_3 r1B(cel.xmax(), cel.ymin(), cel.zmin());
	    Point_3 t1B(cel.xmax(), cel.ymax(), cel.zmin());
	    Point_3 v1B(cel.xmin(), cel.ymax(), cel.zmin());
     
        Point_3 s2B(cel.xmin(), cel.ymin(), cel.zmax());
        Point_3 r2B(cel.xmax(), cel.ymin(), cel.zmax());
        Point_3 t2B(cel.xmax(), cel.ymax(), cel.zmax());
        Point_3 v2B(cel.xmin(), cel.ymax(), cel.zmax());
    
	    //face1
	    Triangle_3 Tri1B(s1B,r1B,v1B);
	    Triangle_3 Tri2B(t1B,r1B,v1B);
	    trianglesB.push_back(Tri1B);
	    trianglesB.push_back(Tri2B);

	    
	    //face2
	    Triangle_3 Tri5B(s2B,r2B,v2B);
	    Triangle_3 Tri6B(t2B,r2B,v2B);
        trianglesB.push_back(Tri5B);
	    trianglesB.push_back(Tri6B);

	    //face3
	    Triangle_3 Tri9B(s2B,s1B,v2B);
	    Triangle_3 Tri10B(v1B,s1B,v2B);
	    trianglesB.push_back(Tri9B);
	    trianglesB.push_back(Tri10B);

	    
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

//Intersection arete/triangle
std::vector<Point_3> intersection_bis(const Segment_3& seg, const Triangle_3& t)
{
  std::vector<Point_3> result;
  Point_3 P;
  Segment_3 s;
  const CGAL::Object& intersec = CGAL::intersection(seg,t);
  if(CGAL::assign(P,intersec)){
    result.push_back(P);
  }
  else if(CGAL::assign(s,intersec)){
    result.push_back(s.operator[](0));
    result.push_back(s.operator[](1));
  }
    
  return result;
}

  
//Nouvelle version de l'intersection triangle/triangle
//But : accelerer le calcul, car nous n'avons pas besoin du type precis de l'intersection, juste de la liste des points extremes de l'intersection
//Reduit a l'intersection segment/triangle de CGAL
std::vector<Point_3> intersection_bis(const Triangle_3& t1, const Triangle_3& t2)
{
  std::vector<Point_3> result;
  //Intersections entre les segments de t1 et le triangle t2
  for(int k=0;k<3;k++){
    int kp = (k+1)%3;
    Point_3 P;
    Segment_3 seg;
    Segment_3 arete(t1.operator[](k),t1.operator[](kp));
    const CGAL::Object& intersec = CGAL::intersection(arete,t2);
    if(CGAL::assign(P,intersec)){
      result.push_back(P);
    }
    else if(CGAL::assign(seg,intersec)){
      result.push_back(seg.operator[](0));
      result.push_back(seg.operator[](1));
    }
  }
  //Intersections entre les segments de t2 et le triangle t1
  for(int k=0;k<3;k++){
    int kp = (k+1)%3;
    Point_3 P;
    Segment_3 seg;
    Segment_3 arete(t2.operator[](k),t2.operator[](kp));
    const CGAL::Object& intersec = CGAL::intersection(arete,t1);
    if(CGAL::assign(P,intersec)){
      result.push_back(P);
    }
    else if(CGAL::assign(seg,intersec)){
      result.push_back(seg.operator[](0));
      result.push_back(seg.operator[](1));
    }
  }
    
  return result;
}

bool inside_tetra(const Tetrahedron &tetra, const Point_3& P){
	
  //bool in = false;
  //in=tetra.has_on_negative_side(P);
	
  //return in;
  /*if(tetra.volume()>0.){
    return tetra.has_on_positive_side(P);
  } else {
    return tetra.has_on_negative_side(P);
    }*/
  return tetra.has_on_bounded_side(P);
}

bool coplanar(std::vector<Point_3>::iterator begin, std::vector<Point_3>::iterator end)
{
  bool test=true;

  bool test_confondus=true;
  Point_3 P1,P2;
  P1 = *begin;
  std::vector<Point_3>::iterator it;
  for(it=begin;it!=end && test_confondus;it++){
    P2 = *it;
    test_confondus = (P1==P2);
  }
  if(test_confondus){
    return true;
  } else {
    bool test_collinear=true;
    Point_3 P3;
    std::vector<Point_3>::iterator iter;
    for(iter=it;iter!=end && test_collinear;iter++){
      P3 = *iter;
      test_collinear = CGAL::collinear(P1,P2,P3);
    }
    if(test_collinear){
      return true;
    } else {
      Point_3 P4;
      for(std::vector<Point_3>::iterator iterat=iter;iterat!=end && test;iterat++){
	P4 = *iterat;
	test = CGAL::coplanar(P1,P2,P3,P4);
      }
    }
  }
  
  return test;
}


#endif
