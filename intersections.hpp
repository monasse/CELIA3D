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

/*! \file
  \authors Maria Adela Puscas and Laurent Monasse
  \brief Inclusion of files and definition of types from the \b CGAL library.
  \warning  <b> Specific coupling file ! </b>
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
#include <CGAL/convex_hull_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/centroid.h>
#include <CGAL/number_utils.h>

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


typedef CGAL::Triangulation_2<Kernel>    Triangulation_2;


typedef Kernel::Point_2                  Point_2;

typedef CGAL::Triangle_2<Kernel>         Triangle_2;
typedef std::vector<Triangle_2>          Triangles_2;
typedef std::vector<Point_2>             Points_2;
typedef Kernel::Segment_2                Segment_2;
typedef CGAL::Bbox_2                     Bbox_2;
typedef Triangles_2::iterator            Triangle2_iterator;
typedef Triangles::iterator              Triangle3_iterator;


/*! \brief Triangulation of the faces of a fluid cubic cell.
   \details Bbox is a 3D bounding box. This allows to use the member functions of class <b> CGAL::Bbox_3  </b>. 
   \param cel Box_3d (fluid cubic cell)
   \param trianglesB list of the triangle splitting of the faces of Box
   \warning Specific coupling procedure !
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
/*! 
  \brief Edge/triangle intersection
  \return std:vector<Point_3>
 */
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


/*! 
  \brief New version of the triangle/triangle intersection
  \details Goal: accelerate the computation, since we do not need the precise intersection type, onl the list of extremal points of the intersection.\n
  Reduces to an edge/triangle intersection in CGAL
*/
std::vector<Point_3> intersection_bis(const Triangle_3& t1, const Triangle_3& t2)
{
  std::vector<Point_3> result;
  //Intersections between the segments of t1 and triangle t2
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
  //Intersections between the segments of t2 and triangle t1
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

/* ! \brief Determine whether a point is inside a tetrahedron
   \details Simple application of CGAL function <b> Tetrahedron_3::has_on_bounded_side <\b>
 */
bool inside_tetra(const Tetrahedron &tetra, const Point_3& P){
  return tetra.has_on_bounded_side(P);
}

/*! \brief Test whether points are coplanar
 */
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
