#include "solide.hpp"
#include "intersections.hpp"
#ifndef SOLIDE_CPP
#define SOLIDE_CPP

//const double eps_relat = numeric_limits<double>::epsilon();
const double eps_relat =0.000001;
Particule::Particule()
{   
	min_x = 0.; 
	min_y = 0.;
	min_z = 0.;
	max_x = 1. ;
	max_y = 1.;
	max_z = 1.;
	
	cube = true;
	
	const Point_3 s1(min_x,min_y,min_z);
	const Point_3 r1(max_x, min_y, min_z);
	const Point_3 t1(max_x, max_y, min_z);
	const Point_3 v1(min_x, max_y, min_z);
	
	
	const Point_3 s2(min_x,min_y,max_z);
	const Point_3 r2(max_x, min_y, max_z);
	const Point_3 t2(max_x, max_y, max_z);
	const Point_3 v2(min_x, max_y, max_z);
	
	
	centre[0] = centroid( s1, s2, v1,v2);
	centre[1] = centroid( r1, r2, t1,t2);
	
	centre[2] = centroid( s1, s2, v1,v2);
	centre[3] = centroid( v1, v2, t1,t2);
	
	centre[4] = centroid( s1, t1, r1,v1);
	centre[5] = centroid( s2, r2, t2,v2);
	
	Point_3 center((min_x + max_x)/2., (min_y + max_y)/2., (min_z + max_z)/2.);
	Vector_3 norm;
	
	//face1

	Triangle_3 Tri1(s1,r1,v1);
	norm= orthogonal_vector(s1,r1,v1);
	Vector_3 vect0(center,s1);
	if((CGAL::to_double(vect0*norm)) >=0.){ normales[0] = norm;}
	else {
		norm = orthogonal_vector(s1,v1,r1);
		normales[0] = norm;
	}
	
	Triangle_3 Tri2(t1,r1,v1);
	norm= orthogonal_vector(t1,r1,v1);
	Vector_3 vect1(center,t1);
	if((CGAL::to_double(vect1*norm)) >=0.){ normales[1] = norm;}
	else {
		norm = orthogonal_vector(t1,v1,r1);
		normales[1] = norm;
	}
	
	
	triangles.push_back(Tri1);
	triangles.push_back(Tri2);
	
	
	//face2
	
	Triangle_3 Tri5(s2,r2,v2);
	norm= orthogonal_vector(s2,r2,v2);
	Vector_3 vect2(center,s2);
	if((CGAL::to_double(vect2*norm)) >=0.){ normales[2] = norm;}
	else {
		norm = orthogonal_vector(s2,v2,r2);
		normales[2] = norm;
	}
	
	Triangle_3 Tri6(t2,r2,v2);
	norm= orthogonal_vector(t2,r2,v2);
	Vector_3 vect3(center,t2);
	if((CGAL::to_double(vect3*norm)) >=0.){ normales[3] = norm;}
	else {
		norm = orthogonal_vector(t2,v2,r2);
		normales[3] = norm;
	}
	
	triangles.push_back(Tri5);
	triangles.push_back(Tri6);
	
	
	//face3
	Triangle_3 Tri9(s2,s1,v2);
	norm= orthogonal_vector(s2,s1,v2);
	Vector_3 vect4(center,s2);
	if((CGAL::to_double(vect4*norm)) >=0.){ normales[4] = norm;}
	else {
		norm = orthogonal_vector(s2,v2,s1);
		normales[4] = norm;
	}
	Triangle_3 Tri10(v1,s1,v2);
	norm= orthogonal_vector(v1,s1,v2);
	Vector_3 vect5(center,v1);
	if((CGAL::to_double(vect5*norm)) >=0.){ normales[5] = norm;}
	else {
		norm = orthogonal_vector(v1,v2,s1);
		normales[5] = norm;
	}
	
	triangles.push_back(Tri9);
	triangles.push_back(Tri10);
	
	
	//face4	
	Triangle_3 Tri13(r2,r1,t2);
	norm= orthogonal_vector(r2,r1,t2);
	Vector_3 vect6(center,r2);
	if((CGAL::to_double(vect6*norm)) >=0.){ normales[6] = norm;}
	else {
		norm = orthogonal_vector(r2,t2,r1);
		normales[6] = norm;
	}
	Triangle_3 Tri14(t1,r1,t2);
	norm= orthogonal_vector(t1,r1,t2);
	Vector_3 vect7(center,t1);
	if((CGAL::to_double(vect7*norm)) >=0.){ normales[7] = norm;}
	else {
		norm = orthogonal_vector(t1,t2,r1);
		normales[7] = norm;
	}
	
	triangles.push_back(Tri13);
	triangles.push_back(Tri14);
	
	
	//face5	
	Triangle_3 Tri17(v2,v1,t2);
	norm= orthogonal_vector(v2,v1,t2);
	Vector_3 vect8(center,v2);
	if((CGAL::to_double(vect8*norm)) >=0.){ normales[8] = norm;}
	else {
		norm = orthogonal_vector(v2,t2,v1);
		normales[8] = norm;
	}
	Triangle_3 Tri18(t1,v1,t2);
	norm= orthogonal_vector(t1,v1,t2);
	Vector_3 vect9(center,t1);
	if((CGAL::to_double(vect9*norm)) >=0.){ normales[9] = norm;}
	else {
		norm = orthogonal_vector(t1,t2,v1);
		normales[9] = norm;
	}
	
	triangles.push_back(Tri17);
	triangles.push_back(Tri18);
	
	
	//face6
	
	Triangle_3 Tri21(s2,s1,r2);
	norm= orthogonal_vector(s2,s1,r2);
	Vector_3 vect10(center,s2);
	if((CGAL::to_double(vect10*norm)) >=0.){ normales[10] = norm;}
	else {
		norm = orthogonal_vector(s2,r2,s1);
		normales[10] = norm;
	}
	Triangle_3 Tri22(r1,s1,r2);
	norm= orthogonal_vector(r1,s1,r2);
	Vector_3 vect11(center,r1);
	if((CGAL::to_double(vect11*norm)) >=0.){ normales[11] = norm;}
	else {
		norm = orthogonal_vector(r1,r2,s1);
		normales[11] = norm;
	}
	
	triangles.push_back(Tri21);
	triangles.push_back(Tri22); 
	
}

Particule::Particule(const double x_min, const double y_min, const double z_min, 
							 const double x_max, const double y_max,const double z_max)
{   
	min_x = x_min; 
	min_y = y_min;
	min_z = z_min;
	max_x = x_max ;
	max_y = y_max;
	max_z = z_max;
	
	cube = true;
	
	const Point_3 s1(min_x,min_y,min_z);
	const Point_3 r1(max_x, min_y, min_z);
	const Point_3 t1(max_x, max_y, min_z);
	const Point_3 v1(min_x, max_y, min_z);
	
	
	const Point_3 s2(min_x,min_y,max_z);
	const Point_3 r2(max_x, min_y, max_z);
	const Point_3 t2(max_x, max_y, max_z);
	const Point_3 v2(min_x, max_y, max_z);
	
	centre[0] = centroid( s1, s2, v1,v2);
	centre[1] = centroid( r1, r2, t1,t2);
	
	centre[2] = centroid( s1, s2, v1,v2);
	centre[3] = centroid( v1, v2, t1,t2);
	
	centre[4] = centroid( s1, t1, r1,v1);
	centre[5] = centroid( s2, r2, t2,v2);
	
	
	Point_3 center((min_x + max_x)/2., (min_y + max_y)/2., (min_z + max_z)/2.);
	Vector_3 norm;
	
	//face1
	
	Triangle_3 Tri1(s1,r1,v1);
	norm= orthogonal_vector(s1,r1,v1);
	Vector_3 vect0(center,s1);
	if((CGAL::to_double(vect0*norm)) >=0.){ normales[0] = norm;}
	else {
		norm = orthogonal_vector(s1,v1,r1);
		normales[0] = norm;
	}
	
	Triangle_3 Tri2(t1,r1,v1);
	norm= orthogonal_vector(t1,r1,v1);
	Vector_3 vect1(center,t1);
	if((CGAL::to_double(vect1*norm)) >=0.){ normales[1] = norm;}
	else {
		norm = orthogonal_vector(t1,v1,r1);
		normales[1] = norm;
	}
	
	
	triangles.push_back(Tri1);
	triangles.push_back(Tri2);
	
	
	//face2
	
	Triangle_3 Tri5(s2,r2,v2);
	norm= orthogonal_vector(s2,r2,v2);
	Vector_3 vect2(center,s2);
	if((CGAL::to_double(vect2*norm)) >=0.){ normales[2] = norm;}
	else {
		norm = orthogonal_vector(s2,v2,r2);
		normales[2] = norm;
	}
	
	Triangle_3 Tri6(t2,r2,v2);
	norm= orthogonal_vector(t2,r2,v2);
	Vector_3 vect3(center,t2);
	if((CGAL::to_double(vect3*norm)) >=0.){ normales[3] = norm;}
	else {
		norm = orthogonal_vector(t2,v2,r2);
		normales[3] = norm;
	}
	
	triangles.push_back(Tri5);
	triangles.push_back(Tri6);
	
	
	//face3
	Triangle_3 Tri9(s2,s1,v2);
	norm= orthogonal_vector(s2,s1,v2);
	Vector_3 vect4(center,s2);
	if((CGAL::to_double(vect4*norm)) >=0.){ normales[4] = norm;}
	else {
		norm = orthogonal_vector(s2,v2,s1);
		normales[4] = norm;
	}
	Triangle_3 Tri10(v1,s1,v2);
	norm= orthogonal_vector(v1,s1,v2);
	Vector_3 vect5(center,v1);
	if((CGAL::to_double(vect5*norm)) >=0.){ normales[5] = norm;}
	else {
		norm = orthogonal_vector(v1,v2,s1);
		normales[5] = norm;
	}
	
	triangles.push_back(Tri9);
	triangles.push_back(Tri10);
	
	
	//face4	
	Triangle_3 Tri13(r2,r1,t2);
	norm= orthogonal_vector(r2,r1,t2);
	Vector_3 vect6(center,r2);
	if((CGAL::to_double(vect6*norm)) >=0.){ normales[6] = norm;}
	else {
		norm = orthogonal_vector(r2,t2,r1);
		normales[6] = norm;
	}
	Triangle_3 Tri14(t1,r1,t2);
	norm= orthogonal_vector(t1,r1,t2);
	Vector_3 vect7(center,t1);
	if((CGAL::to_double(vect7*norm)) >=0.){ normales[7] = norm;}
	else {
		norm = orthogonal_vector(t1,t2,r1);
		normales[7] = norm;
	}
	
	triangles.push_back(Tri13);
	triangles.push_back(Tri14);
	
	
	//face5	
	Triangle_3 Tri17(v2,v1,t2);
	norm= orthogonal_vector(v2,v1,t2);
	Vector_3 vect8(center,v2);
	if((CGAL::to_double(vect8*norm)) >=0.){ normales[8] = norm;}
	else {
		norm = orthogonal_vector(v2,t2,v1);
		normales[8] = norm;
	}
	Triangle_3 Tri18(t1,v1,t2);
	norm= orthogonal_vector(t1,v1,t2);
	Vector_3 vect9(center,t1);
	if((CGAL::to_double(vect9*norm)) >=0.){ normales[9] = norm;}
	else {
		norm = orthogonal_vector(t1,t2,v1);
		normales[9] = norm;
	}
	
	triangles.push_back(Tri17);
	triangles.push_back(Tri18);
	
	
	//face6
	
	Triangle_3 Tri21(s2,s1,r2);
	norm= orthogonal_vector(s2,s1,r2);
	Vector_3 vect10(center,s2);
	if((CGAL::to_double(vect10*norm)) >=0.){ normales[10] = norm;}
	else {
		norm = orthogonal_vector(s2,r2,s1);
		normales[10] = norm;
	}
	Triangle_3 Tri22(r1,s1,r2);
	norm= orthogonal_vector(r1,s1,r2);
	Vector_3 vect11(center,r1);
	if((CGAL::to_double(vect11*norm)) >=0.){ normales[11] = norm;}
	else {
		norm = orthogonal_vector(r1,r2,s1);
		normales[11] = norm;
	}
	
	triangles.push_back(Tri21);
	triangles.push_back(Tri22); 


}


Particule::Particule(const double x_min, const double y_min, const double z_min, 
							 const double x_max, const double y_max,const double z_max, 
							 const Point_3 s1, const Point_3 r1, const Point_3 t1,const Point_3 v1,
							 const Point_3 s2, const Point_3 r2, const Point_3 t2, const Point_3 v2)
{   
	min_x = x_min; 
	min_y = y_min;
	min_z = z_min;
	max_x = x_max ;
	max_y = y_max;
	max_z = z_max;
	
	cube = false;
	
	centre[0] = centroid( s1, s2, v1,v2);
	centre[1] = centroid( r1, r2, t1,t2);
	
	centre[2] = centroid( s1, s2, v1,v2);
	centre[3] = centroid( v1, v2, t1,t2);
	
	centre[4] = centroid( s1, t1, r1,v1);
	centre[5] = centroid( s2, r2, t2,v2);
	Point_3 center((min_x + max_x)/2., (min_y + max_y)/2., (min_z + max_z)/2.);
	Vector_3 norm;
	
	//face1
	
	Triangle_3 Tri1(s1,r1,v1);
	norm= orthogonal_vector(s1,r1,v1);
	Vector_3 vect0(center,s1);
	if((CGAL::to_double(vect0*norm)) >=0.){ normales[0] = norm;}
	else {
		norm = orthogonal_vector(s1,v1,r1);
		normales[0] = norm;
	}
	//std::cout<<CGAL::to_double(vect0*norm)<<std::endl;
	
	Triangle_3 Tri2(t1,r1,v1);
	norm= orthogonal_vector(t1,r1,v1);
	Vector_3 vect1(center,t1);
	if((CGAL::to_double(vect1*norm)) >=0.){ normales[1] = norm;}
	else {
		norm = orthogonal_vector(t1,v1,r1);
		normales[1] = norm;
	}
	//std::cout<<CGAL::to_double(vect1*norm)<<std::endl;
	
	
	triangles.push_back(Tri1);
	triangles.push_back(Tri2);
	
	
	//face2
	
	Triangle_3 Tri5(s2,r2,v2);
	norm= orthogonal_vector(s2,r2,v2);
	Vector_3 vect2(center,s2);
	if((CGAL::to_double(vect2*norm)) >=0.){ normales[2] = norm;}
	else {
		norm = orthogonal_vector(s2,v2,r2);
		normales[2] = norm;
	}
	//std::cout<<CGAL::to_double(vect2*norm)<<std::endl;
	
	Triangle_3 Tri6(t2,r2,v2);
	norm= orthogonal_vector(t2,r2,v2);
	Vector_3 vect3(center,t2);
	if((CGAL::to_double(vect3*norm)) >=0.){ normales[3] = norm;}
	else {
		norm = orthogonal_vector(t2,v2,r2);
		normales[3] = norm;
	}
	//std::cout<<CGAL::to_double(vect3*norm)<<std::endl;
	
	triangles.push_back(Tri5);
	triangles.push_back(Tri6);
	
	
	//face3
	Triangle_3 Tri9(s2,s1,v2);
	norm= orthogonal_vector(s2,s1,v2);
	Vector_3 vect4(center,s2);
	if((CGAL::to_double(vect4*norm)) >=0.){ normales[4] = norm;}
	else {
		norm = orthogonal_vector(s2,v2,s1);
		normales[4] = norm;
	}
	//std::cout<<CGAL::to_double(vect4*norm)<<std::endl;
	
	Triangle_3 Tri10(v1,s1,v2);
	norm= orthogonal_vector(v1,s1,v2);
	Vector_3 vect5(center,v1);
	if((CGAL::to_double(vect5*norm)) >=0.){ normales[5] = norm;}
	else {
		norm = orthogonal_vector(v1,v2,s1);
		normales[5] = norm;
	}
	//std::cout<<CGAL::to_double(vect5*norm)<<std::endl;
	
	triangles.push_back(Tri9);
	triangles.push_back(Tri10);
	
	
	//face4	
	Triangle_3 Tri13(r2,r1,t2);
	norm= orthogonal_vector(r2,r1,t2);
	Vector_3 vect6(center,r2);
	if((CGAL::to_double(vect6*norm)) >=0.){ normales[6] = norm;}
	else {
		norm = orthogonal_vector(r2,t2,r1);
		normales[6] = norm;
	}
	//std::cout<<CGAL::to_double(vect6*norm)<<std::endl;
	
	Triangle_3 Tri14(t1,r1,t2);
	norm= orthogonal_vector(t1,r1,t2);
	Vector_3 vect7(center,t1);
	if((CGAL::to_double(vect7*norm)) >=0.){ normales[7] = norm;}
	else {
		norm = orthogonal_vector(t1,t2,r1);
		normales[7] = norm;
	}
	//std::cout<<CGAL::to_double(vect7*norm)<<std::endl;
	
	triangles.push_back(Tri13);
	triangles.push_back(Tri14);
	
	
	//face5	
	Triangle_3 Tri17(v2,v1,t2);
	norm= orthogonal_vector(v2,v1,t2);
	Vector_3 vect8(center,v2);
	if((CGAL::to_double(vect8*norm)) >=0.){ normales[8] = norm;}
	else {
		norm = orthogonal_vector(v2,t2,v1);
		normales[8] = norm;
	}
	//std::cout<<CGAL::to_double(vect8*norm)<<std::endl;
	
	Triangle_3 Tri18(t1,v1,t2);
	norm= orthogonal_vector(t1,v1,t2);
	Vector_3 vect9(center,t1);
	if((CGAL::to_double(vect9*norm)) >=0.){ normales[9] = norm;}
	else {
		norm = orthogonal_vector(t1,t2,v1);
		normales[9] = norm;
	}
	//std::cout<<CGAL::to_double(vect9*norm)<<std::endl;
	
	triangles.push_back(Tri17);
	triangles.push_back(Tri18);
	
	
	//face6
	
	Triangle_3 Tri21(s2,s1,r2);
	norm= orthogonal_vector(s2,s1,r2);
	Vector_3 vect10(center,s2);
	if((CGAL::to_double(vect10*norm)) >=0.){ normales[10] = norm;}
	else {
		norm = orthogonal_vector(s2,r2,s1);
		normales[10] = norm;
	}
	//std::cout<<CGAL::to_double(vect10*norm)<<std::endl;
	
	Triangle_3 Tri22(r1,s1,r2);
	norm= orthogonal_vector(r1,s1,r2);
	Vector_3 vect11(center,r1);
	if((CGAL::to_double(vect11*norm)) >=0.){ normales[11] = norm;}
	else {
		norm = orthogonal_vector(r1,r2,s1);
		normales[11] = norm;
	}
	//std::cout<<CGAL::to_double(vect11*norm)<<std::endl;
	triangles.push_back(Tri21);
	triangles.push_back(Tri22); 

}
//Destructeur
Particule::~Particule(){
	
}



void Particule::Affiche(){
	
//	std::cout<<" volume of solide := "<<volume()<<std::endl;
// 	std::cout<<" Point min := "<< min_x<<std::endl;
// 	std::cout<<" Point max := "<< max_x<<std::endl;
}




double Particule::volume(){
	
	double vol = 0.;
	std::vector<Point_3> Points_poly; 
	
	for(int l= 0; l<triangles.size(); l++)
	{
		Points_poly.push_back(triangles[l].operator[](0));
		Points_poly.push_back(triangles[l].operator[](1));
		Points_poly.push_back(triangles[l].operator[](2));
	}
	
	Finite_cells_iterator cit;
	Triangulation T(Points_poly.begin(), Points_poly.end());
	
	for (cit = T.finite_cells_begin(); cit!= T.finite_cells_end(); cit++){
		vol+= CGAL::to_double(T.tetrahedron( cit).volume());
	}
	
	return vol;
}


void Solide::impression(int n){ //Sortie au format vtk
int nb_part = solide.size();

int nb_triangles = 0.;
for(int it=0; it<nb_part; it++){
	nb_triangles += solide[it].triangles.size();
}

const char* solidevtk;
{
	std::ostringstream oss;
	oss << "resultats/solide" << n << ".vtk";
	string s = oss.str();
	//cout << s << endl;
	solidevtk = s.c_str();
}

//Ouverture des flux en donne en ecriture
std::ofstream vtk(solidevtk,ios::out);
if(vtk)
{
	// cout <<"ouverture de xt.vtk reussie" << endl;
} else {
	cout <<"ouverture de solide" << n << ".vtk rate" << endl;
}
//Initialisation du fichier vtk
vtk << "# vtk DataFile Version 3.0" << endl;
vtk << "#Simulation Euler" << endl;
vtk << "ASCII" << endl;
vtk<<"\n";
vtk << "DATASET UNSTRUCTURED_GRID" << endl;
vtk << "POINTS " << 3*nb_triangles << " DOUBLE" << endl;

for(int it=0; it<nb_part; it++){
	for(int l= 0; l<solide[it].triangles.size(); l++){
		vtk << solide[it].triangles[l].operator[](0).operator[](0) << " " << solide[it].triangles[l].operator[](0).operator[](1) << " " << solide[it].triangles[l].operator[](0).operator[](2) << endl;
		vtk << solide[it].triangles[l].operator[](1).operator[](0) << " " << solide[it].triangles[l].operator[](1).operator[](1) << " " << solide[it].triangles[l].operator[](1).operator[](2) << endl;
		vtk << solide[it].triangles[l].operator[](2).operator[](0) << " " << solide[it].triangles[l].operator[](2).operator[](1) << " " << solide[it].triangles[l].operator[](2).operator[](2) << endl;
	}
}
vtk<<"\n";
vtk << "CELLS " << nb_triangles << " " << 4*nb_triangles<< endl;
int num=0;
for(int it=0; it<nb_part; it++){
	for(int l= 0; l<solide[it].triangles.size(); l++){
		vtk << 3 << " " << 3*num << " " << 3*num+1 << " " << 3*num+2 << endl;
		num++;
	}
}
vtk << "\n";
vtk << "CELL_TYPES " << nb_triangles << endl;
for(int l= 0; l<nb_triangles; l++)
{
	vtk << 5 << endl;
}
vtk << "\n";
vtk << "CELL_DATA " << nb_triangles << endl;
//Deplacement
vtk << "SCALARS pression double 1" << endl;
vtk << "LOOKUP_TABLE default" << endl;
for(int it=0; it<nb_part; it++){
	for(int l= 0; l<solide[it].triangles.size(); l++)
	{
		vtk << 0. << endl;
	}
}
}
	
bool inside_convex_polygon(const Particule& S, const Point_3& P){
	
	bool in = false;
	
	if((S.min_x - P.x())<= eps_relat && (S.min_y - P.y())<= eps_relat && 
		 (S.min_z - P.z())<= eps_relat && (S.max_x - P.x())>=-eps_relat && 
		 (S.max_y - P.y())>=-eps_relat && (S.max_z - P.z())>=-eps_relat )
	{
		if(S.cube) {in = true;}
		else{
			in = true;
			for(int l= 0; l<S.triangles.size(); l++){
				Point_3 vertex = S.triangles[l].operator[](0);
				Vector_3 vect(P,vertex);
				if((CGAL::to_double(vect*S.normales[l])) < 0.){in = false;}
			}
		}
	}
	
	return in;
}	

bool box_inside_convex_polygon(const Particule& S, const Bbox& cell){
	
	bool in = false;
	
	if ((S.min_x - cell.xmin()) <= eps_relat && (S.min_y - cell.ymin() <= eps_relat) && 
		  (S.min_z - cell.zmin()) <= eps_relat && (S.max_x - cell.xmax() >=-eps_relat) && 
		  (S.max_y - cell.ymax()) >=-eps_relat && (S.max_z - cell.zmax()>=-eps_relat) ) 
	{
		if(S.cube) { return S.cube;}
		
		else{
			
	  in = true;
		
		Point_3 p1(cell.xmin(),cell.ymin(),cell.zmin());
		in = inside_convex_polygon(S,p1);
		if(!in) {return in;}
		
		Point_3 p2(cell.xmax(),cell.ymax(),cell.zmax());
		in = inside_convex_polygon(S,p2);
		if(!in) {return in;}
		
		Point_3 p3(cell.xmax(),cell.ymin(),cell.zmin());
		in = inside_convex_polygon(S,p3);
		if(!in) {return in;}
		
		Point_3 p4(cell.xmax(),cell.ymin(),cell.zmax());
		in = inside_convex_polygon(S,p4);
		if(!in) {return in;}
		
		Point_3 p5(cell.xmax(),cell.ymax(),cell.zmin());
		in = inside_convex_polygon(S,p5);
		if(!in) {return in;}
		
		Point_3 p6(cell.xmin(),cell.ymax(),cell.zmin());
		in = inside_convex_polygon(S,p6);
		if(!in) {return in;}
		
		Point_3 p7(cell.xmin(),cell.ymax(),cell.zmax());
		in = inside_convex_polygon(S,p7);
		if(!in) {return in;}
		
		Point_3 p8(cell.xmin(),cell.ymin(),cell.zmax());
		in = inside_convex_polygon(S,p8);
		if(!in) {return in;}
		
		}
	}
	
	return in;
}


bool inside_box(const Bbox& cell, const Point_3& P){
	
	bool in = false;
	
	if((cell.xmin() - P.x())<= eps_relat && (cell.ymin() - P.y())<= eps_relat &&
		 (cell.zmin() - P.z())<= eps_relat && (cell.xmax() - P.x())>=-eps_relat &&
		 (cell.ymax() - P.y())>=-eps_relat && (cell.zmax() - P.z())>=-eps_relat )
	  { in = true; }
	
	return in;
}


Solide::Solide(std::vector<Particule> & Part){
	
	for(int i=0; i<Part.size(); i++){
		solide.push_back(Part[i]);
	}
	
}

Solide::~Solide(){   
}


void Solide::Affiche(){
	
	for(int i=0; i<solide.size(); i++){
		solide[i].Affiche();
	}
}

#endif
