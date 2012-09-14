#include <iostream>
#include <ctime>

//#include "fluide.hpp"
#include "fluide.cpp"
#include "couplage.cpp"
//#include "parametres.hpp"
#include "parametres.cpp"
#include "solide.cpp"


using namespace std;          // espace de nom standard


int main()
{		cout << "test" << endl;
	char tempsF[]="resultats/tempsF.dat";
	
	//Ouverture des flux en donne en ecriture
	std::ofstream out(tempsF,ios::out);
	if(out)
	{
		// cout <<"ouverture de xt.vtk reussie" << endl;
	} else {
		cout <<"ouverture de .dat rate" << endl;
	}
	
	Grille Fluide;
	
	
	double t=0.;


	
	std::ifstream maillage("maillage.dat",ios::in);
	if(maillage)
	{
		// cout <<"ouverture de xt.vtk reussie" << endl;
	} else {
		cout <<"ouverture de maillage.dat ratee" << endl;
	}

	
	//Recuperation du maillage solide
	int Npoint;
	string sp;
	maillage >> sp >> Npoint;
	const int nb_points = Npoint;

	vector<Point_3> Points(nb_points);

	for(int i=0;i<nb_points;i++){
	  double x,y,z;
	  maillage >> x >> y >> z;
	  Points[i] = Point_3(x,y,z);
	}

	int Npart;
	string sP;
	maillage >> sP >> Npart;
	const int nb_particule = Npart;
	
	vector<Particule> P(nb_particule);

	bool points_particules[nb_points][nb_particule];

	for(int i=0;i<nb_particule;i++){
	  for(int j=0;j<nb_points;j++){
		//Remise a zero du tableau
		points_particules[j][i] = false;
	  }
	}
	for(int i=0;i<nb_particule;i++){
	  int Nfaces;
	  bool fixe;
	  double X,Y,Z,u,v,w,theta,phi,psi,xmin,ymin,zmin,xmax,ymax,zmax;
	  string s;
	  maillage >> s >> Nfaces >> fixe;
	  maillage >> s >> X >> Y >> Z;
	  maillage >> s >> u >> v >> w;
	  maillage >> s >> theta >> phi >> psi;
	  xmin = X;
	  ymin = Y;
	  zmin = Z;
	  const int nb_faces = Nfaces;
	  std::vector<Face> Faces(nb_faces);
	  for(int j=0;j<nb_faces;j++){
		int Nvertex;
		maillage >> Nvertex;
		const int nb_vertex = Nvertex;
		std::vector<Vertex> Vertex(nb_vertex);
		for(int k=0;k<nb_vertex;k++){
		  int p;
		  maillage >> p;
		  Vertex[k].pos = Points[p];
		  Vertex[k].num = p;
		  points_particules[p][i] = true;
		  double x = CGAL::to_double(Points[p].operator[](0));
		  double y = CGAL::to_double(Points[p].operator[](1));
		  double z = CGAL::to_double(Points[p].operator[](2));
		  xmin = min(x,xmin);
		  xmax = max(x,xmax);
		  ymin = min(y,ymin);
		  ymax = max(y,ymax);
		  zmin = min(z,zmin);
		  zmax = max(z,zmax);
		}
		int voisin;
		maillage >> voisin;
		Faces[j] = Face::Face(Vertex, voisin);
	  }
	  
	  P[i] = Particule::Particule(xmin, ymin, zmin, xmax, ymax, zmax, Faces);
	  
	}
	//Boucle de mise a jour des particules sur les sommets du maillage
	for(int i=0;i<P.size();i++){
	  for(int j=0;j<P[i].size();j++){
		for(int k=0;k<P[i].faces[j].size();k++){
		  for(int l=0;l<P.size();l++){
			if(points_particules[P[i].faces[j].vertex[k].num][l]){
			  P[i].faces[j].vertex[k].particules.push_back(l);
			}
		  }
		}
	  }
	}
	
	
	
	
	

	
	Solide S(P);
	//S.Affiche();


	Fluide.init();
	Fluide.parois(S);
	Fluide.BC();
	

	int iter=0;
	
	clock_t start,end;
	
	start =clock();

	int kimp = 0; //Numero de suivi de l'impression
	double next_timp = dtimp; //Instant de la prochaine impression
	Fluide.impression(kimp);
	S.impression(kimp);
	kimp++;
	
	
	for (int n=0; (t<T) && n<Nmax; n++){
		
	  if(t>next_timp){
		Fluide.impression(kimp);
		S.impression(kimp);
		kimp++;
		next_timp += dtimp;
	  }
		cout<<"Energie: "<<Fluide.Energie()<<"  "<<"Masse : "<<"  "<< Fluide.Masse()<<endl;
		double dt = Fluide.pas_temps(t, T);
		//Fluide.affiche("avant Solve");
		Fluide.Solve(dt, t, n);
		//Fluide.affiche("Solve");
		Fluide.modif_fnum(dt);
		//Fluide.affiche("modif_fnum");
		Fluide.mixage();
		//Fluide.affiche("mixage");
		Fluide.fill_cel(S);
		//Fluide.affiche("fill_cell");
		Fluide.BC();
		//Fluide.affiche("BC");
		out<< n << " temps actuel "<<t<<" dt "<<dt<<"\n";
		cout<<"iteration="<<n<< " dt="<<dt<<" t="<<t<<endl;
		t+= dt;
		iter++;
		
	}
	end=clock();
	
	Fluide.impression(kimp);
	S.impression(kimp);
	
	out<< "Temps final  "<< t<<endl;
	out<<"nb iter= "<< iter<<endl;    
	out <<"Temps de calcul " <<(double) (end-start)/CLOCKS_PER_SEC << endl;     
	cout<<"nb iter= "<< iter<<endl;    
	cout <<"Temps de calcul " <<(double) (end-start)/CLOCKS_PER_SEC << endl;    
	
	return 0;
	
	
}
