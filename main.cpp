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
{	
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
	const int nb_particule = 4;
	
	vector<Solide> S(nb_particule);
	



   //rampe
	Point_3 p1_s1(-1.,-1.10, -0.5), p1_r1(0.25, -1.10, -0.5), p1_t1(0.25, 0.15, -0.5), p1_v1(-1, 0.15, -0.5);
	Point_3 p1_s2(-1.,-1.10, 1.5), p1_r2(0.25, -1.10, 1.5), p1_t2(0.25, 0.15, 1.5), p1_v2(-1, 0.15, 1.5);
	S[0] = Solide(-1.,-1.10, -0.5, 0.25, 0.15, 1.5, p1_s1, p1_r1, p1_t1, p1_v1, p1_s2, p1_r2, p1_t2, p1_v2);
	//S[0].Affiche();
	Point_3 p2_s1(0.25, -1.10, -0.5), p2_r1(1., -1.10, -0.5), p2_t1(1., 0.5829, -0.5), p2_v1(0.25, 0.15, -0.5);
	Point_3 p2_s2(0.25, -1.10, 1.5), p2_r2(1., -1.10, 1.5), p2_t2(1., 0.5829, 1.5), p2_v2(0.25, 0.15, 1.5);
	S[1] = Solide(0.25, -1.10, -0.5, 1., 0.5829, 1.5,p2_s1, p2_r1, p2_t1, p2_v1, p2_s2, p2_r2, p2_t2, p2_v2);
	//S[1].Affiche();
	Point_3 p3_s1(1., -1.10, -0.5), p3_r1(2., -1.10, -0.5), p3_t1(2., 1.1602, -0.5), p3_v1(1., 0.5829, -0.5);
	Point_3 p3_s2(1., -1.10, 1.5), p3_r2(2., -1.10, 1.5), p3_t2(2., 1.1602, 1.5), p3_v2(1., 0.5829, 1.5);
	S[2] = Solide(1., -1.10, -0.5, 2., 1.1602, 1.5,p3_s1, p3_r1, p3_t1, p3_v1, p3_s2, p3_r2, p3_t2, p3_v2);
	//S[2].Affiche();
	Point_3 p4_s1(2., -1.10, -0.5), p4_r1(3.4543, -1.10, -0.5), p4_t1(3.4543, 2., -0.5), p4_v1(2., 1.1602, -0.5);
	Point_3 p4_s2(2., -1.10, 1.5), p4_r2(3.4543, -1.10, 1.5), p4_t2(3.4543, 2., 1.5), p4_v2(2., 1.1602, 1.5);
	S[3] = Solide(2., -1.10, -0.5, 3.4543, 2., 1.5,p4_s1, p4_r1, p4_t1, p4_v1, p4_s2, p4_r2, p4_t2, p4_v2);
	//S[3].Affiche();
	


	Fluide.init();
	Fluide.parois(S);
	Fluide.BC();
	

	int iter=0;
	
	clock_t start,end;
	
	start =clock();

	int kimp = 0; //Numero de suivi de l'impression
	double next_timp = dtimp; //Instant de la prochaine impression
	Fluide.impression(kimp);
	kimp++;
	
	
	for (int n=0; (t<T) && n<Nmax; n++){
		
	  if(t>next_timp){
		Fluide.impression(kimp);
		kimp++;
		next_timp += dtimp;
	  }
		cout<<"Energie: "<<Fluide.Energie()<<"  "<<"Masse : "<<"  "<< Fluide.Masse()<<endl;
		double dt = Fluide.pas_temps(t, T);
		Fluide.Solve(dt, t, n);
		Fluide.modif_fnum(dt);
		Fluide.mixage();
		Fluide.fill_cel(S);
		Fluide.BC();
		out<< n << " temps actuel "<<t<<" dt "<<dt<<"\n";
		cout<<"iteration="<<n<< " dt="<<dt<<" t="<<t<<endl;
		t+= dt;
		iter++;
		
	}
	end=clock();
	
	Fluide.impression(iter);
	
	out<< "Temps final  "<< t<<endl;
	out<<"nb iter= "<< iter<<endl;    
	out <<"Temps de calcul " <<(double) (end-start)/CLOCKS_PER_SEC << endl;     
	cout<<"nb iter= "<< iter<<endl;    
	cout <<"Temps de calcul " <<(double) (end-start)/CLOCKS_PER_SEC << endl;    
	
	return 0;
	
	
}
