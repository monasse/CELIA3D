#include <iostream>
#include <ctime>

//#include "fluide.hpp"
#include "fluide.cpp"
#include "couplage.cpp"
//#include "parametres.hpp"
#include "parametres.cpp"
#include "solide.cpp"


using namespace std;          // espace de nom standard


int main(){
  char tempsF[]="resultats/tempsF.dat";
	
	//Ouverture des flux en donne en ecriture
	std::ofstream out(tempsF,ios::out);
	if(out)
	{
		// cout <<"ouverture de xt.vtk reussie" << endl;
	} else {
		cout <<"ouverture de .dat rate" << endl;
	}

	char energie[]="resultats/energie.dat";
	
	//Ouverture des flux en donne en ecriture
	std::ofstream ener(energie,ios::out);
	if(ener)
	{
		// cout <<"ouverture de xt.vtk reussie" << endl;
	} else {
		cout <<"ouverture de .dat rate" << endl;
	}
	
	double t=0.;


	Solide S;
	S.init("maillage.dat"); //Initialisation du solide a partir du fichier "maillage.dat"
	//S.Affiche();

	Grille Fluide;
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
	
	double E0 = Fluide.Energie()+S.Energie();
	double E0S= S.Energie();

	S.Forces_internes();
	
	CGAL::Timer user_time;
	for (int n=0; (t<T) && n<Nmax; n++){
		
		if(t>next_timp){
			Fluide.impression(kimp);
			S.impression(kimp);
			kimp++;
			next_timp += dtimp;
		}
	  //cout<<"Energie: "<< Fluide.Energie()+S.Energie() << " Solide:" << S.Energie() <<"  "<<"Masse : "<<"  "<< Fluide.Masse() <<endl;
		cout<<"Energie: "<< Fluide.Energie() << " Solide:" << S.Energie() <<"  "<<"Masse : "<<"  "<< Fluide.Masse() <<endl;
	  ener << t << " " << Fluide.Energie()+S.Energie() << " " << S.Energie() << " " << Fluide.Energie()+S.Energie()-E0 << " " << S.Energie()-E0S << endl;
	  double dt = min(Fluide.pas_temps(t, T),S.pas_temps(t,T));
		//Fluide.affiche("avant Solve");
		Fluide.Solve(dt, t, n);
		//Fluide.affiche("Solve");
		S.solve_position(dt);
		S.Forces_internes();
		S.solve_vitesse(dt);
		Fluide.parois(S);
		//S.Affiche();
		Fluide.modif_fnum(dt);
		//Fluide.affiche("modif_fnum");
		user_time.start();
		Fluide.swap(dt,S);
		cout << "Temps swap: " << user_time.time() << " seconds." << endl;
		user_time.reset();
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
	
// 	//test 8 nov
// 	Fluide.affiche();
// 	//fin test 8 nov
	
	return 0;
	
	
}
