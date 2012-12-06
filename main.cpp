#include <iostream>
#include <ctime>
#include "fluide.cpp"
#include "couplage.cpp"
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
	
	CGAL::Timer user_time, user_time2;
	for (int n=0; (t<T) && n<Nmax; n++){
		
		if(t>next_timp){
			Fluide.impression(kimp);
			S.impression(kimp);
			kimp++;
			next_timp += dtimp;
		}
	  //cout<<"Energie: "<< Fluide.Energie()+S.Energie() << " Solide:" << S.Energie() <<"  "<<"Masse : "<<"  "<< Fluide.Masse() <<endl;
		cout<<"Energie Fluide: "<< Fluide.Energie() << " Energie Solide:" << S.Energie() <<"  "<<"Masse : "<<"  "<< Fluide.Masse() <<endl;
	  ener << t << " " << Fluide.Energie()+S.Energie() << " " << S.Energie() << " " << Fluide.Energie()+S.Energie()-E0 << " " << S.Energie()-E0S << endl;
	  //double dt = min(Fluide.pas_temps(t, T),S.pas_temps(t,T));
		double dt = 0.005;
		//Fluide.affiche("avant Solve");
		user_time2.start();
		Fluide.Solve(dt, t, n);
		cout << "Temps calcul flux: " << user_time2.time() << " seconds." << endl;
		user_time2.reset();
		//Fluide.affiche("Solve");
		S.solve_position(dt);
		S.Forces_internes();
		S.solve_vitesse(dt);
		Fluide.parois(S);
		//S.Affiche();
		user_time.start();
		int n0=0.,n1=0., m=0.;
		Fluide.swap_2d(dt,S,n0,n1,m);
		cout << "Temps swap: " << user_time.time() << " seconds." << endl;
		user_time.reset();
		cout<<"Triangles en n "<<n0<<" Triangles en n+1 "<<n1<<" Triangles sous maillage "<<m<<endl;
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
