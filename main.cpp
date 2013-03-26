#include <iostream>
#include <ctime>
#include "fluide.cpp"
#include "solide.cpp"
#include "couplage.cpp"
#include "parametres.cpp"

using namespace std;          // espace de nom standard

int main(){
  char temps_it[]="resultats/temps.dat";
  char temps_reprise[]="resultats/temps_reprise.dat";
  char solide_center[]="resultats/solide_center.dat"; 
  //En cas de reprise
  double temps[numrep+1];
  if(rep){
    reprise();
    std::ifstream in(temps_reprise,ios::in);
    if(!in){
      cout <<"ouverture de temps_reprise.dat rate" << endl;
    }
    for(int i=0;i<numrep+1;i++){
      in >> temps[i];
    }
    //R�cup�ration de l'energie
    int result = system("cp resultats/energie.dat resultats/energie_reprise.dat");
  }
  std::ifstream in_energie("resultats/energie_reprise.dat",ios::in);
  
	
	//Ouverture des flux en donne en ecriture
	std::ofstream temps_iter(temps_it,ios::out);
	std::ofstream sorties_reprise(temps_reprise,ios::out);
	std::ofstream center(solide_center,ios::out);
	if(temps_iter)
	{
		// cout <<"ouverture de 'temps.dat' reussie" << endl;
	} else {
		cout <<"ouverture de 'temps.dat' rate" << endl;
	}
	if(center)
	{
		// cout <<"ouverture de 'temps.dat' reussie" << endl;
	} else {
		cout <<"ouverture de 'solide_center.dat' rate" << endl;
	}
	if(rep){
	  for(int i=0;i<numrep+1;i++){
	    sorties_reprise << temps[i] << endl;
	  }
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
	double dE0rep,dE0Srep,dm0;
	if(rep){
	  double t_ener = 0.;
	  double E,Es,dE,dEs,dm;
	  for(int i=0;t_ener<temps[numrep];i++){
	    in_energie >>t_ener >> E >> Es >> dE >> dEs >> dm;
	    if(t_ener<=temps[numrep]){
	      ener << t_ener << " " << E << " " << Es << " " << dE << " " << dEs << " " << dm << endl;
	      dE0rep = dE;
	      dE0Srep = dEs;
	      dm0 = dm;
	    }
	  }
	  int result = system("rm resultats/energie_reprise.dat");
	}
	
	
	double t=0., dt=0.;
	if(rep){
	  t = temps[numrep];
	}
	Solide S;
	S.init("maillage4.dat"); //Initialisation du solide a partir du fichier "maillage.dat"
	Grille Fluide;
	Fluide.init();
	Fluide.parois(S, dt);
	Fluide.BC();

	int iter=0;	
	clock_t start,end;
	start =clock();

	int kimp = 0; //Numero de suivi de l'impression
	double next_timp = dtimp; //Instant de la prochaine impression
	if(rep){
	  kimp = numrep;
	  next_timp = t+dtimp;
	} else {
	  Fluide.impression(kimp);
	  S.impression(kimp);
	  sorties_reprise << t << endl;
	}
	kimp++;
	
	double E0 = Fluide.Energie()+S.Energie();
	double E0S= S.Energie();
	if(rep){
	  E0 -= dE0rep;
	  E0S -= dE0Srep;
	}
	double masse = Fluide.Masse();
	if(rep){
	  masse -= dm0;
	}
	S.Forces_internes();
	int nb_part = S.size();
	CGAL::Timer user_time, user_time2;
	for (int n=0; (t<T) && n<Nmax; n++){
		
		if(t>next_timp){
			Fluide.impression(kimp);
			S.impression(kimp);
			sorties_reprise << t << endl;
			kimp++;
			next_timp += dtimp;
		}
	  //cout<<"Energie: "<< Fluide.Energie()+S.Energie() << " Solide:" << S.Energie() <<"  "<<"Masse : "<<"  "<< Fluide.Masse() <<endl;
		cout<<"Energie Fluide: "<< Fluide.Energie() << " Energie Solide:" << S.Energie() <<"  "<<"Masse : "<<"  "<< Fluide.Masse() <<endl;
		ener << t << " " << Fluide.Energie()+S.Energie() << " " << S.Energie() << " " << Fluide.Energie()+S.Energie()-E0 << " " << S.Energie()-E0S <<" "<<Fluide.Masse() - masse <<endl;
		center<< t << " "<<S.solide[nb_part-1].x0.operator[](0) + S.solide[nb_part-1].Dx.operator[](0) << " "<<S.solide[nb_part-1].x0.operator[](1) + S.solide[nb_part-1].Dx.operator[](1) << " "<<S.solide[nb_part-1].x0.operator[](2) + S.solide[nb_part-1].Dx.operator[](2) <<endl;

		
		dt = min(Fluide.pas_temps(t, T),S.pas_temps(t,T));
		//Fluide.affiche("avant Solve");
		user_time2.start();
		Fluide.Solve(dt, t, n);
		cout << "Temps calcul flux: " << user_time2.time() << " seconds." << endl;
		user_time2.reset();
		//Fluide.affiche("Solve");
		Fluide.Forces_fluide(S,dt);
		S.solve_position(dt);
		S.Forces_internes();
		S.solve_vitesse(dt);
		Fluide.parois(S,dt);
		//S.Affiche();
		user_time.start();
		int n0=0.,n1=0., m=0.;
		Fluide.swap_2d(dt,S,n0,n1,m);
		cout << "Temps swap: " << user_time.time() << " seconds." << endl;
		user_time.reset();
		//cout<<"Triangles en n "<<n0<<" Triangles en n+1 "<<n1<<" Triangles sous maillage "<<m<<endl;
		Fluide.modif_fnum(dt);
		//Fluide.affiche("modif_fnum");
		Fluide.mixage();
		//Fluide.affiche("mixage");
		Fluide.fill_cel(S);
		//Fluide.affiche("fill_cell");
		Fluide.BC();
		//Fluide.affiche("BC");
		temps_iter<< n << " "<<t<<" "<<dt<<endl;
		cout<<"iteration="<<n<< " dt="<<dt<<" t="<<t<<endl;
		t+= dt;
		iter++;

	}
	end=clock();
	
	Fluide.impression(kimp);
	S.impression(kimp);
	
	temps_iter<< "Temps final  "<< t<<endl;
	temps_iter<<"nb iter= "<< iter<<endl;    
	temps_iter<<"Temps de calcul " <<(double) (end-start)/CLOCKS_PER_SEC << endl;     
	cout<<"nb iter= "<< iter<<endl;    
	cout <<"Temps de calcul " <<(double) (end-start)/CLOCKS_PER_SEC << endl;  

	return 0;
}
