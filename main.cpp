#include <iostream>
#include <ctime>
#include "fluide.cpp"
#include "solide.cpp"
#include "couplage.cpp"
#include "parametres.cpp"


/*!
 \mainpage 
 \section intro Introduction
 Code de simulation 3d d'interaction fluide-structure (CELIA_3d).
 La m&eacute;thode num&eacute;rique de couplage utilise une approche de type Volumes Finis de l'interaction fluide-structure.\n
 Le solide est mod&eacute;lis&eacute; selon une description de type &quot;&eacute;l&eacute;ments discrets" (Mka3d).\n
 Le fluide est trait&eacute; par une m&eacute;thode de Volumes Finis, le calcul des flux num&eacute;riques utilise le sch&eacute;ma OSMP (CHORUS).

 \section description Description du projet
 
 Pour utiliser le logiciel:
 - installer la <b> libraire CGAL-4.0 </b>
 - cgal_create_cmake_script
 - cmake .
 - make: compilation 
 - ./main: execution 
 
 Ce que vous devez remplir pour lancer une simulation:
 
 - \a parametres.hpp: d&eacute;finir les param&egrave;tres du probl&egrave;me; 
 - \a parametres.cpp: d&eacute;finir l'&eacute;tat initial du fluide: densit&eacute;, pression, vitesse;
 - cr&eacute;ation du fichier <b> "maillage.dat" </b>d&eacute;finissant le maillage pour le solide.\note Attention, respectez la syntaxe expliqu&eacute;e dans le fichier \b README. \n
 
 
 
 \remark Les proc&eacute;dures concernant le fluide se trouvent dans les fichiers fluide.hpp et fluide.cpp. Celles concernant le solide dans solide.hpp et solide.cpp. Les proc&eacute;dures r&eacute;alisant le couplage sont regroup&eacute;es dans les fichiers couplage.cpp, intersections.hpp et intersections.cpp. Les fichiers parametres.hpp et parametres.cpp sont d&eacute;dies aux d&eacute;finitions des param&eacute;tr&eacute;s du probl&egrave;me (param&eacute;tr&eacute;s physique, li&eacute;s aux maillages fluide et solide, au couplage, etc.). La r&eacute;solution du probl&egrave;me est r&eacute;alis&eacute;e  dans le fichier principal main.cpp. 
 
 
 \authors Maria Adela PUSCAS et Laurent MONASSE
 */
/*!
*  \file main.cpp
*  \brief Fonction principale. Initialisation du probl&egrave;me et r&eacute;solution.
*  \author M. A. PUSCAS et L. MONASSE
*  \version 1.0
*  \date 06/04/2013
*
*/


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
    //Recuperation de l'energie
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
	S.Init("maillage.dat"); //Initialisation du solide a partir du fichier "maillage.dat"
	Grille Fluide;
	Fluide.Init();
	Fluide.Parois(S, dt);
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
	  Fluide.Impression(kimp);
	  S.Impression(kimp);
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
	int nb_iter_implicit=0;
	CGAL::Timer user_time, user_time2, user_time3,user_time4,user_time5;
	double temps_flux=0., temps_solide_f_int=0., temps_couplage=0., temps_swap=0., temps_intersections=0., temps_semi_implicit=0.;
	for (int n=0; (t<T) && n<Nmax; n++){
		
		if(t>next_timp){
			Fluide.Impression(kimp);
			S.Impression(kimp);
			sorties_reprise << t << endl;
			kimp++;
			next_timp += dtimp;
		}
		cout<<"iteration="<<n<< " dt="<<dt<<" t="<<t<<endl;
		
	  //cout<<"Energie: "<< Fluide.Energie()+S.Energie() << " Solide:" << S.Energie() <<"  "<<"Masse : "<<"  "<< Fluide.Masse() <<endl;
		cout<<"Energie Fluide: "<< Fluide.Energie() << " Energie Solide:" << S.Energie() <<"  "<<"Masse : "<<"  "<< Fluide.Masse() <<endl;
		ener << t << " " << Fluide.Energie()+S.Energie() << " " << S.Energie() << " " << Fluide.Energie()+S.Energie()-E0 << " " << S.Energie()-E0S <<" "<<Fluide.Masse() - masse <<endl;
		center<< t << " "<<S.solide[nb_part-1].x0.operator[](0) + S.solide[nb_part-1].Dx.operator[](0) << " "<<S.solide[nb_part-1].x0.operator[](1) + S.solide[nb_part-1].Dx.operator[](1) << " "<<S.solide[nb_part-1].x0.operator[](2) + S.solide[nb_part-1].Dx.operator[](2) <<endl;

		
		dt = min(Fluide.pas_temps(t, T),S.pas_temps(t,T));
		//Fluide.affiche("avant Solve");
		user_time2.start();
		Fluide.Solve(dt, t, n);
		//cout << "Temps calcul flux: " << user_time2.time() << " seconds." << endl;
		temps_flux += CGAL::to_double(user_time2.time());
		user_time2.reset();
		//Fluide.affiche("Solve");
		S.Forces_internes();
		temps_solide_f_int += CGAL::to_double(user_time3.time());
		user_time3.reset();
		if(explicite){ //algo de couplage explicit
		Fluide.Forces_fluide(S,dt);
		S.Solve_position(dt);
		user_time3.start();
		S.Solve_vitesse(dt);
		user_time4.start();
		Fluide.Parois(S,dt);
		temps_intersections += CGAL::to_double(user_time4.time());
		user_time4.reset();
		//S.Affiche();
		}
		
		else{ //algo de couplage semi-implicit
			//semi-implicit
			int kmax=50;
			double erreur = 1.;
			user_time5.start();
			Solide Sk = S;
			Solide Sk_prev = Sk;
			int k;
			for(k=0;(erreur>1.e-19) && (k<kmax) ;k++){
				Fluide.Forces_fluide(Sk,dt);
				Copy_f_m(S,Sk); //attention c'est tres important d'appeller cette fonction car sinon on va ecraser les valeurs Sk.F_f et Sk.M_f!!!!!
				Sk_prev = Sk; 
				Sk = S ; 
				Sk.Solve_position(dt);
				Sk.Solve_vitesse(dt);
				Fluide.Parois(Sk,dt);
				erreur = Error(Sk, Sk_prev);
				//cout<<" erreur := "<<erreur<<endl;
			}//fin boucle 
			S=Sk;
			temps_semi_implicit += CGAL::to_double(user_time4.time());
			user_time5.reset();
			cout<<"nb iteration semi-implicit "<<k<<endl;
			nb_iter_implicit += k;
			//semi-implicit	
		}
		user_time.start();
		Fluide.Swap_2d(dt,S);
		//cout << "Temps swap: " << user_time.time() << " seconds." << endl;
		temps_swap += CGAL::to_double(user_time.time());
		user_time.reset();
		//cout<<"Triangles en n "<<n0<<" Triangles en n+1 "<<n1<<" Triangles sous maillage "<<m<<endl;
		Fluide.Modif_fnum(dt);
		//Fluide.affiche("modif_fnum");
		Fluide.Mixage();
		//Fluide.affiche("mixage");
		Fluide.Fill_cel(S);
		//Fluide.affiche("fill_cell");
		Fluide.BC();
		//Fluide.affiche("BC");
		//temps_iter<< n << " "<<t<<" "<<dt<<endl;
		t+= dt;
		iter++;

	}
	end=clock();
	
	Fluide.Impression(kimp);
	S.Impression(kimp);
	
	temps_iter<< "Temps final  "<< t<<endl;
	temps_iter<<"Nb iter= "<< iter<<endl;    
	temps_iter<<"Temps calcul " <<(double) (end-start)/CLOCKS_PER_SEC << endl; 
	temps_iter<<"Temps calcul flux " <<temps_flux<< endl; 
	temps_iter<<"Temps calcul forces internes " <<temps_solide_f_int<< endl; 
	temps_iter<<"Temps calcul swap " <<temps_swap << endl; 
	temps_iter<<"Temps calcul intersections " << temps_intersections<< endl; 
	temps_iter<<"Temps semi-implicit" << temps_semi_implicit<< endl; 
	temps_iter<<"Temps couplage " << temps_swap + temps_intersections + temps_semi_implicit<< endl; 
	cout<<"Nb iter= "<< iter<<endl;    
	cout <<"Temps calcul " <<(double) (end-start)/CLOCKS_PER_SEC << endl;  

	return 0;
}
