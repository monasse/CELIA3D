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
 - cr&eacute;ation du fichier <b> "maillage.dat" </b>d&eacute;finissant le maillage pour le solide.\n
 
 
 
 
Cr&eacute;ation du fichier "maillage.dat": \n
- Premi&egrave;re ligne le mot cl&eacute; \b "POINTS" et le nombre de sommets du solide. \n
- Lignes suivantes : l'ensemble des coordonn&eacute;es de sommets. Chaque sommet sur une ligne.\n
- Le mot cl&eacute; \b "SOLIDE" et le nombre de particules solides du syst&egrave;me .\n
- Les caract&eacute;ristiques de chaque particule: \n
 - le mot cl&eacute; \b "PARTICULE" suivi du nombre de faces et la condition sur le mouvement de la particule: 0 ou 1. Si 0, la particule peut bouger librement si 1, la particule est fix&eacute;e &agrave; sa position initiale.\n
 - le mot cl&eacute; \b "POSITION" et les coordonn&eacute;es du centre de la particule. \n
 - le mot cl&eacute; \b "VITESSE" et les composantes de la vitesse initiale de la particule. \n
 - le mot cl&eacute; \b "VITROT"  et les composantes de la vitesse de rotation initiale de la particule. \n
- Description g&eacute;om&eacute;trique de chacune des particules face par face: 
 - nombre de sommets de la face, les num&eacute;ros de sommets (dans le m&ecirc;me ordre que dans le listing 
 des coordonn&eacute;es de sommets) d&eacute;crivant la face et un dernier nombre, part, indiquent le num&eacute;ro 
 de la particule (dans le m&ecirc;me ordre que dans le listing de particules) &eacute;ventuellement 
 pr&eacute;sente de l'autre côt&eacute; de la face. Si aucune particule n'est pr&eacute;sente de l'autre côt&eacute; 
 de la face (donc si la particule est en contact avec le fluide), on a part =-1. \n
 
\remark L'ordre de ces sommets donne l'orientation de la normale sortante &agrave; la particule,  
 les sommets sont num&eacute;rot&eacute;s de 0 &agrave; nombre de sommets-1 et les particules sont 
 num&eacute;rot&eacute;s de 0 &agrave; nombre de particules-1 comme dans un tableau C++. \n
 
 
Exemple: \n
 
			POINTS 8 \n
			0.4 0.4 0.4 \n
			0.9 0.4 0.4 \n
			0.9 0.6 0.4 \n
			0.4 0.6 0.4 \n
			0.4 0.4 0.6 \n
			0.9 0.4 0.6 \n
			0.9 0.6 0.6 \n
			0.4 0.6 0.6 \n
			
			SOLIDE 1 \n
			
			PARTICULE 6 0 \n
			POSITION  1 1 1 \n
			VITESSE 1. 0. 0.\n
			VITROT 0. 0. 0.\n
			4 0 3 2 1 -1 \n
			4 7 4 5 6 -1 \n
			4 0 4 7 3 -1 \n
			4 1 2 6 5 -1 \n
			4 4 0 1 5 -1 \n
			4 3 7 6 2 -1 \n
 
 
 Les r&eacute;sultats sont &eacute;crits dans le dossier \b resultats. Certains fichiers sont remplis &agrave; 
 chaque pas de temps :
 les fichiers energie.dat, solide_center.dat et temps_reprise.dat donnent respectivement l'&eacute;volution
 de l'&eacute;nergie, la position du centre du solide au cours du temps et le temps de la simulation. 
 Les fichiers fluide*.vtk et solide*.vtk sont &eacute;crits un nombre limit&eacute;
 de fois au cours du calcul; fluide*.vtk et solide*.vtk donnent l'&eacute;tat du fluide et la position du solide,
 peuvent &ecirc;tre lus avec Paraview. Le fichier temps.dat est &eacute;crit &agrave; la fin de la simulation 
 donnent le temps de calcul (temps machine) des differents sous-proc&eacute;dures  du code. \n
 Il est possible d'effectuer une reprise sur un calcul interrompu &agrave; partir des sorties fluide*.vtk
 et solide*.vtk: il suffit de changer le flag de reprise bool rep=false en bool rep=true dans le fichier parametres.h
 et d'indiquer &agrave; partir de quel point de reprise on souhaite reprendre &agrave; l'aide de la variable int numrep.
 Par exemple, pour reprendre avec fluide1.vtk et solide1.vtk (c'est-a-dire &agrave; partir de la premi&egrave;re 
 sortie), on utilisera int numrep=1.
 
 
 \remark Les proc&eacute;dures concernant le fluide se trouvent dans les fichiers fluide.hpp et fluide.cpp. Celles concernant le solide dans solide.hpp et solide.cpp. Les proc&eacute;dures r&eacute;alisant le couplage sont regroup&eacute;es dans les fichiers couplage.cpp, intersections.hpp et intersections.cpp. Les fichiers parametres.hpp et parametres.cpp sont d&eacute;dies aux d&eacute;finitions des param&eacute;tr&eacute;s du probl&egrave;me (param&eacute;tr&eacute;s physique, li&eacute;s aux maillages fluide et solide, au couplage, etc.). La r&eacute;solution du probl&egrave;me est r&eacute;alis&eacute;e  dans le fichier principal main.cpp. 
 
 
 \authors Maria Adela PUSCAS et Laurent MONASSE
 */


/*!
*  \file main.cpp
*  \brief Fonction principale. Initialisation du probl&egrave;me et r&eacute;solution.
*/

#include <iostream>
#include <ctime>
#include "fluide.cpp"
#include "solide.cpp"
#include "couplage.cpp"
#include "parametres.cpp"

using namespace std;          // espace de nom standard

/*!
* \fn int main()
*\brief Initialisation du probl&egrave;me et r&eacute;solution:

- Initialisation du solide via la fonction \a Solide.Init(const char*) et du fluide via la fonction  Grille.Init().
- R&eacute;solution du probl&egrave;me:
 - R&eacute;solution des &eacute;quations fluides via la fonction \a Grille.Solve(const double, double, int).
 - Calcul des forces internes via la fonction \a Solide.Forces_internes().
 -  Calcul des forces (\a Particule.Ff) et moments fluides (\a Particule.Mf) exerc&eacute;s sur le solide via la fonction \a Grille.Forces_fluide(Solide&, const double).
 - Mise &agrave; jour de la position du solide via la fonction \a Solide.Solve_position(double).
 - Calcul de la vitesse du solide via la fonction \a Solide.Solve_vitesse(double).
 - Intersection de la grille fluide avec le solide via la fonction \a Grille.Parois(Solide&, double).
 - Calcul de la quantit&eacute; balay&eacute;e par le solide via la fonction \a Grille.Swap_2d(double, Solide&).
 - Modification des flux fluide via la fonction \a Grille.Modif_fnum(double).
 - M&eacute;lange conservatif de petites cellules coup&eacute;es via la fonction \a Grille.Mixage().
 - Remplissage des cellules fictives via la fonction \a Grille.Fill_cel(Solide&).
 - Imposition des conditions aux limites via la fonction \a Grille.BC().
 
*\return int
*/
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
