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
#include "fluide1d.hpp"
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
	char porte_rot_e[]="resultats/porte_e.dat"; 
	char porte_rot_w[]="resultats/porte_w.dat"; 
	char door_press[]="resultats/door.dat"; 
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
	std::ofstream porte_e(porte_rot_e,ios::out);
	std::ofstream porte_w(porte_rot_w,ios::out);
	std::ofstream door(door_press,ios::out);
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
	  double E,Es,dE,dEs,dm, tmp;
	  for(int i=0;t_ener<temps[numrep];i++){
	    in_energie >>t_ener >> E >> Es >> dE >> dEs >> dm >> tmp;
	    if(t_ener<=temps[numrep]){
	      ener << t_ener << " " << E << " " << Es << " " << dE << " " << dEs << " " << dm << 0 << endl;
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
	//S.Affiche();
	Grille Fluide;
	Fluide.Init();
	Fluide.Parois_particles(S,dt);
	//Fluide.parois_cellule_vide(S); //test 25 octobre 2013
	//Fluide.Fill_cel(S); //test 4 nov 2013
	Fluide.BC();
	
	double volume_initial= 0.;
	
	for (int count=0; count<S.size(); count++){
		volume_initial +=S.solide[count].volume();
		//cout<<"position du centre de la particule "<<S.solide[count].x0 + S.solide[count].Dx<<endl;
		
	}
	//probleme 1D
	ofstream xt("resultats/xt.dat",ios::out | ios::trunc);
	Grille1D Fluide1D;
	Fluide1D.init();
	Fluide1D.cond_lim();
	//Fluide1D.affiche();
	//fin probleme1D
	
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
	CGAL::Timer user_time, user_time2, user_time3,user_time4,user_time5,user_time_total;
	double temps_flux=0., temps_solide_f_int=0., temps_couplage=0., temps_swap=0., temps_intersections=0., temps_semi_implicit=0., temps_explicit=0., temps_solide_vitesse=0.,temps_modif_fnum=0.,temps_mixage=0.,temps_fill_cel=0.,temps_BC=0.,temps_total=0.;
	double variation_masse= 0.;
	double variation_energy= 0.;
	double variation_volume = 0.;
	double volume_solide = 0.;
	//double flux=0.;
	
	for (int n=0; (t<T) && n<Nmax; n++){
		user_time_total.start();

		cout<<"iteration="<<n<< " dt="<<dt<<" t="<<t<<endl;
		//probleme1D
		if(t<=0.0075 && couplage1d){
			//Fluide1D.affiche();
			dt =Fluide1D.pas_temps(t, T);
			Fluide1D.solve_fluid( dt,t);
			Fluide1D.cond_lim();
			Fluide1D.impression_1d(t,n);
			t+= dt;
			iter++;
			//Sortie des parametres du fluide vers gnuplot
			ofstream fluide_temp("resultats/fluide1d_temp.txt",ios::out | ios::trunc);
			for(int i=marge; i<N+marge; i++){
				Cellule1D c = Fluide1D.grille[i];
				fluide_temp << c.x << " " << t << " " << c.rho << " " << c.u << " " << c.p <<endl;
			}
			
// 			//Ouverture gnuplot
// 			char* path1;
// 			char* path2;
// 			
// 			path1 = "/usr/bin/gnuplot";
// 			path2 = "/u/cermics/l/puscasa/Bureau/ouverture_porte/resultats/fluide1d_temp.txt";
// 			FILE* gp = popen(path1,"w");
// 			if(gp == NULL){
// 				cout << "Oops, I can't find %s." << "gnuplot" << endl;
// 				exit(EXIT_FAILURE);
// 				getchar();
// 			}
// 			
// 			fprintf(gp,"set xrange [0:4];\n");
// 			fprintf(gp,"p '");
// 			fprintf(gp,path2);
// 			fprintf(gp,"' u 1:($3)");
// 			fprintf(gp," w l;\n");
// 			fflush(gp);
			
		}
		else{
			
			if(t>next_timp){
			  Fluide.Impression(kimp);
			  S.Impression(kimp);
			  if(couplage1d){
			    Fluide1D.impression(t, kimp);
			  }
			  sorties_reprise << t << endl;
			  kimp++;
			  next_timp += dtimp;
			}
			cout<<"Energie Fluide: "<< Fluide.Energie() << " Energie Solide:" << S.Energie() <<"  "<<"Masse : "<<"  "<< Fluide.Masse() <<endl;
			ener << t << " " << Fluide.Energie()+S.Energie() << " " << S.Energie() << " " << Fluide.Energie()+S.Energie()-E0 << " " << S.Energie()-E0S <<" "<<Fluide.Masse() - masse <<endl;
			center<< t << " "<<S.solide[nb_part-1].x0.operator[](0) + S.solide[nb_part-1].Dx.operator[](0) << " "<<S.solide[nb_part-1].x0.operator[](1) + S.solide[nb_part-1].Dx.operator[](1) << " "<<S.solide[nb_part-1].x0.operator[](2) + S.solide[nb_part-1].Dx.operator[](2) <<endl;
			cout<<"Variation Energie: "<< Fluide.Energie() +  S.Energie() - E0<<" Variation Masse : "<< Fluide.Masse() - masse<<endl;
			//Pas de temps
			dt = min(Fluide.pas_temps(t, T),S.pas_temps(t,T));
			//Fluide.affiche("avant Solve");
			//cout<<"avant solve"<<endl;
			Fluide.affiche();
			user_time2.start();
			//Fluide1D.solve_fluid( dt,t);
			//Fluide1D.cond_lim();
			double tab[marge][3];
			if(couplage1d){
			  for(int iter_tab=0; iter_tab<marge; iter_tab++){
			    tab[iter_tab][0] = Fluide1D.grille[N+iter_tab].rho;
			    tab[iter_tab][1] = Fluide1D.grille[N+iter_tab].imp;
			    tab[iter_tab][2] = Fluide1D.grille[N+iter_tab].rhoE;
			  }
			}
			//cout<<" cd limites 1d pour la grille 3d avant solve "<<tab[0][2]<< " "<<tab[1][2]<< " "<<tab[2][2]<< " "<<tab[3][2]<<" "<<tab[4][2]<<" "<<tab[5][2]<<endl;
			cout <<"BC_couplage Variation Masse : "<< Fluide.Masse() - masse<<endl;
			if(couplage1d){
			  Fluide.BC_couplage(tab,couplage1d);
			  //conditions aux limites couplage 1D pour la grille 1d
			  vector< vector < double> > tab_1d;
			  tab_1d.resize(marge, std::vector<double>(0));
			  for(int iter_tab=0; iter_tab<marge; iter_tab++){
			    tab_1d[iter_tab].resize(3, 0.);
			  }
			  cout <<"BC_couplage_1d Variation Masse : "<< Fluide.Masse() - masse<<endl;
			  Fluide.BC_couplage_1d(tab_1d);
			  //cout<<" cd limites 3d pour la grille 1d "<<tab_1d[0][2]<< " "<<tab_1d[1][2]<< " "<<tab_1d[2][2]<< " "<<tab_1d[3][2]<<" "<<tab_1d[4][2]<<" "<<tab_1d[5][2]<<endl;
			  cout <<"Fluide1D.cond_lim_couplage Variation Masse : "<< Fluide.Masse() - masse<<endl;
			  Fluide1D.cond_lim_couplage(tab_1d);
			}
			//fin conditions aux limites couplage 1D pour la grille 1d
			//cout<<" cd limites 1d pour la grille 3d avant solve "<<tab[0][0]<< " "<<tab[1][0]<< " "<<tab[2][0]<< " "<<tab[3][0]<<" "<<tab[4][0]<<" "<<tab[5][0]<<endl;
			//cout << "x fluide3d=" << Fluide.grille[marge][23][23].x << " " << Fluide.grille[marge+1][23][23].x << " " << Fluide.grille[marge+2][23][23].x << " " << Fluide.grille[marge+3][23][23].x << " " << Fluide.grille[marge+4][23][23].x << " " << Fluide.grille[marge+5][23][23].x << endl;
			//cout << "y fluide3d=" << Fluide.grille[marge][23][23].y << " " << Fluide.grille[marge+1][23][23].y << " " << Fluide.grille[marge+2][23][23].y << " " << Fluide.grille[marge+3][23][23].y << " " << Fluide.grille[marge+4][23][23].y << " " << Fluide.grille[marge+5][23][23].y << endl;
			//cout << "z fluide3d=" << Fluide.grille[marge][23][23].z << " " << Fluide.grille[marge+1][23][23].z << " " << Fluide.grille[marge+2][23][23].z << " " << Fluide.grille[marge+3][23][23].z << " " << Fluide.grille[marge+4][23][23].z << " " << Fluide.grille[marge+5][23][23].z << endl;
			//cout << "p marge fluide 1d=" << Fluide1D.grille[N+marge].p << " " << Fluide1D.grille[N+marge+1].p << " " << Fluide1D.grille[N+marge+2].p << " " << Fluide1D.grille[N+marge+3].p << " " << Fluide1D.grille[N+marge+4].p << " " << Fluide1D.grille[N+marge+5].p << endl;
			//cout << "p       fluide 3d=" << Fluide.grille[marge][23][23].p << " " << Fluide.grille[marge+1][23][23].p << " " << Fluide.grille[marge+2][23][23].p << " " << Fluide.grille[marge+3][23][23].p << " " << Fluide.grille[marge+4][23][23].p << " " << Fluide.grille[marge+5][23][23].p << endl;
			//cout << "p marge fluide 3d=" << Fluide.grille[0][23][23].p << " " << Fluide.grille[1][23][23].p << " " << Fluide.grille[2][23][23].p << " " << Fluide.grille[3][23][23].p << " " << Fluide.grille[4][23][23].p << " " << Fluide.grille[5][23][23].p << endl;
			//cout << "p       fluide 1d=" << Fluide1D.grille[N].p << " " << Fluide1D.grille[N+1].p << " " << Fluide1D.grille[N+2].p << " " << Fluide1D.grille[N+3].p << " " << Fluide1D.grille[N+4].p << " " << Fluide1D.grille[N+5].p << endl;
			//cout << "rho marge fluide 1d=" << Fluide1D.grille[N+marge].rho << " " << Fluide1D.grille[N+marge+1].rho << " " << Fluide1D.grille[N+marge+2].rho << " " << Fluide1D.grille[N+marge+3].rho << " " << Fluide1D.grille[N+marge+4].rho << " " << Fluide1D.grille[N+marge+5].rho << endl;
			//cout << "rho       fluide 3d=" << Fluide.grille[marge][23][23].rho << " " << Fluide.grille[marge+1][23][23].rho << " " << Fluide.grille[marge+2][23][23].rho << " " << Fluide.grille[marge+3][23][23].rho << " " << Fluide.grille[marge+4][23][23].rho << " " << Fluide.grille[marge+5][23][23].rho << endl;
			//cout << "rho marge fluide 3d=" << Fluide.grille[0][23][23].rho << " " << Fluide.grille[1][23][23].rho << " " << Fluide.grille[2][23][23].rho << " " << Fluide.grille[3][23][23].rho << " " << Fluide.grille[4][23][23].rho << " " << Fluide.grille[5][23][23].rho << endl;
			//cout << "rho       fluide 1d=" << Fluide1D.grille[N].rho << " " << Fluide1D.grille[N+1].rho << " " << Fluide1D.grille[N+2].rho << " " << Fluide1D.grille[N+3].rho << " " << Fluide1D.grille[N+4].rho << " " << Fluide1D.grille[N+5].rho << endl;
			/*for(int i=0;i<Nx+2*marge;i++){
			  for(int j=1;j<Ny+2*marge;j++){
			    for(int k=1;k<Nz+2*marge;k++){
			      Cellule c = Fluide.grille[i][j][k];
			      Cellule cg = Fluide.grille[i][j-1][k];
			      Cellule cb = Fluide.grille[i][j][k-1];
			      if(c.y<0.2 && c.y>0.1 && c.z>0.083 && c.z<0.2 && abs(c.alpha-1.)>eps){
				if(cg.y<0.2 && cg.y>0.1 && cg.z>0.083 && abs(cg.alpha-1.)>eps && abs(cg.p-c.p)>0.1){
				  //cout << "cellules differentes cg=" << cg.x << " " << cg.y << " " << cg.z << " " << cg.alpha << " " << cg.p << " c=" << c.x << " " << c.y << " " << c.z << " " << c.alpha << " " << c.p << endl;
				  //getchar();
				}
				if(cb.y<0.2 && cb.y>0.1 && cb.z>0.083 && abs(cb.alpha-1.)>eps && abs(cb.p-c.p)>0.1){
				  //cout << "cellules differentes cb=" << cb.x << " " << cb.y << " " << cb.z << " " << cb.alpha << " " << cb.p << " c=" << c.x << " " << c.y << " " << c.z << " " << c.alpha << " " << c.p << endl;
				  //getchar();
				}
			      }
			      
			    }
			  }
			  }*/
			
			cout <<"Fluide.Solve Variation Masse : "<< Fluide.Masse() - masse<<endl;
			Fluide.Solve(dt, t, n, tab, couplage1d, S);
			if(couplage1d){
			  cout <<"Fluide1D.solve_fluid Variation Masse : "<< Fluide.Masse() - masse<<endl;
			  Fluide1D.solve_fluid(dt,t);
			  cout <<"Fluide1D.cond_lim Variation Masse : "<< Fluide.Masse() - masse<<endl;
			  Fluide1D.cond_lim();
			}
			//cout << "Temps calcul flux: " << user_time2.time() << " seconds." << endl;
			//cout<<" cd limites 1d pour la grille 3d apres solve "<<tab[0][2]<< " "<<tab[1][2]<< " "<<tab[2][2]<< " "<<tab[3][2]<<" "<<tab[4][2]<<" "<<tab[5][2]<<endl;
			//cout<<"apres solve"<<endl;
			Fluide.affiche();
			temps_flux += CGAL::to_double(user_time2.time());
			user_time2.reset();
			if(explicite){ //algo de couplage explicit
			  user_time4.start();
			  cout <<"Forces_fluide Variation Masse : "<< Fluide.Masse() - masse<<endl;
			  Fluide.Forces_fluide(S,dt);
			  cout <<"Solve_position Variation Masse : "<< Fluide.Masse() - masse<<endl;
			  S.Solve_position(dt);
			  temps_explicit += CGAL::to_double(user_time4.time());
			  user_time4.reset();
			  //S.Solve_vitesse(dt);
			  // 		S.Affiche();
			  user_time4.start();
			  cout <<"Parois_particles Variation Masse : "<< Fluide.Masse() - masse<<endl;
			  Fluide.Parois_particles(S,dt);
			  //Fluide.parois_cellule_vide(S);
			  temps_intersections += CGAL::to_double(user_time4.time());
			  user_time4.reset();
			}
			
			else{ //algo de couplage semi-implicit
				//semi-implicit
				int kmax=20;
				double erreur = 1.;
				user_time5.start();
				Solide Sk = S;
				Solide Sk_prev = Sk;
				int k;
				for(k=0;(erreur>1.e-10) && (k<kmax) ;k++){
				  cout <<"Forces_fluide semiimpl Variation Masse : "<< Fluide.Masse() - masse<<endl;
				  Fluide.Forces_fluide(Sk,dt);
				  Copy_f_m(S,Sk); //attention c'est tres important d'appeller cette fonction car sinon on va ecraser les valeurs Sk.F_f et Sk.M_f!!!!!
				  Sk_prev = Sk; 
				  Sk = S ;
				  cout <<"Solve_position semiimpl Variation Masse : "<< Fluide.Masse() - masse<<endl;
				  Sk.Solve_position(dt);
				  //S.Solve_vitesse(dt);
				  cout <<"Parois_particles semiimpl Variation Masse : "<< Fluide.Masse() - masse<<endl;
				  Fluide.Parois_particles(Sk,dt);
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
			user_time3.start();
			cout <<"Forces_internes Variation Masse : "<< Fluide.Masse() - masse<<endl;
			S.Forces_internes();
			temps_solide_f_int += CGAL::to_double(user_time3.time());
			user_time3.reset();
			cout <<"Solve_vitesse Variation Masse : "<< Fluide.Masse() - masse<<endl;
			user_time3.start();
			S.Solve_vitesse(dt);
			temps_solide_vitesse += CGAL::to_double(user_time3.time());
			user_time3.reset();
			user_time.start();
			cout <<"Swap_2d Variation Masse : "<< Fluide.Masse() - masse<<endl;
			Fluide.Swap_2d(dt,S);
			//cout<<"apres swap"<<endl;
			//Fluide.affiche();
			temps_swap += CGAL::to_double(user_time.time());
			user_time.reset();
			cout <<"Modif_fnum Variation Masse : "<< Fluide.Masse() - masse<<endl;
			user_time.start();
			Fluide.Modif_fnum(dt);
			temps_modif_fnum += CGAL::to_double(user_time.time());
			user_time.reset();
			//cout<<"apres modif flux"<<endl;
			//Fluide.affiche();
			//cout<<"Masse : "<<"  "<< Fluide.Masse() - masse <<endl;
			//Fluide.Mixage();
			//Fluide.Mixage_cible(); //test 25 nov. 2013
			user_time.start();
			bool test_fini = false;
			for(int count=1;!test_fini && count<100;count++){
			  cout <<"Mixage_cible2 iteration " << count << " Variation Masse : "<< Fluide.Masse() - masse<<endl;
			  test_fini = Fluide.Mixage_cible2(); //test 14 mars 2014
			  cout << "iterations de Mixage_cible2=" << count << endl;
			}
			temps_mixage += CGAL::to_double(user_time.time());
			user_time.reset();
			//cout<<"Masse : "<<"  "<< Fluide.Masse() - masse <<endl;
			//Fluide.affiche("mixage");
			user_time.start();
			cout <<"Fill_cel Variation Masse : "<< Fluide.Masse() - masse<<endl;
			Fluide.Fill_cel(S);
			temps_fill_cel += CGAL::to_double(user_time.time());
			user_time.reset();
			//Fluide.affiche("fill_cell");
			//cout<<"apres fill"<<endl;
			//Fluide.affiche();
			user_time.start();
			cout <<"BC Variation Masse : "<< Fluide.Masse() - masse<<endl;
			Fluide.BC();
			cout <<"BC_couplage Variation Masse : "<< Fluide.Masse() - masse<<endl;
			Fluide.BC_couplage(tab,couplage1d);
			temps_BC += CGAL::to_double(user_time.time());
			user_time.reset();
			
			//cout<<"apres BC"<<endl;
			//Fluide.affiche();
			//Fluide.affiche("BC");
			//temps_iter<< n << " "<<t<<" "<<dt<<endl;
			t+= dt;
			iter++;
			variation_masse += Fluide.Masse() - masse;
			variation_energy += Fluide.Energie()+S.Energie()-E0;
			volume_solide = 0.;
			for (int count=0; count<S.size(); count++){
				volume_solide +=S.solide[count].volume();
			}
			//Fluide.affiche();
			cout<<"volume solide particules "<<volume_solide<<endl;
			variation_volume += volume_solide - volume_initial;
		}

		
		temps_total += CGAL::to_double(user_time_total.time());
		user_time_total.reset();
		cout << "############## COUTS #################" << endl;
		cout << "Fluide Solve=    " << 100*temps_flux/temps_total << "%     " << temps_flux/(n+1.) << endl;
		cout << "Solide solve=    " << 100*temps_explicit/temps_total << "%     " << temps_explicit/(n+1.)  << endl;
		cout << "Parois particles=" << 100*temps_intersections/temps_total << "%     " << temps_intersections/(n+1.) << endl;
		cout << "Forces internes= " << 100*temps_solide_f_int/temps_total << "%     " << temps_solide_f_int/(n+1.) << endl;
		cout << "Solve vitesse=   " << 100*temps_solide_vitesse/temps_total << "%     " << temps_solide_vitesse/(n+1.) << endl;
		cout << "Swap2d=          " << 100*temps_swap/temps_total << "%     " << temps_swap/(n+1.) << endl;
		cout << "Modif fnum=      " << 100*temps_modif_fnum/temps_total << "%     " << temps_modif_fnum/(n+1.) << endl;
		cout << "Mixage=          " << 100*temps_mixage/temps_total << "%     " << temps_mixage/(n+1.) << endl;
		cout << "Fill cel=        " << 100*temps_fill_cel/temps_total <<"%     " << temps_fill_cel/(n+1.) << endl;
		cout << "BC=              " << 100*temps_BC/temps_total << "%     " << temps_BC/(n+1.) << endl;
		cout << "Reste=           " << 100-100*(temps_flux+temps_explicit+temps_intersections+temps_solide_f_int+temps_solide_vitesse+temps_swap+temps_modif_fnum+temps_mixage+temps_fill_cel+temps_BC)/temps_total << "%     " << (temps_total-(temps_flux+temps_explicit+temps_intersections+temps_solide_f_int+temps_solide_vitesse+temps_swap+temps_modif_fnum+temps_mixage+temps_fill_cel+temps_BC))/(n+1.) << endl;
		cout << "########################################" << endl;
		
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
	temps_iter<<" variation masse "<< variation_masse<<endl;
	temps_iter<<" variation energy "<<variation_energy<<endl;
	temps_iter<<" variation volume "<<variation_volume<<endl;
	cout<<"Nb iter= "<< iter<<endl;    
	cout <<"Temps calcul " <<(double) (end-start)/CLOCKS_PER_SEC << endl;  
	cout<<" variation masse "<< variation_masse<<endl;
	cout<<" variation energy "<<variation_energy<<endl;
	cout<<" variation volume "<<variation_volume<<endl;

	cout << "Fluide Solve=" << 100*temps_flux/temps_total << "%" << endl;
	cout << "Solide solve=" << 100*temps_explicit/temps_total << "%"  << endl;
	cout << "Parois particles=" << 100*temps_intersections/temps_total << "%" << endl;
	cout << "Forces internes=" << 100*temps_solide_f_int/temps_total << "%" << endl;
	cout << "Solve vitesse=" << 100*temps_solide_vitesse/temps_total << "%" << endl;
	cout << "Swap2d=" << 100*temps_swap/temps_total << "%" << endl;
	cout << "Modif fnum=" << 100*temps_modif_fnum/temps_total << "%" << endl;
	cout << "Mixage=" << 100*temps_mixage/temps_total << "%" << endl;
	cout << "Fill cel=" << 100*temps_fill_cel/temps_total <<"%" << endl;
	cout << "BC=" << 100*temps_BC/temps_total << "%" << endl;
	cout << "Reste=" << 100-100*(temps_flux+temps_explicit+temps_intersections+temps_solide_f_int+temps_solide_vitesse+temps_swap+temps_modif_fnum+temps_mixage+temps_fill_cel+temps_BC)/temps_total << "%" << endl;
	

	return 0;
}
