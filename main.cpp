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

/*!
  \mainpage 
  \section intro Introduction
  Simulation code for 3D fluid-structure interaction (CELIA3D).\n
  The numerical coupling method uses a Finite Volume cut-cell approach for fluid-structure interaction.\n
  The solid is modeled with a Discrete Element Method type formulation (Mka3d).\n
  The fluid is handled with a Finite Volume method. The numerical flux scheme is scheme OSMP.

  \section description Description of the project
 
  Using the software:
  - install <b> library CGAL-4.0 </b>
  - cgal_creat_cmake_script
  - cmake .
  - make: compile
  - ./main: execute 
 
  Parameters to be filled in before launching a simulation:
 
  - \a parametres.hpp: define the parameters of the problem; 
  - \a parametres.cpp: define the initial state of the fluid: density, pressure, velocity;
  - create file <b> "maillage.dat" </b> to define the solid mesh.\n
 
 
 
 
  Creation of file "maillage.dat": \n
  - First line: key word \b "POINTS"  and the number of solid vertices. \n
  - Following lines: vertices coordinates. one line for one vertex.\n
  - Key word \b "SOLIDE"  and the number of solid particles.\n
  - Characteristics for each particle: \n
  - Key word \b "PARTICULE" followed by the number of faces and a flag indicating whether the particle is free to move (0) or fixed to its initial position (1).\n
  - Key word \b "POSITION"  and coordinates of the center of mass of the particle. \n
  - Key word \b "VITESSE" and components of the initial velocity of the particle. \n
  - Key word \b "VITROT" and components of the initial angular momentum of the particle. \n
  - Geometric description  of each particle face: 
  - Number of face vertices, list of indices of the vertices following direct orientation with regard to the exterior normal vector to the face, and index of the particle sharing the face (in case no particle shares the face and the face is a fluid interface, this number is set to -1). \n
 
  \remark Vertex indices start at 0 as in C++. \n
 
 
  Example for a cube: \n
 
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
 
 
  The results are written in directory \b resultats. Some files are updated at each time-step:
  files energie.dat and temps_reprise.dat give respectively the evolution of energy and the simulation time.
  Files fluide*.vtk and solide*.vtk are written a limited number of times in the span of the simulation. fluide*.vtk and solide*.vtk give respectively the state of the fluid and the position of the solid. They can be read using Paraview.
  File temps.dat gives the cpu cost at the end of the simulation. \n
  It is possible to restart interrupted simulations from recovery files fluide*.vtk and solide*.vtk. It suffices to change recovery flag bool rep = false to bool rep=true in file parametres.h and indicate the recovery point with int numrep.
 
 
  \remark Procedures dealing with the fluid are in files fluide.hpp and
  fluide.cpp. Procedures dealing with the solid are in fildes
  solide.hpp and solide.cpp. Coupling procedures are gathered in files
  couplage.cpp, intersections.hpp and intersections.cpp. FIles
  parametres.hpp and parametres.cpp are dedicated to definitions of
  problem parameters. 
 
 
  \authors Maria Adela Puscas and Laurent Monasse
*/


/*!
 *  \file
 *  \brief Main function. Initialization of the problem and resolution.
 */

#include <iostream>
#include <ctime>
#include "fluide.cpp"
#include "solide.cpp" 
#include "couplage.cpp"
#include "parametres.cpp"
using namespace std;          

/*!\brief Initialization of the problem and resolution:

 - Initialization of the solid and the fluid using respectively functions \a Solide.Init(const char*) and \a Grille.Init().
 - Resolution of the problem:
 - Resolution of the fluid equations using function \a Grille.Solve(const double, double, int).
 - Computation of internal forces using function \a Solide.Forces_internes().
 -  Computation of fluid forces (\a Particule.Ff) and torques (\a Particule.Mf) acting on the solid using function \a Grille.Forces_fluide(Solide&, const double).
 - Update of the position of the solid using function \a Solide.Solve_position(double).
 - Computation of the solid velocity using function \a Solide.Solve_vitesse(double).
 - Intersection of the fluid grid with the solid using function \a Grille.Parois(Solide&, double).
 - Computation of the fluid quantity swept by the solid using function \a Grille.Swap_2d(double, Solide&).
 - Modification of fluid fluxes using function \a Grille.Modif_fnum(double).
 - Conservative mixing of small cut-cells using function \a Grille.Mixage().
 - Filling of ghost cells using function \a Grille.Fill_cel(Solide&).
 - Imposing boundary conditions using function \a Grille.BC().
 
 \return int
 */
int main(){
  char temps_it[]="resultats/temps.dat";
  char temps_reprise[]="resultats/temps_reprise.dat";
  //In case of recovery
  double temps[numrep+1];
  if(rep){
    reprise();
    std::ifstream in(temps_reprise,ios::in);
    if(!in){
      cout <<"Opening of 'temps_reprise.dat' failed" << endl;
    }
    for(int i=0;i<numrep+1;i++){
      in >> temps[i];
    }
    //Energy recovery
    int result = system("cp resultats/energie.dat resultats/energie_reprise.dat");
  }
  std::ifstream in_energie("resultats/energie_reprise.dat",ios::in);
  
  
  //Open output fluxes
  std::ofstream temps_iter(temps_it,ios::out);
  std::ofstream sorties_reprise(temps_reprise,ios::out);
  if(!temps_iter){
    cout <<"Opening of 'temps.dat' failed" << endl;
  }
  if(rep){
    for(int i=0;i<numrep+1;i++){
      sorties_reprise << temps[i] << endl;
    }
  }
  
  
  char energie[]="resultats/energie.dat";
	
//Open output fluxes
  std::ofstream ener(energie,ios::out);
  if(!ener){
    cout <<"Opening of  'energie.dat' failed" << endl;
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
  S.Init("maillage.dat"); //Initialization of solid
  Grille Fluide;
  Fluide.Init();
  Fluide.Parois_particles(S,dt);
  Fluide.BC();
	
  double volume_initial= 0.;
	
  for (int count=0; count<S.size(); count++){
    volume_initial +=S.solide[count].volume();
  }
  	
  int iter=0;	
  clock_t start,end;
  start =clock();

  int kimp = 0; //Output index
  double next_timp = dtimp; //Next output time
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
	
  for (int n=0; (t<T) && n<Nmax; n++){
    user_time_total.start();

    cout<<"iteration="<<n<< " dt="<<dt<<" t="<<t<<endl;
    			
    			
    if(t>next_timp){
      Fluide.Impression(kimp);
      S.Impression(kimp);
      sorties_reprise << t << endl;
      kimp++;
      next_timp += dtimp;
    }
    cout<<"Fluid energy: "<< Fluide.Energie() << " Solid energy:" << S.Energie() <<"  "<<"Fluid mass : "<<"  "<< Fluide.Masse() <<endl;
    ener << t << " " << Fluide.Energie()+S.Energie() << " " << S.Energie() << " " << Fluide.Energie()+S.Energie()-E0 << " " << S.Energie()-E0S <<" "<<Fluide.Masse() - masse <<endl;
    cout<<"Variation Energie: "<< Fluide.Energie() +  S.Energie() - E0<<" Variation Masse : "<< Fluide.Masse() - masse<<endl;
    //Time step
    dt = min(Fluide.pas_temps(t, T),S.pas_temps(t,T));

    Fluide.affiche();
    user_time2.start();

    Fluide.Solve(dt, t, n, S);

    Fluide.affiche();
    temps_flux += CGAL::to_double(user_time2.time());
    user_time2.reset();
    if(explicite){ //Explicit coupling algorithm
      user_time4.start();
      cout <<"Forces_fluide Mass Variation : "<< Fluide.Masse() - masse<<endl;
      Fluide.Forces_fluide(S,dt);
      cout <<"Solve_position Mass Variation : "<< Fluide.Masse() - masse<<endl;
      S.Solve_position(dt);
      temps_explicit += CGAL::to_double(user_time4.time());
      user_time4.reset();
      user_time4.start();
      cout <<"Parois_particles Mass Variation : "<< Fluide.Masse() - masse<<endl;
      Fluide.Parois_particles(S,dt);
      temps_intersections += CGAL::to_double(user_time4.time());
      user_time4.reset();
    }
			
    else{ //Semi-implicit coupling algorithm
      //semi-implicit
      int kmax=20;
      double erreur = 1.;
      user_time5.start();
      Solide Sk = S;
      Solide Sk_prev = Sk;
      int k;
      for(k=0;(erreur>1.e-10) && (k<kmax) ;k++){
	cout <<"Forces_fluide semiimpl Mass Variation : "<< Fluide.Masse() - masse<<endl;
	Fluide.Forces_fluide(Sk,dt);
	Copy_f_m(S,Sk); //This function is important to prevent erasing Sk.F_f et Sk.M_f!!!!!
	Sk_prev = Sk; 
	Sk = S ;
	cout <<"Solve_position semiimpl Mass Variation : "<< Fluide.Masse() - masse<<endl;
	Sk.Solve_position(dt);
	cout <<"Parois_particles semiimpl Mass Variation : "<< Fluide.Masse() - masse<<endl;
	Fluide.Parois_particles(Sk,dt);
	erreur = Error(Sk, Sk_prev);
      }//end semi-implicit fixed point loop on surfaces
      S=Sk;
      temps_semi_implicit += CGAL::to_double(user_time4.time());
      user_time5.reset();
      cout<<"number of semi-implicit iterations: "<<k<<endl;
      nb_iter_implicit += k;
      //semi-implicit	
    }
    user_time3.start();
    cout <<"Forces_internes Mass Variation : "<< Fluide.Masse() - masse<<endl;
    S.Forces_internes();
    temps_solide_f_int += CGAL::to_double(user_time3.time());
    user_time3.reset();
    cout <<"Solve_vitesse Mass Variation : "<< Fluide.Masse() - masse<<endl;
    user_time3.start();
    S.Solve_vitesse(dt);
    temps_solide_vitesse += CGAL::to_double(user_time3.time());
    user_time3.reset();
    user_time.start();
    cout <<"Swap_2d Mass Variation : "<< Fluide.Masse() - masse<<endl;
    Fluide.Swap_2d(dt,S);
    temps_swap += CGAL::to_double(user_time.time());
    user_time.reset();
    cout <<"Modif_fnum Mass Variation : "<< Fluide.Masse() - masse<<endl;
    user_time.start();
    Fluide.Modif_fnum(dt);
    temps_modif_fnum += CGAL::to_double(user_time.time());
    user_time.reset();
    user_time.start();
    bool test_fini = false;
    for(int count=1;!test_fini && count<100;count++){
      cout <<"Mixage_cible2 iteration " << count << " Mass Variation : "<< Fluide.Masse() - masse<<endl;
      test_fini = Fluide.Mixage_cible2(); 
      cout << "iterations of Mixage_cible2=" << count << endl;
    }
    temps_mixage += CGAL::to_double(user_time.time());
    user_time.reset();
    user_time.start();
    cout <<"Fill_cel Mass Variation : "<< Fluide.Masse() - masse<<endl;
    Fluide.Fill_cel(S);
    temps_fill_cel += CGAL::to_double(user_time.time());
    user_time.reset();
    user_time.start();
    cout <<"BC Variation Masse : "<< Fluide.Masse() - masse<<endl;
    Fluide.BC();
    cout <<"BC_couplage Variation Masse : "<< Fluide.Masse() - masse<<endl;
    
    temps_BC += CGAL::to_double(user_time.time());
    user_time.reset();
			
    t+= dt;
    iter++;
    variation_masse += Fluide.Masse() - masse;
    variation_energy += Fluide.Energie()+S.Energie()-E0;
    volume_solide = 0.;
    for (int count=0; count<S.size(); count++){
      volume_solide +=S.solide[count].volume();
    }
    //Fluide.affiche();
    cout<<"volume solid particles "<<volume_solide<<endl;
    variation_volume += volume_solide - volume_initial;

		
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
    cout << "Rest=           " << 100-100*(temps_flux+temps_explicit+temps_intersections+temps_solide_f_int+temps_solide_vitesse+temps_swap+temps_modif_fnum+temps_mixage+temps_fill_cel+temps_BC)/temps_total << "%     " << (temps_total-(temps_flux+temps_explicit+temps_intersections+temps_solide_f_int+temps_solide_vitesse+temps_swap+temps_modif_fnum+temps_mixage+temps_fill_cel+temps_BC))/(n+1.) << endl;
    cout << "########################################" << endl;
		
  }
  end=clock();
  Fluide.Impression(kimp);
  S.Impression(kimp);
	
  temps_iter<< "Final time  "<< t<<endl;
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
