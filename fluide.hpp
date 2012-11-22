#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <sstream>
#include <cassert>
#include "parametres.hpp"
#include "solide.hpp"



#ifndef FLUIDE_HPP
#define FLUIDE_HPP

using namespace std;


//Definition de la classe Cellule 
class Cellule {

public :

 //Constructeur 
  Cellule();
  Cellule(double x, double y, double z); // surcharge du constructeur
  Cellule(double x, double y, double z, double dx, double dy, double dz); // surcharge du constructeur
  
  Cellule & operator=(const Cellule &cell); // operateur = surcharge pour l'affectation

  //Destructeur
  ~Cellule();

  //Fonction testant si un point est dans la cellule
  bool is_in_cell(double x,double y, double z);

  void Affiche ();  //fonction auxilaire utile pour les test :)

protected:
  //test 9 nov
	double delta_w;
	//fin test 9 nov
  double x;        //position du centre de la cellule 
  double y;
  double z; 

  double dx;       //Taille de la cellule
  double dy;
  double dz;

  double rho;       //densite dans la cellule
  double rho1;      //densite dans la cellule en t-dt

  double u;        //vitesse dans la cellule selon x
  double v;        //vitesse dans la cellule selon y
  double w;        //vitesse dans la cellule selon z

  double p;        //pression dans la cellule
  double p1;        //pression dans la cellule en t-dt

  double pdtx;     //Pression efficace selon la direction x pendant le pas de temps dt
  double pdty;     //Pression efficace selon la direction y pendant le pas de temps dt
  double pdtz;     //Pression efficace selon la direction z pendant le pas de temps dt

  double impx;      //impulsion selon x
  double impy;      //impulsion selon y 
  double impz;      //impulsion selon z

  double impx0;      //impulsion selon x avant calcul des flux
  double impy0;      //impulsion selon y avant calcul des flux
  double impz0;      //impulsion selon y avant calcul des flux

  double rhoE;      //densite d'energie dans la cellule
  double rhoE0;     //densite d'energie dans la cellule avant calcul des flux
  double rho0;      //densite dans la cellule avant calcul des flux



  double lambda[5];     //Valeurs propres du systeme
  double rp[5];         //Variables pour le limiteur de flux TVD decentre a gauche 
  double rm[5];        //Variables pour le limiteur de flux TVD decentre a droite 

  double Mrho;      //Terme de mixage de densite pour les petites cellules
  double Mimpx;     //Terme de mixage d'impulsion pour les petites cellules
  double Mimpy;     //Terme de mixage d'impulsion pour les petites cellules
  double Mimpz;     //Terme de mixage d'impulsion pour les petites cellules
  double MrhoE;     //Terme de mixage d'energie pour les petites cellules

  double cells;      //Volume des cellules melangees avec la cellule si elle est cible
  double alpha;      //proportion de la cellule occupee par du solide a l'instant t
  double alpha0;    //proportion de la cellule occupee par du solide a l'instant t-dt

  double kappai;    //Taux d'occupation des faces par du solide
  double kappaj;
  double kappak;

  double kappai0;    //Taux d'occupation des faces par du solide au temps t-dt
  double kappaj0;
  double kappak0;

  int proche;         //Indicateur qui vaut 0 loin de la paroi, 1 pr�s de la paroi
  int proche1;       //Indicateur de proche au temps t-dt

  
  double fluxi[5];     //Flux a droite de la cellule
  double fluxj[5];     //Flux en haut de la cellule
  double fluxk[5];     //Flux en bas de la cellule
  
  double flux_modif[5]; 
	double phi_x;   //necessaire pour le calcul du flux à la parois
  double phi_y;
  double phi_z;     
  double phi_v;  
	
  double xi;            //Position du centre de la face effective 
  double yj;
  double zk;

    
   
  double S;                  //Entropie physique
  double ve[5];             //Vecteur des variables entropiques

  double fex;               //Flux d'entropie suivant x
  double fey;               //Flux d'entropie suivant y
  double fez;               //Flux d'entropie suivant z
    

  double Qci[5];            //Correction d'entropie en i
  double Qcj[5];            //Correction d'entropie en j
  double Qck[5];            //Correction d'entropie en k

  double dtfxi[5];          //flux suivant x en i fois le pas de temps dt
  double dtfyj[5];          //flux suivant y en j fois le pas de temps dt
  double dtfzk[5];          //flux suivant z en k fois le pas de temps dt

  double delw[5];           //delta V en i+1/2
  double delwnu[5];         // facteur devant le limitateur en i+1/2

  double cf2[5];            //Corrections pour l'ordre superieur 
  double cf3[5]; 
  double cf4[5]; 
  double cf5[5]; 
  double cf6[5]; 
  double cf7[5]; 
  double cf8[5]; 
  double cf9[5]; 
  double cf10[5]; 
  double cf11[5]; 


  double psic0[5];           //Corrections centrees pour l'ordre superieur
  double psic1[5]; 
  double psic2[5]; 
  double psic3[5]; 
  double psic4[5]; 

  double psid0[5];            //Corrections decentrees pour l'ordre superieur
  double psid1[5]; 
  double psid2[5]; 
  double psid3[5]; 
  double psid4[5]; 

  double vpr[5][5];            //Matrice des vecteurs propres du systeme


  double psic0r[5];           //Corrections centrees pour l'ordre superieur mises dans la base des vecteurs propres 
  double psic1r[5]; 
  double psic2r[5]; 
  double psic3r[5]; 
  double psic4r[5]; 

  double psid0r[5];            //Corrections decentrees pour l'ordre superieur mises dans la base des vecteurs propres 
  double psid1r[5]; 
  double psid2r[5]; 
  double psid3r[5]; 
  double psid4r[5]; 

  double psid[5]; 
  double am[5];                 //Variables de mesure de la monotonicite
  double am1[5];  
 
  int ordre;                    //ordre du schema pour la cellule 
  double co[11];               //stockage des coefficients d'ordre 
  
  friend class Grille;         //la class Cellule est visible seulement a travers la class Grille

};

class Grille
{

 public:

  //Constructeur 
  Grille();
  //surcharge du Constructeur 
 Grille(int Nx0,int Ny0, int Nz0, double dx0, double x0, double dy0, double y0, double dz0, double z0);

  //Destructeur 
  ~Grille();

  void affiche(); //fonction auxilaire utile pour les test :)

  void affiche(string r); //fonction auxilaire utile pour les test :)

  //Acces a une cellule i,j, k
   Cellule cellule(int i, int j, int k);
  //Acces a la cellule contenant le point (x,y,z)
   Cellule in_cell(Point_3 p);
	//Acces a la cellule contenant le point (x,y,z)
	void in_cell(Point_3 p, int &i, int& j, int& k, bool& interieur);
  //Sous-programme de definition des conditions initiales 
  void init();

  //Sous-programme de calcul du pas de temps 
  double pas_temps(double t, double T);
  
  //Sous-programme pour mettre des conditions aux limites 
  void BC();

   //Calcul des quantites conservatives totales 
   double Masse();
   double Impulsionx();
   double Impulsiony();
   double Impulsionz();
   double Energie();
 
  //Sous-programme d'impression des resultats 
  void impression(int n);
	   
  //Resolution des equations pour le fluide
 // void solve_fluid(double dt, double t);

  //Resolution des equations pour le fluide dans la direction x
  void solve_fluidx(const double dt);
  //Resolution des equations pour le fluide dans la direction y
  void solve_fluidy(const double dt);
  //Resolution des equations pour le fluide dans la direction z
  void solve_fluidz(const double dt);
	//Fonction qui melange les cellules a densite ou pression negative 
	void melange(const double dt);
  //Calcul du flux dans la direction x entre les cellules 
  void fnumx( const double sigma, double t);
  //Calcul du flux dans la direction y entre les cellules 
  void fnumy(const double sigma, double t);
  //Calcul du flux dans la direction z entre les cellules 
  void fnumz( const double sigma, double t);


  //Correction d'entropie selon x
  void corentx(double sigma);
  //Correction d'entropie selon y
  void corenty(double sigma);
  //Correction d'entropie selon z
  void corentz(double sigma);
   
  void Solve(const double dt, double t, int n);
  
  void parois(Solide& S);  // Mise a jour de Kappai,kappaj,kappak et alpha   
  void modif_fnum(const double dt);  //Modification du flux
  void mixage(); //Procedure de mixage pour le cellules avec c.alpha>0.5
  void fill_cel(Solide& S); // Remplissage de cellules fantomes (c.alpha = 1.)
	//void fill_cel_old(Solide& S); // Remplissage de cellules fantomes (c.alpha = 1.) old version
	void swap_modification_flux(Triangles& T3d_prev, Triangles& T3d_n, const double dt);
	void cells_intersection_face(int& in,int& jn,int& kn,int& in1,int& jn1,int& kn1, std::vector<Bbox>& box_cells, std::vector<Cellule>& Cells);
	void swap(const double dt, Solide& S,int& n, int &n1, int& m);
private :

  double x;          //position de l'origine
  double y;
  double z;
  double dx;          //pas d'espace
  double dy;
  double dz;
	vector< vector< vector<Cellule > > > grille;
	//Cellule grille[Nx+2*marge][Ny+2*marge][Nz+2*marge];   //tableau des cellules 
    
};

#endif
