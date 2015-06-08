/*!
 *  \file fluide.hpp
 *  \brief D&eacute;finition des classes Cellule et Grille n&eacute;cessaires pour la r&eacute;solution du fluide.
 * Les membres sp&eacute;cifiques au couplage sont pr&eacute;c&egrave;des d'un "warning".
 */

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


//! D&eacute;finition de la classe Cellule 
class Cellule {

public :

 // Constructeur 
  Cellule();
  Cellule(double x, double y, double z); // Surcharge du constructeur
  Cellule(double x, double y, double z, double dx, double dy, double dz); // Surcharge du constructeur
  
  Cellule & operator=(const Cellule &cell); // operateur = Surcharge pour l'affectation

  // Destructeur
  ~Cellule();

  //Fonction testant si un point est dans la cellule
  bool is_in_cell(double x,double y, double z);

  void Affiche ();  // Fonction auxilaire utile pour les tests 
	
// 	double getpdtx()const {return pdtx;} // return pdtx dans la cellule
// 	double getpdty()const {return pdty;} // return pdty dans la cellule
// 	double getpdtz()const {return pdtz;} // return pdtz dans la cellule
  
  //! (x,y,z) Position du centre de la cellule.   
  double x;       
  double y;
  double z; 
	
	int i;
	int j;
	int k;

  //! (dx, dy,dz) Taille de la cellule.
  double dx;       
  double dy;
  double dz;

	double rho;       //!< Densit&eacute; dans la cellule.
	double rho1;      //!< Densit&eacute; dans la cellule en t-dt.

  double u;        //!< Vitesse dans la cellule selon x.
  double v;        //!< Vitesse dans la cellule selon y.
  double w;        //!< Vitesse dans la cellule selon z.

double p;         //!< Pression dans la cellule.
double p1;        //!< Pression dans la cellule en t-dt.
    

double impx;      //!< Impulsion selon x.
double impy;      //!< Impulsion selon y.
double impz;      //!< Impulsion selon z.

double impx0;      //!< Impulsion selon x avant calcul des flux.
double impy0;      //!< Impulsion selon y avant calcul des flux.
double impz0;      //!< Impulsion selon z avant calcul des flux.

double rhoE;      //!< Densit&eacute; d'&eacute;nergie dans la cellule.
double rhoE0;     //!< Densit&eacute; d'&eacute;nergie dans la cellule avant calcul des flux.
double rho0;      //!< Densit&eacute; dans la cellule avant calcul des flux.



double lambda[5];     //!< Valeurs propres du syst&egrave;me.
double rp[5];         //!< Variables pour le limiteur de flux TVD d&eacute;centre &agrave; gauche. 
double rm[5];        //!< Variables pour le limiteur de flux TVD d&eacute;centre &agrave; droite. 
   
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
    double pdtx;     //!< Pression efficace selon x pendant le pas de temps dt. 
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
    double pdty;     //!< Pression efficace selon y pendant le pas de temps dt. 
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
    double pdtz;     //!< Pression efficace selon z pendant le pas de temps dt. 
    
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */

		double Mrho;      //!< Terme de mixage de densit&eacute; pour les petites cellules.  
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double Mimpx;     //!< Terme de mixage d'impulsion pour les petites cellules.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double Mimpy;     //!< Terme de mixage d'impulsion pour les petites cellules.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double Mimpz;     //!< Terme de mixage d'impulsion pour les petites cellules.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double MrhoE;     //!< Terme de mixage d'&eacute;nergie pour les petites cellules.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double cells;      //!< Volume des cellules m&eacute;lang&eacute;es avec la cellule si elle est cible.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double alpha;      //!< Proportion de la cellule occup&eacute;e par du solide &agrave; l'instant t.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double alpha0;    //!< Proportion de la cellule occup&eacute;e par du solide &agrave; l'instant t-dt.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double kappai;    //!< Taux d'occupation des faces par du solide selon x.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double kappaj; //!< Taux d'occupation des faces par du solide selon y.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double kappak; //!< Taux d'occupation des faces par du solide selon z.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double kappai0;    //!< Taux d'occupation des faces par du solide selon x au temps t-dt.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double kappaj0; //!< Taux d'occupation des faces par du solide selon y au temps t-dt.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double kappak0; //!< Taux d'occupation des faces par du solide selon z au temps t-dt.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		int proche;         //!< Indicateur qui vaut 0 loin de la paroi, 1 pr&egrave;s de la paroi.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		int proche1;       //!< Indicateur qui vaut 0 loin de la paroi, 1 pr&egrave;s de la paroi au temps t-dt.
		
		/*! 
		* \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
		*/
    bool vide;
  
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double flux_modif[5]; //!< Modification flux pour les cellules coup&eacute;es.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double delta_w[5];    //!< Quantit&eacute; balay&eacute;e.
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double phi_x;        //!< Flux &agrave; la parois. 
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double phi_y; //!< Flux &agrave; la parois .
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double phi_z; //!< Flux &agrave; la parois. 
    /*! 
     * \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
     */
		double phi_v;  //!< Flux &agrave; la parois. 
		
		double cible_alpha;      //!< Terme de mixage de densit&eacute; pour les petites cellules.  
		/*! 
		* \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
		*/
		double cible_rho;      //!< Terme de mixage de densit&eacute; pour les petites cellules.  
		/*! 
		* \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
		*/
		double cible_impx;     //!< Terme de mixage d'impulsion pour les petites cellules.
		/*! 
		* \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
		*/
		double cible_impy;     //!< Terme de mixage d'impulsion pour les petites cellules.
		/*! 
		* \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
		*/
		double cible_impz;     //!< Terme de mixage d'impulsion pour les petites cellules.
		/*! 
		* \warning  <b>  Param&egrave;tre sp&eacute;cifique  au  couplage! </b>
		*/
		double cible_rhoE;     //!< Terme de mixage d'&eacute;nergie pour les petites cellules.
		
		double cible_i;      //!< Index de la cellule cible 
		double cible_j;
		double cible_k;
	
	double xi;      //!< Position du centre de la face effective. 
  double yj;
  double zk;

	double fluxi[5];     //!< Flux &agrave; droite de la cellule.
	double fluxj[5];     //!< Flux en haut de la cellule.
	double fluxk[5];     //!< Flux en bas de la cellule.
   
   double S;                  //!< Entropie physique.
   double ve[5];             //!< Vecteur des variables entropiques.

double fex;               //!< Flux d'entropie suivant x.
double fey;               //!< Flux d'entropie suivant y.
double fez;               //!< Flux d'entropie suivant z.
    

double Qci[5];            //!< Correction d'entropie en i.
double Qcj[5];            //!< Correction d'entropie en j.
double Qck[5];            //!< Correction d'entropie en k.

double dtfxi[5];          //!< Flux suivant x en i fois le pas de temps dt.
double dtfyj[5];          //!< Flux suivant y en j fois le pas de temps dt.
double dtfzk[5];          //!< Flux suivant z en k fois le pas de temps dt.

double delw[5];           //!< Delta V en i+1/2.
double delwnu[5];         //!< Facteur devant le limitateur en i+1/2.

double cf2[5];            //!< Corrections pour l'ordre sup&eacute;rieur. 
  double cf3[5]; 
  double cf4[5]; 
  double cf5[5]; 
  double cf6[5]; 
  double cf7[5]; 
  double cf8[5]; 
  double cf9[5]; 
  double cf10[5]; 
  double cf11[5]; 


	double psic0[5];           //!< Corrections centr&eacute;es pour l'ordre superieur.
  double psic1[5]; 
  double psic2[5]; 
  double psic3[5]; 
  double psic4[5]; 

	double psid0[5];            //!< Corrections decentr&eacute;es pour l'ordre superieur.
  double psid1[5]; 
  double psid2[5]; 
  double psid3[5]; 
  double psid4[5]; 

	double vpr[5][5];            //!< Matrice des vecteurs propres du syst&egrave;me.


double psic0r[5];           //!< Corrections centr&eacute;es pour l'ordre sup&eacute;rieur mises dans la base de vecteurs propres. 
  double psic1r[5]; 
  double psic2r[5]; 
  double psic3r[5]; 
  double psic4r[5]; 

	double psid0r[5];            //!< Corrections decentr&eacute;es pour l'ordre sup&eacute;rieur mises dans la base de vecteurs propres. 
  double psid1r[5]; 
  double psid2r[5]; 
  double psid3r[5]; 
  double psid4r[5]; 

  double psid[5]; 
	double am[5];                 //!< Variables de mesure de la monotonicit&eacute;.
  double am1[5];  
 
	int ordre;                    //!< Ordre du sch&eacute;ma pour la cellule.
	double co[11];               //!< Stockage des coefficients d'ordre. 
  
  friend class Grille;         //!< La classe Cellule est visible uniquement &agrave; travers la classe Grille.

};

//! D&eacute;finition de la classe Grille
class Grille
{

 public:

  // Constructeur 
  Grille();
  // surcharge du Constructeur 
 Grille(int Nx0,int Ny0, int Nz0, double dx0, double x0, double dy0, double y0, double dz0, double z0);

  // Destructeur 
  ~Grille();

  void affiche(); // Fonction auxilaire utile pour les tests 

  void affiche(string r); // Fonction auxilaire utile pour les tests

  // Acces a une cellule i,j, k
   Cellule cellule(int i, int j, int k);
  // Acces a la cellule contenant le point p
   Cellule in_cell(const Point_3& p);
	// Acces a la cellule contenant le point p
   void in_cell(const Point_3& p, int &i, int& j, int& k, bool& interieur);
  // Sous-programme de definition des conditions initiales 
  void Init();

  // Sous-programme de calcul du pas de temps 
  double pas_temps(double t, double T);
  
  // Sous-programme pour mettre des conditions aux limites 
  void BC();
  void BC_couplage(double tab[marge][3], bool couplage1d);
	void BC_couplage_1d(vector< vector < double> > &tab_1d);


   // Calcul des quantites conservatives totales 
   double Masse();
   double Impulsionx();
   double Impulsiony();
   double Impulsionz();
   double Energie();
 
  // Sous-programme d'impression des resultats 
  void Impression(int n);
	   
  // Resolution des equations pour le fluide dans la direction x
  void solve_fluidx(const double dt);
  // Resolution des equations pour le fluide dans la direction y
  void solve_fluidy(const double dt);
  // Resolution des equations pour le fluide dans la direction z
  void solve_fluidz(const double dt);
  // Fonction qui melange les cellules a densite ou pression negative 
  void melange(const double dt);
  // Calcul du flux dans la direction x entre les cellules 
  void fnumx( const double sigma, double t, double tab[marge][3], bool couplage1d);
  // Calcul du flux dans la direction y entre les cellules 
  void fnumy(const double sigma, double t);
  // Calcul du flux dans la direction z entre les cellules 
  void fnumz( const double sigma, double t);


  //Correction d'entropie selon x
  void corentx(double sigma);
  //Correction d'entropie selon y
  void corenty(double sigma);
  //Correction d'entropie selon z
  void corentz(double sigma);
  void Solve(const double dt, double t, int n, double tab[marge][3], bool couplage1d, Solide& S); 
    
  void Forces_fluide(Solide& S, const double dt); // Calcul des Forces fluides et Moments fluides exerces sur le solide	
  void Modif_fnum(const double dt);  // Modification du flux
  void Mixage(); // Procedure de mixage pour le cellules avec c.alpha>0.5
  void Fill_cel(Solide& S); // Remplissage de cellules fantomes (c.alpha = 1.)
  void swap_face(const Triangles& T3d_prev, const Triangles& T3d_n, const double dt,  Particule & P, double & volume_test); // Sous-programme de calcul des quantitees balayees
  void swap_face_inexact(const Triangle_3& Tr_prev, const Triangle_3& Tr, const Triangles& T3d_prev, const Triangles& T3d_n, const double dt,  Particule & P, double & volume_test); // Sous-programme de calcul inexact des quantitees balayees
  void cells_intersection_face(int& in,int& jn,int& kn,int& in1,int& jn1,int& kn1, std::vector<Bbox>& box_cells, std::vector<Cellule>& Cells); 
  void Swap_2d(const double dt, Solide& S); // Sous-programme de calcul des quantitees balayees
  void Swap_3d(const double dt, Solide& S); // Sous-programme de calcul des quantitees balayees
	
	Cellule voisin_fluide(const Cellule &c,  bool &target);
	Cellule voisin_mixt(const Cellule &c,  bool &target);
	Cellule voisin(const Cellule &c);
	Cellule cible(const Cellule &c, std::vector< std::vector<int> > & tab_cible );
	void Mixage_cible();
  bool Mixage_cible2();
	void Parois_particles(Solide& S,double dt);
  std::vector<Point_3> intersection(Triangle_3 t1, Triangle_3 t2);
  //void Parois_tetra(Solide& S,double dt);
	void parois_cellule_vide(Solide& S);
//private :

	double x;          //!< Position de l'origine.
  double y;
  double z;
	double dx;          //!< Pas d'espace.
  double dy;
  double dz;
	vector< vector< vector<Cellule > > > grille; //!< Maillage fluide.
	//Cellule grille[Nx+2*marge][Ny+2*marge][Nz+2*marge];   //tableau des cellules 

};

#endif
