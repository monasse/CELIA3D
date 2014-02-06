#include <cstdlib> 
#include <iostream> 
#include <stdio.h> 
#include <fstream> 
#include <cmath>


using namespace std;


const int marge_1d = 6;   
const int N=471;
const double Dx=  0.008333333;
//Conditions aux limites
const string CL1D = "outflow";        

inline double sign1d(const double x)
{
	return (x < 0.) ? -1. : 1. ;
}
double rho_1D( double x){
	
	double rho = 1.;
	if(x<0.75 ){rho = 1.62865;} 
	else {rho = 1.1845467;}
	
	return rho;
}	

double u_1D( double x){
	
	double u=0.;

	
	return u;
}

double p_1D( double x){
	
	double p = 101325.; 
	if(x<0.75){p = 158235.;} 
	else {p = 101325.;}
	
	return p;
} 



//Definition de la classe Cellule 
class Cellule1D 
{ 
	public : 
		
		double x;        //position du centre de la cellule 
		double rho;      //densit� dans la cellule 
		double u;        //vitesse dans la cellule 
		double p;        //pression dans la cellule 
		double imp;      //impulsion 
		double rhoE;     //densit� d'�nergie dans la cellule 

		double lambda[3];//Valeurs propres du syst�me	    
		double rp[3];    //Variables pour le limiteur de flux TVD d�centr� �� gauche 
		double rm[3];    //Variables pour le limiteur de flux TVD d�centr� �� droite 
		double fluxd[3]; //Flux �� droite de la cellule 
		double delw[3];  //delta V en i+1/2 
		double delwnu[3];//facteur devant le limiteur en i+1/2 
		double cf2[3];   //Corrections pour l'ordre sup�rieur 
		double cf3[3]; 
		double cf4[3]; 
		double cf5[3]; 
		double cf6[3]; 
		double cf7[3]; 
		double cf8[3]; 
		double cf9[3]; 
		double cf10[3]; 
		double cf11[3]; 
		double psic0[3]; //Corrections centr�es pour l'ordre sup�rieur 
		double psic1[3]; 
		double psic2[3]; 
		double psic3[3]; 
		double psic4[3]; 
		double psid0[3]; //Corrections d�centr�es pour l'ordre sup�rieur 
		double psid1[3]; 
		double psid2[3]; 
		double psid3[3]; 
		double psid4[3]; 
		double vpr[3][3]; //Matrice des vecteurs propres du syst�me 
		double psic0r[3]; //Corrections centr�es pour l'ordre sup�rieur mises dans la base des vecteurs propres 
		double psic1r[3]; 
		double psic2r[3]; 
		double psic3r[3]; 
		double psic4r[3]; 
		double psid0r[3]; //Corrections d�centr�es pour l'ordre sup�rieur mises dans la base des vecteurs propres 
		double psid1r[3]; 
		double psid2r[3]; 
		double psid3r[3]; 
		double psid4r[3]; 
		double psid[3]; 
		double am[3];    //Variables de mesure de la monotonicit� 
		double am1[3];  
		
		int ordre;       //ordre du sch�ma pour la cellule 
		double co[11];   //stockage des coefficients d'ordre 
		double S;
		double ve[3];
		double fe;
		double Qcd[3];
	

		
		//Constructeur 
		Cellule1D (){ 
			x=1.;         
			rho = 1.;     
			u=0.;         
			p=1.;         
			imp=rho*u;    
			rhoE=rho*u*u/2.+p/(gam-1.);  
       
			lambda[0] = -1.; 
			lambda[1] = 0.; 
			lambda[2] = 1.;	    
			rp[0] = 1.; 
			rp[1] = 1.; 
			rp[2] = 1.; 
			rm[0] = 1.; 
			rm[1] = 1.; 
			rm[2] = 1.; 
			fluxd[0] = 0.; 
			fluxd[1] = 0.; 
			fluxd[2] = 0.; 
			delw[0] = 1.; 
			delw[1] = 1.; 
			delw[2] = 1.;   
			delwnu[0] = 1.; 
			delwnu[1] = 1.; 
			delwnu[2] = 1.; 
			for(int l=0;l<3;l++){ 
				cf2[l] = 0.; 
				cf3[l] = 0.; 
				cf4[l] = 0.; 
				cf5[l] = 0.; 
				cf6[l] = 0.; 
				cf7[l] = 0.; 
				cf8[l] = 0.; 
				cf9[l] = 0.; 
				cf10[l] = 0.; 
				cf11[l] = 0.; 
				psic0[l] = 0.; 
				psic1[l] = 0.; 
				psic2[l] = 0.; 
				psic3[l] = 0.; 
				psic4[l] = 0.; 
				psid0[l] = 0.; 
				psid1[l] = 0.; 
				psid2[l] = 0.; 
				psid3[l] = 0.; 
				psid4[l] = 0.; 
				for(int m=0;m<3;m++){ 
					vpr[l][m] = 0.; 
				}      
				psic0r[l] = 0.; 
				psic1r[l] = 0.; 
				psic2r[l] = 0.; 
				psic3r[l] = 0.; 
				psic4r[l] = 0.; 
				psid0r[l] = 0.;  
				psid1r[l] = 0.; 
				psid2r[l] = 0.; 
				psid3r[l] = 0.; 
				psid4r[l] = 0.; 
				psid[l] = 0.; 
				am[l] = 0.; 
				am1[l] = 0.; 
			} 
			ordre = ordremax; 
			for(int l=0; l<ordremax;l++){ 
				co[l] = 1.; 
			} 
			
		} 
		
		Cellule1D (double x0){ 
			x = x0; 
			rho = 1.;     
			u=0.;         
			p=1.;         
			imp=rho*u;    
			rhoE=rho*u*u/2.+p/(gam-1.);     
   
			lambda[0] = -1.; 
			lambda[1] = 0.; 
			lambda[2] = 1.;	    
			rp[0] = 1.; 
			rp[1] = 1.; 
			rp[2] = 1.; 
			rm[0] = 1.; 
			rm[1] = 1.; 
			rm[2] = 1.; 
			fluxd[0] = 0.; 
			fluxd[1] = 0.; 
			fluxd[2] = 0.; 
			delw[0] = 1.; 
			delw[1] = 1.; 
			delw[2] = 1.;   
			delwnu[0] = 1.; 
			delwnu[1] = 1.; 
			delwnu[2] = 1.; 
			for(int l=0;l<3;l++){ 
				cf2[l] = 0.; 
				cf3[l] = 0.; 
				cf4[l] = 0.; 
				cf5[l] = 0.; 
				cf6[l] = 0.; 
				cf7[l] = 0.; 
				cf8[l] = 0.; 
				cf9[l] = 0.; 
				cf10[l] = 0.; 
				cf11[l] = 0.; 
				psic0[l] = 0.; 
				psic1[l] = 0.; 
				psic2[l] = 0.; 
				psic3[l] = 0.; 
				psic4[l] = 0.; 
				psid0[l] = 0.; 
				psid1[l] = 0.; 
				psid2[l] = 0.; 
				psid3[l] = 0.; 
				psid4[l] = 0.; 
				for(int m=0;m<3;m++){ 
					vpr[l][m] = 0.; 
				}      
				psic0r[l] = 0.; 
				psic1r[l] = 0.; 
				psic2r[l] = 0.; 
				psic3r[l] = 0.; 
				psic4r[l] = 0.; 
				psid0r[l] = 0.;  
				psid1r[l] = 0.; 
				psid2r[l] = 0.; 
				psid3r[l] = 0.; 
				psid4r[l] = 0.; 
				psid[l] = 0.; 
				am[l] = 0.; 
				am1[l] = 0.; 
			} 
			ordre = ordremax; 
			for(int l=0; l<ordremax;l++){ 
				co[l] = 1.; 
			} 

		} 
			
		//Destructeur 
		~Cellule1D (){ 
		}           
}; 
//D�finition du noyau d'interpolation B-spline 
double M4(double s){ 
	double M = 0.; 
	if((s>=0.) && (s<=1.)){ 
		M = 1.-5.*s*s/2.+3.*s*s*s/2.; 
	} 
	else if((s>1.) && (s<=2.)){ 
		M = 1./2.*(2.-s)*(2.-s)*(1.-s); 
	} 
	return M; 
} 


//D�finition de la classe Grille 
class Grille1D 
{ 
	public:
		double x0;          //position de l'origine 
		double dx;          //pas d'espace 
		Cellule1D grille[N+2*marge];   //tableau pointant vers les cellules en dynamique 
		
		//Constructeur 
		Grille1D (){ 
			dx = Dx; 
			x0 = 0.; 
			for(int i=0;i<N+2*marge;i++){ 
				Cellule1D  c(x0+dx/2.+(i-marge)*dx); 
				grille[i] = c; 
			} 
		} 
		
		Grille1D (int N,double dx1, double x1){ 
			x0 = x1;
			dx = dx1; 
			for(int i=0;i<N+2*marge;i++){ 
				Cellule1D  c(x1+dx1/2.+(i-marge)*dx1); 
				grille[i] = c;
			} 
		} 
		
		//Destructeur 
		~Grille1D(){ 
		} 
		

 		Cellule1D cellule1D(int i){ 
			return grille[i]; 
		} 
		
		void init(){
		  if(rep){
		    std::ostringstream oss;
		    oss << "resultats/xt" << numrep << ".dat";
		    string s = oss.str();
		    const char* nom = s.c_str();
		    std::ifstream init(nom,std::ios::in);
		    string dump;
		    double x;
		    for(int i=marge;i<N+marge;i++){
		      Cellule1D c = grille[i];
		      init >> x >> c.rho >> c.p >> c.u;
		      c.imp = c.rho*c.u; 
		      c.rhoE= c.rho*c.u*c.u/2.+c.p/(gam-1.);
		      grille[i]=c;
		    }
		    cond_lim();
		  } else {
		    for(int i=0;i<N+2*marge;i++){ 
		      Cellule1D c = grille[i]; 
		      c.rho = rho_1D(c.x); 
		      c.u   = u_1D(c.x); 
		      c.p   = p_1D(c.x); 
		      c.imp = c.rho*c.u; 
		      c.rhoE= c.rho*c.u*c.u/2.+c.p/(gam-1.); 
		      grille[i] = c; 
		    }
		  }
		}
		
		void affiche(){ 
			cout<<"pb 1d"<<endl;
			int cellules=0;
			for(int i=marge;i<N+marge;i++){ 
				Cellule1D c = grille[i]; 
				cellules +=1;
				if(c.p >2){
					cout<<"c.x "<<c.x<<" c.p "<<c.p<<endl;
				}
			} 
			cout<<"nb celllules "<<cellules<<endl;
		}
		
		//Sous-programme de calcul du pas de temps 
		double pas_temps(double t, double T){ 
			double dt = 10000.;
			for(int i=marge;i<N+marge;i++){ 
				  Cellule1D c = grille[i]; 
					double c2 = gam*c.p/c.rho;
					double dt1 = 0.5*dx/(sqrt(c2)+abs(c.u));
					dt = min(dt,dt1); 
			} 
			dt = min(dt,T-t); 
			return dt; 
		}
		void cond_lim_couplage(vector< vector < double> > &tab){
			int j=0;
			for(int i=N+marge;i<N+2*marge;i++){
				Cellule1D c = grille[i];
					c.rho = tab[i-N-marge][0] ;
					c.imp = tab[i-N-marge][1];
					c.rhoE = tab[i-N-marge][2];
					c.u = c.imp/c.rho; 
					c.p = (gam-1.)*(c.rhoE-c.rho*c.u*c.u/2.); 
					grille[i] = c;
					//j++;
			}
		}
		void cond_lim(){ 
			//CL1D �� gauche 
			for(int i=0;i<marge;i++){
				Cellule1D c = grille[i];
				Cellule1D cm = grille[2*marge-i-1]; //Cellule miroir
				Cellule1D cp = grille[N+i];         //Cellule p�riodique
				//if(CL1D == "miroir"){
					// Conditions miroir
					c.rho = cm.rho;
					c.u   = -cm.u;
					c.p   = cm.p;
					c.imp = -cm.imp;
					c.rhoE= cm.rhoE;
					/*}
				else if(CL1D=="periodique"){
					//Conditions p�riodiques
					c.rho = cp.rho;
					c.u   = cp.u;
					c.p   = cp.p;
					c.imp = cp.imp;
					c.rhoE= cp.rhoE;
				}
				else if(CL1D=="outflow"){
					//Conditions  transmisibles
					c.rho = cm.rho;
					c.u   = cm.u;
					c.p   = cm.p;
					c.imp = cm.imp;
					c.rhoE= cm.rhoE;
					}*/
				grille[i] = c;
			}
			//CL1D  droite
			for(int i=N+marge;i<N+2*marge;i++){
				Cellule1D c = grille[i];
				Cellule1D cm = grille[2*N+2*marge-i-1]; //Cellule miroir
				Cellule1D cp = grille[i-N];              //Cellule periodique
				if(CL1D == "miroir"){
					//Conditions miroir
					c.rho = cm.rho;
					c.u   = -cm.u;
					c.p   = cm.p;
					c.imp = -cm.imp;
					c.rhoE= cm.rhoE;
				}
				else if(CL1D=="periodique"){
					//Conditions periodiques
					c.rho = cp.rho;
					c.u   = cp.u;
					c.p   = cp.p;
					c.imp = cp.imp;
					c.rhoE= cp.rhoE;
				}
				else if(CL1D=="outflow"){
					//Conditions  transmisibles
					c.rho = cm.rho;
					c.u   = cm.u;
					c.p   = cm.p;
					c.imp = cm.imp;
					c.rhoE= cm.rhoE;
				}
				grille[i] = c; 
			}	
			
		}
		

		void solve_fluid( double dt, double t){ 
			double sigma = dt/dx; 
			

			fnum1D(sigma,t);
			for(int i=marge;i<N+marge;i++){ 
				Cellule1D c = grille[i]; 
				Cellule1D cg = grille[i-1]; //Cellule � gauche 
				double dw1 = -sigma*(c.fluxd[0]-cg.fluxd[0]); 
				double dw2 = -sigma*(c.fluxd[1]-cg.fluxd[1]); 
				double dw3 = -sigma*(c.fluxd[2]-cg.fluxd[2]); 
				c.rho += dw1; 
				c.imp += dw2; 
				c.rhoE += dw3; 
			  c.u = c.imp/c.rho; 
			  c.p = (gam-1)*(c.rhoE-c.rho*c.u*c.u/2.); 
						
				grille[i] = c; 
		}
		cond_lim(); 
}
void fnum1D( double sigma, double t){ 
	
	
	//Initialisation du flux au flux centr� 
	for(int i=0; i<N+2*marge-1; i++){ 
		Cellule1D c = grille[i];      //Cellule de r�f�rence pour le calcul du flux � droite 
		Cellule1D cd = grille[i+1];   //Cellule de droite 
		Cellule1D cg = grille[i-1];   //Cellule de gauche

			//Calcul d'indicateurs de l'ordre 
			for(int l=0;l<c.ordre;l++){ 
				c.co[l]=1.; 
			} 
			for(int l=c.ordre;l<11;l++){ 
				c.co[l]=0.; 
			} 
			c.fluxd[0] = (c.imp+cd.imp)/2.; 
			c.fluxd[1] = (c.rho*c.u*c.u+c.p+cd.rho*cd.u*cd.u+cd.p)/2.; 
			c.fluxd[2] = ((c.rhoE+c.p)*c.u+(cd.rhoE+cd.p)*cd.u)/2.; 


		
		grille[i] = c; 
		
	} 
	
	//Boucle sur les cellules : calcul des variables pr�liminaires au calcul des flux (variables de limitation) 
	for(int i=0; i<N+2*marge-1; i++){ 
		//D�finition des deux cellules encadrant le flux en i+1/2 : la cellule c � gauche (en i) et la cellule cd � droite (en i+1) 
		Cellule1D cd = grille[i+1]; 
		Cellule1D c = grille[i]; 
 
			//Calcul des variables de Roe 
			double roe = sqrt(cd.rho/c.rho); 
			double rhor = roe*c.rho; 
			double ur = (roe*cd.u+c.u)/(1.+roe); 
			double Hr = (roe*(cd.rho*cd.u*cd.u/2.+cd.p*gam/(gam-1.))/cd.rho+(c.rho*c.u*c.u/2.+c.p*gam/(gam-1.))/c.rho)/(1.+roe); 
			double ur2 = ur*ur; 
			double cr2 = (gam-1.)*(Hr-ur2/2.); 
			
			//Test sur la vitesse du son 
			if(cr2<eps){ 
				cout << "i=" << i << " vitesse du son n�gative : c2=" << cr2 << endl; 
				cout << "t=" << t << endl; 
				cout << "c.p=" << c.p << endl; 
				cout << "c.rho=" << c.rho << endl; 
				cout << "c.u=" << c.u << endl; 
				cout << "cd.p=" << cd.p << endl; 
				cout << "cd.rho=" << cd.rho << endl; 
				cout << "cd.u=" << cd.u << endl; 
				cout << "ur=" << ur << endl; 
				cout << "ur2=" << ur2 << endl; 
				cout << "Hr=" << Hr << endl;
				getchar(); 
			} 
			
			//Vitesse du son 
			double cr = sqrt(cr2); 
			
			//Valeurs propres 
			c.lambda[0] = ur-cr; 
			c.lambda[1] = ur; 
			c.lambda[2] = ur+cr; 
			
			//Calcul des diff�rences entre Wd et Wg 
			double drho = cd.rho - c.rho; 
			double du = cd.u - c.u; 
			double dp = cd.p - c.p; 
			
			//Calcul des deltaV (diff�rences entre Wd et Wg dans la base des vecteurs propres du syst�me) 
			double ros2c = rhor/cr/2.; 
			c.delw[0] = dp/cr2/2. - ros2c*du; 
			c.delw[1] = drho - dp/cr2; 
			c.delw[2] = dp/cr2/2. + ros2c*du; 
			
			//Calcul de la correction compl�te dans la base des vecteurs propres 
			double xnu[3]; 
			for(int l=0;l<3;l++){ 
				xnu[l]  = sigma*abs(c.lambda[l]); 
				c.delwnu[l] = abs(c.lambda[l])*(1.-xnu[l])*c.delw[l]; 
				//Calcul des coefficients correctifs pour l'ordre sup�rieur 
				c.cf2[l]  = c.co[1]*abs(c.lambda[l])*(1.-xnu[l]); 
				c.cf3[l]  = c.co[2]*c.cf2[l]*(1.+xnu[l])/3.; 
				c.cf4[l]  = c.co[3]*c.cf3[l]*(xnu[l]-2.)/4.; 
				c.cf5[l]  = c.co[4]*c.cf4[l]*(xnu[l]+2.)/5.; 
				c.cf6[l]  = c.co[5]*c.cf5[l]*(xnu[l]-3.)/6.; 
				c.cf7[l]  = c.co[6]*c.cf6[l]*(xnu[l]+3.)/7.; 
				c.cf8[l]  = c.co[7]*c.cf7[l]*(xnu[l]-4.)/8.; 
				c.cf9[l]  = c.co[8]*c.cf8[l]*(xnu[l]+4.)/9.; 
				c.cf10[l] = c.co[9]*c.cf9[l]*(xnu[l]-5.)/10.; 
				c.cf11[l] = c.co[10]*c.cf10[l]*(xnu[l]+5.)/11.; 
			} 
			
			//Calcul des correctifs centr�s 
			for(int l=0;l<3;l++){ 
				c.psic0[l] = (c.cf2[l]-2.*c.cf4[l]+6.*c.cf6[l]-20.*c.cf8[l]+70.*c.cf10[l])*c.delw[l]; 
				c.psic1[l] = (c.cf4[l]-4.*c.cf6[l]+15.*c.cf8[l]-56.*c.cf10[l])*c.delw[l]; 
				c.psic2[l] = (c.cf6[l]-6.*c.cf8[l]+28.*c.cf10[l])*c.delw[l]; 
				c.psic3[l] = (c.cf8[l]-8.*c.cf10[l])*c.delw[l]; 
				c.psic4[l] = (c.cf10[l])*c.delw[l]; 
				//Calcul des correctifs d�centr�s 
				c.psid0[l] = (126.*c.cf11[l]-35.*c.cf9[l]+10.*c.cf7[l]-3.*c.cf5[l]+c.cf3[l])*c.delw[l]; 
				c.psid1[l] = (84.*c.cf11[l]-21.*c.cf9[l]+5.*c.cf7[l]-c.cf5[l])*c.delw[l]; 
				c.psid2[l] = (36.*c.cf11[l]-7.*c.cf9[l]+c.cf7[l])*c.delw[l]; 
				c.psid3[l] = (9.*c.cf11[l]-c.cf9[l])*c.delw[l]; 
				c.psid4[l] = (c.cf11[l])*c.delw[l]; 
			} 
			
			//Calcul des vecteurs propres �� gauche 
			c.vpr[0][0] = 1.; 
			c.vpr[1][0] = c.lambda[0]; 
			c.vpr[2][0] = Hr-ur*cr; 
			c.vpr[0][1] = 1.; 
			c.vpr[1][1] = c.lambda[1]; 
			c.vpr[2][1] = ur2/2.; 
			c.vpr[0][2] = 1.; 
			c.vpr[1][2] = c.lambda[2]; 
			c.vpr[2][2] = Hr+ur*cr; 
			
			//Calcul des corrections dans la base des vecteurs propres du syst�me 
			for(int l=0;l<3;l++){ 
				c.psic0r[l] = 0.; 
				c.psic1r[l] = 0.; 
				c.psic2r[l] = 0.; 
				c.psic3r[l] = 0.; 
				c.psic4r[l] = 0.; 
				c.psid0r[l] = 0.; 
				c.psid1r[l] = 0.; 
				c.psid2r[l] = 0.; 
				c.psid3r[l] = 0.; 
				c.psid4r[l] = 0.; 
			} 
			for(int m=0;m<3;m++){ 
				for(int l=0;l<3;l++){ 
					c.psic0r[m] += c.psic0[l]*c.vpr[m][l]; 
					c.psic1r[m] += c.psic1[l]*c.vpr[m][l]; 
					c.psic2r[m] += c.psic2[l]*c.vpr[m][l]; 
					c.psic3r[m] += c.psic3[l]*c.vpr[m][l]; 
					c.psic4r[m] += c.psic4[l]*c.vpr[m][l]; 
					c.psid0r[m] += c.psid0[l]*c.vpr[m][l]; 
					c.psid1r[m] += c.psid1[l]*c.vpr[m][l]; 
					c.psid2r[m] += c.psid2[l]*c.vpr[m][l]; 
					c.psid3r[m] += c.psid3[l]*c.vpr[m][l]; 
					c.psid4r[m] += c.psid4[l]*c.vpr[m][l]; 
				} 
			} 
			grille[i] = c; 

	} 
	//Fin de la boucle de calcul des valeurs initiales 
	
	//Boucle de calcul des indicateurs de monotonicit� 
	for(int l=0;l<3;l++){ 
		for(int i=1;i<N+2*marge-1;i++){ 
			Cellule1D c = grille[i]; 
			Cellule1D cg = grille[i-1]; 

				c.am[l] = c.lambda[l]*c.delw[l]-cg.lambda[l]*cg.delw[l]; 

		} 
		//Calcul de dj^m4 
		for(int i=1;i<N+2*marge-2;i++){ 
			Cellule1D c = grille[i]; 
			Cellule1D cd = grille[i+1]; 
			 
				double z1 = 4.*c.am[l]-cd.am[l]; 
				double z2 = 4.*cd.am[l]-c.am[l]; 
				double z3 = c.am[l]; 
				double z4 = cd.am[l]; 
				c.am1[l] = (sign1d(z1)+sign1d(z2))/2.*abs((sign1d(z1)+sign1d(z3))/2.)*(sign1d(z1)+sign1d(z4))/2.*min(abs(z1),min(abs(z2),min(abs(z3),abs(z4)))); 
				grille[i] = c; 
		}
	} 
	
	//Boucle de calcul de r+ et r- 
	for(int l=0;l<3;l++){ 
		for(int i=marge;i<N+2*marge-4;i++){ 
			Cellule1D c = grille[i]; 
			Cellule1D cd = grille[i+1]; 
			Cellule1D cg = grille[i-1]; 
 
				c.rp[l] = sign1d(c.delw[l])*sign1d(cg.delw[l])*(abs(cg.delw[l])+eps)/(abs(c.delw[l])+eps); 
				c.rm[l] = sign1d(c.delw[l])*sign1d(cd.delw[l])*(abs(cd.delw[l])+eps)/(abs(c.delw[l])+eps); 
				//Corrections d'ordre superieur 
				Cellule1D cg2 = grille[i-2]; 
				Cellule1D cg3 = grille[i-3]; 
				Cellule1D cg4 = grille[i-4]; 
				Cellule1D cg5 = grille[i-5]; 
				Cellule1D cd2 = grille[i+2]; 
				Cellule1D cd3 = grille[i+3]; 
				Cellule1D cd4 = grille[i+4]; 
				c.psid[l] = -c.psid0[l]+cg.psid0[l]+cd.psid1[l]-cg2.psid1[l]-cd2.psid2[l]+cg3.psid2[l]+cd3.psid3[l]-cg4.psid3[l]-cd4.psid4[l]+cg5.psid4[l]; 
				grille[i] = c; 

		} 
	} 
	
	//Boucle de calcul des flux 
	for(int i=marge;i<N+marge;i++){ 
		//Cellule de r�f�rence 
		Cellule1D c = grille[i]; 
		//Cellules voisines 
		Cellule1D cg = grille[i-1]; 
		Cellule1D cg2 = grille[i-2]; 
		Cellule1D cg3 = grille[i-3]; 
		Cellule1D cg4 = grille[i-4]; 
		Cellule1D cd = grille[i+1]; 
		Cellule1D cd2 = grille[i+2]; 
		Cellule1D cd3 = grille[i+3]; 
		Cellule1D cd4 = grille[i+4]; 

			//Flux TVD a ajouter au flux centr� d�j�� calcul� 
			double tvd[3]; 
			double psict[3]; 
			
			//Initialisation 
			for(int l=0;l<3;l++){ 
				tvd[l] = 0.; 
				//Partie centr�e du sch�ma 
				psict[l] = c.psic0r[l] + cg.psic1r[l] + cd.psic1r[l] + cg2.psic2r[l] + cd2.psic2r[l] + cg3.psic3r[l] + cd3.psic3r[l] + cg4.psic4r[l] + cd4.psic4r[l]; 
			} 
			
			//Limiteur 
			double psic; 
			for(int l=0;l<3;l++){ 
				psic = c.psic0[l] + cg.psic1[l] + cd.psic1[l] + cg2.psic2[l] + cd2.psic2[l] + cg3.psic3[l] + cd3.psic3[l] + cg4.psic4[l] + cd4.psic4[l]; 
				
				//Partie d�centr�e du sch�ma suivant le sign1de des valeurs propres 
				double r; 
				double xnum; 
				double xnume; 
				int is; 
				double psi; 
				if(c.lambda[l]>0.){ 
					r = c.rp[l]; 
					xnum = sigma*abs(cg.lambda[l]); 
					xnume = max(xnum,eps); 
					is = 1; 
					psi = psic+c.psid[l]; 
				} else { 
					r = c.rm[l]; 
					xnum = sigma*abs(cd.lambda[l]); 
					xnume = max(xnum,eps); 
					is = -1; 
					psi = psic-cd.psid[l]; 
				} 
				
				double xnu = sigma*abs(c.lambda[l]); 
				xnu = max(xnu,eps); 
				
				//Calcul du limiteur TVD psitvd 
				psi = (double) sign1d(c.delwnu[l])*psi/(abs(c.delwnu[l]+eps)); 
				double psimax1 = 2.*r*(1.-xnume)/(xnu*(1.-xnu)); 
				double psimax2 = 2./(1.-xnu); 
				double psitvd = max(0.,min(psi,min(psimax1,psimax2))); 
				
				//Crit�re de monotonicit� 
				if((c.delwnu[l] != 0.) && (abs(psi-psitvd)>eps)){ 
					double dfo = psi*c.delwnu[l]/2.; 
					double dabsf = psimax2*c.delwnu[l]/2.; 
					double dful = psimax1*c.delwnu[l]/2.; 
					double dfmd = dabsf/2.-c.am1[l]/2.; 
					Cellule1D camont = grille[i-is];   //Cellule en amont d�centr�e 
					double dflc = dful/2.+((1.-xnume)/xnu)*camont.am1[l]/2.; 
					double dfmin = max(min(0.,min(dabsf,dfmd)),min(0.,min(dful,dflc))); 
					double dfmax = min(max(0.,max(dabsf,dfmd)),max(0.,max(dful,dflc))); 
					if((dfmin-dfo)*(dfmax-dfo)>0.){ 
						psi = psitvd; 
					} 
				} 
				
				//A d�commenter si on veut utiliser uniquement le TVD et pas le MP 
				//psi = psitvd; 
				
				//A d�commenter si on veut utiliser le sch�ma sans TVD ni MP 
				//psi = 0.; 
				
				double ctvd = psi*c.delwnu[l]/2.-abs(c.lambda[l])*c.delw[l]/2.; 
				for(int m=0;m<3;m++){ 
					tvd[m] += ctvd*c.vpr[m][l]; 
				} 
			} 
			//Fin du calcul de la correction TVD 
			
			// Calcul final du flux �� droite de la cellule c 
			for(int l=0;l<3;l++){
				c.fluxd[l] += tvd[l]; 
			} 
			
			grille[i] = c; 

	} 
	//Fin de la boucle sur les cellules 
	
	corent1D(sigma);

}
void corent1D(double sigma){
	
	//Initialisation des variables
	//Cellule cel;
	for(int i=0;i<N+2*marge;i++){
		Cellule1D cel = grille[i];
		cel.S = log(cel.p) - gam*log(cel.rho);
		cel.ve[0] = (1.-gam)/cel.p*cel.rhoE-(cel.S-gam-1.);
		cel.ve[1] = (gam-1.)/cel.p*cel.imp;
		cel.ve[2] = (1.-gam)*cel.rho/cel.p;
		cel.fe = -cel.imp*cel.S;
		cel.Qcd[0] = cel.Qcd[1] = cel.Qcd[2] = 0.;
		grille[i] = cel;
	}
	//Calcul du correcteur d'entropie selon x
	
	double df0 = 0., df1=0., df2=0.;
	double F0 = 0., F1=0., F2=0.;
	
	//Cellule c, cd;
	for(int i=0;i<N+2*marge-1;i++){
		Cellule1D c = grille[i];
		Cellule1D cd = grille[i+1];
		double alpha = 0.;
		//Calcul de pe
		double pe = (cd.ve[0]-c.ve[0])*(cd.rho-c.rho);
		pe += (cd.ve[1]-c.ve[1])*(cd.imp-c.imp);
		pe += (cd.ve[2]-c.ve[2])*(cd.rhoE-c.rhoE);
		
		//Calcul des differences du flux dans les deux cellules adjacentes
		df0 = (cd.imp-c.imp);
		df1 = (cd.rho*cd.u*cd.u+cd.p)-(c.rho*c.u*c.u+c.p);
		df2 = (cd.rhoE*cd.u+cd.p*cd.u)-(c.rhoE*c.u+c.p*c.u);
		
		//Calcul du flux centre
		F0 = 1./2.*(cd.imp+c.imp);
		F1 = 1./2.*((cd.rho*cd.u*cd.u+cd.p)+(c.rho*c.u*c.u+c.p));
		F2 = 1./2.*((cd.rhoE*cd.u+cd.p*cd.u)+(c.rhoE*c.u+c.p*c.u));
		
		//Calcul de qef
		double qef = cd.fe - c.fe;
		qef -= 0.5*(cd.ve[0]+c.ve[0])*df0;
		qef -= 0.5*(cd.ve[1]+c.ve[1])*df1;
		qef -= 0.5*(cd.ve[2]+c.ve[2])*df2;
		
		//Calcul de q-q*
		double qmqet = qef;
		qmqet += (cd.ve[0]-c.ve[0])*(c.fluxd[0]-F0);
		qmqet += (cd.ve[1]-c.ve[1])*(c.fluxd[1]-F1);
		qmqet += (cd.ve[2]-c.ve[2])*(c.fluxd[2]-F2);
		qmqet *= -2.*sigma;
		
		//Calcul de alpha
		//if(qmqet<0. && pe>eps){
			if(pe>eps){
				alpha = 2.*max(qef,0.)/pe;
				//alpha = max(-qmqet,0.)/pe;
			}
			
			//Calcul du correcteur d'entropie a droite
			
			c.Qcd[0] = alpha*(cd.rho-c.rho);
			c.Qcd[1] = alpha*(cd.imp-c.imp);
			c.Qcd[2] = alpha*(cd.rhoE-c.rhoE);
			
			for(int l=0;l<3;l++){
				c.fluxd[l] -= c.Qcd[l]; //modification flux
			}
			grille[i] = c;
	}
	
	//Modification des flux
	
}
//Sous-programme d'impression des r�sultats 

void impression(double t, int n){ 
	const char* fichier;
	{
		std::ostringstream oss;
		oss << "resultats/xt" << n << ".dat";
		string s = oss.str();
		fichier = s.c_str();
	}
	std::ofstream xt(fichier,std::ios::out | std::ios::trunc);
	for(int i=marge; i<N+marge; i++){ 
		Cellule1D c = grille[i]; 
		xt << c.x << " " << c.rho << " "<<c.p<<" "<<c.u <<endl; 
	}
}	
void impression_1d(double t, int n){ 
	const char* fichier;
	{
		std::ostringstream oss;
		oss << "resultats/fluide1d_temp" << n << ".dat";
		string s = oss.str();
		fichier = s.c_str();
	}
	std::ofstream xt(fichier,std::ios::out | std::ios::trunc);
	for(int i=marge; i<N+marge; i++){ 
		Cellule1D c = grille[i]; 
		xt << c.x << " " << c.rho << " "<<c.p<<" "<<c.u <<endl; 
	}
}

}; 

			
			
