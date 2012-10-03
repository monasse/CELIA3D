#include <iostream> 
#include <stdio.h> 
#include <fstream> 
#include <math.h> 
#include "fluide.hpp"
#include "abs.hpp"


#ifndef FLUIDE_CPP
#define FLUIDE_CPP 

using namespace std;  // espace de nom standard


//Definition des methodes de la classe Cellule 

//Constructeur 

Cellule::Cellule()
{
    
    x = y = z = 1.;
    
    dx = dy = dz = 1.;
    
    rho = rho1 = 1.;
    
    u = v = w = 0.;
    
    p = p1 = 1.;
    
    pdtx = pdty = pdtz = 0.;
    
    impx = rho*u; impy = rho*v; impz = rho*w;
    
    rhoE=rho*u*u/2. + rho*v*v/2.  + rho*w*w/2. + p/(gam-1.);
    
    Mrho = Mimpx = Mimpy = Mimpz = MrhoE = 0.;
    
    rho0 = impx0 = impy0 = impz0 = rhoE0 = 0.;
    
    cells = alpha = alpha1 = 0.;
    
    kappai = kappaj = kappak = kappai1 = kappaj1 = kappak1 = 0.;
    
    proche = proche1 = 0;
    
    lambda[0] = lambda[1] = lambda[2] = lambda[3] = lambda[4] = 1.; 
    
    xi = yj = zk = 1.;
    
    phi_x = phi_y = phi_z =0.;
    
    S = log(p) - gam*log(rho);
    
    ve[0] = (1.-gam)/p*rhoE-(S-gam-1.);
    ve[1] = (gam-1.)/p*impx;
    ve[2] = (gam-1.)/p*impy;
    ve[3] = (gam-1.)/p*impz;
    ve[4] = (1-gam)*rho/p;
    
    fex = -impx*S; fey = -impy*S; fez = -impz*S;
    
    
    for(int l=0;l<5;l++){
        
        rp[l] = rm[l] = 1.; 
        
        fluxi[l] = fluxj[l] = fluxk[l] = flux_modif[l] = 0.; 
        
        Qci[l] = Qcj[l] = Qck[l] = 0.; 
        
        delw[l] = delwnu[l] = 1.; 
        
        cf2[l] = cf3[l] = cf4[l] = cf5[l] = cf6[l] = cf7[l] = cf8[l] = cf9[l] = cf10[l] = cf11[l] = 0.;
        
        
        psic0[l] = psic1[l] = psic2[l] = psic3[l] = psic4[l] = 0.; 
        
        psid0[l] = psid1[l] = psid2[l] = psid3[l] = psid4[l] = 0.; 
        
        psic0r[l] = psic1r[l] = psic2r[l] = psic3r[l] = psic4r[l] = 0.; 
        
        psid0r[l] = psid1r[l] = psid2r[l] = psid3r[l] = psid4r[l] = 0.; 
        
        psid[l] = 0.; 
        am[l] = am1[l] = 0.;
        
        for(int m=0;m<5;m++){ 
            vpr[l][m] = 0.; 
        }
    } 
    
    ordre = ordremax; 
    
    for(int l=0; l<ordre;l++){ 
        co[l] = 1.; 
    } 
    
}


Cellule::Cellule(double x0, double y0, double z0)
{
    
    x=x0; y=y0; z=z0;
    
    
    dx = dy = dz = 1.;
    
    rho = rho1 = 1.;
    
    u = v = w = 0.;
    
    p = p1 = 1.;
    
    pdtx = pdty = pdtz = 0.;
    
    impx = rho*u; impy = rho*v; impz = rho*w;
    
    rhoE=rho*u*u/2. + rho*v*v/2.  + rho*w*w/2. + p/(gam-1.);
    
    Mrho = Mimpx = Mimpy = Mimpz = MrhoE = 0.;
    
    rho0 = impx0 = impy0 = impz0 = rhoE0 = 0.;
    
    cells = alpha = alpha1 = 0.;
    
    kappai = kappaj = kappak = kappai1 = kappaj1 = kappak1 = 0.;
    
    proche = proche1 = 0;
    
    lambda[0] = lambda[1] = lambda[2] = lambda[3] = lambda[4] = 1.; 
    
    xi = yj = zk = 1.;
    
    phi_x = phi_y = phi_z =0.;
    
    S = log(p) - gam*log(rho);
    
    ve[0] = (1.-gam)/p*rhoE-(S-gam-1.);
    ve[1] = (gam-1.)/p*impx;
    ve[2] = (gam-1.)/p*impy;
    ve[3] = (gam-1.)/p*impz;
    ve[4] = (1-gam)*rho/p;
    
    fex = -impx*S; fey = -impy*S; fez = -impz*S;
    
    
    for(int l=0;l<5;l++){
        
        rp[l] = rm[l] = 1.; 
        
        fluxi[l] = fluxj[l] = fluxk[l] = flux_modif[l] = 0.; 
        
        Qci[l] = Qcj[l] = Qck[l] = 0.; 
        
        delw[l] = delwnu[l] = 1.; 
        
        cf2[l] = cf3[l] = cf4[l] = cf5[l] = cf6[l] = cf7[l] = cf8[l] = cf9[l] = cf10[l] = cf11[l] = 0.;
        
        
        psic0[l] = psic1[l] = psic2[l] = psic3[l] = psic4[l] = 0.; 
        
        psid0[l] = psid1[l] = psid2[l] = psid3[l] = psid4[l] = 0.; 
        
        psic0r[l] = psic1r[l] = psic2r[l] = psic3r[l] = psic4r[l] = 0.; 
        
        psid0r[l] = psid1r[l] = psid2r[l] = psid3r[l] = psid4r[l] = 0.; 
        
        psid[l] = 0.; 
        am[l] = am1[l] = 0.;
        
        for(int m=0;m<5;m++){ 
            vpr[l][m] = 0.; 
        }
    } 
    
    ordre = ordremax; 
    
    for(int l=0; l<ordre;l++){ 
        co[l] = 1.; 
    } 
    
}



Cellule::Cellule(double x0, double y0, double z0, double dx0, double dy0, double dz0)
{
    
    x = x0; y = y0; z = z0;
    
    dx = dx0 ; dy = dy0 ; dz = dz0 ;
    
    
    rho = rho1 = 1.;
    
    u = v = w = 0.;
    
    p = p1 = 1.;
    
    pdtx = pdty = pdtz = 0.;
    
    impx = rho*u; impy = rho*v; impz = rho*w;
    
    rhoE=rho*u*u/2. + rho*v*v/2.  + rho*w*w/2. + p/(gam-1.);
    
    Mrho = Mimpx = Mimpy = Mimpz = MrhoE = 0.;
    
    rho0 = impx0 = impy0 = impz0 = rhoE0 = 0.;
    
    cells = alpha = alpha1 = 0.;
    
    kappai = kappaj = kappak = kappai1 = kappaj1 = kappak1 = 0.;
    
    proche = proche1 = 0;
    
    lambda[0] = lambda[1] = lambda[2] = lambda[3] = lambda[4] = 1.; 
    
    xi = yj = zk = 1.;
    
    phi_x = phi_y = phi_z =0.;
    
    S = log(p) - gam*log(rho);
    
    ve[0] = (1.-gam)/p*rhoE-(S-gam-1.);
    ve[1] = (gam-1.)/p*impx;
    ve[2] = (gam-1.)/p*impy;
    ve[3] = (gam-1.)/p*impz;
    ve[4] = (1.-gam)*rho/p;
    
    fex = -impx*S; fey = -impy*S; fez = -impz*S;
    
    
    for(int l=0;l<5;l++){
        
        rp[l] = rm[l] = 1.; 
        
        fluxi[l] = fluxj[l] = fluxk[l] = 0.; 
        
        Qci[l] = Qcj[l] = Qck[l] = 0.; 
        
        delw[l] = delwnu[l] = 1.; 
        
        cf2[l] = cf3[l] = cf4[l] = cf5[l] = cf6[l] = cf7[l] = cf8[l] = cf9[l] = cf10[l] = cf11[l] = 0.;
        
        
        psic0[l] = psic1[l] = psic2[l] = psic3[l] = psic4[l] = 0.; 
        
        psid0[l] = psid1[l] = psid2[l] = psid3[l] = psid4[l] = 0.; 
        
        psic0r[l] = psic1r[l] = psic2r[l] = psic3r[l] = psic4r[l] = 0.; 
        
        psid0r[l] = psid1r[l] = psid2r[l] = psid3r[l] = psid4r[l] = 0.; 
        
        psid[l] = 0.; 
        am[l] = am1[l] = 0.;
        
        for(int m=0;m<5;m++){ 
            vpr[l][m] = 0.; 
        }
    } 
    
    ordre = ordremax; 
    
    for(int l=0; l<ordre;l++){ 
        co[l] = 1.; 
    } 
    
}


Cellule & Cellule:: operator=(const Cellule &c){
    
    assert(this != &c);		
    x = c.x; y = c.y; z = c.z;
    
    dx = c.dx ; dy = c.dy ; dz = c.dz ;
    
    rho = c.rho; rho1 = c.rho1;
    
    u=c.u; v=c.v; w=c.w;
    
    p=c.p; p1 = c.p1;
    
    pdtx = c.pdtx; pdty = c.pdty; pdtz = c.pdtz;
    
    impx = c.impx; impy = c.impy; impz = c.impz;  
    
    phi_x = c.phi_x; phi_y = c.phi_y; phi_z = c.phi_z;
    
    rhoE=c.rhoE;
    
    Mrho = c.Mrho; Mimpx = c.Mimpx; Mimpy = c.Mimpy; Mimpz = c.Mimpz; MrhoE = c.MrhoE;
    
    rho0 = c.rho0; impx0 = c.impx0; impy0 = c.impy0; impz0 = c.impz0; rhoE0 = c.rhoE0;
    
    cells = c.cells; alpha = c.alpha; alpha1 = c.alpha1;
    
    kappai = c.kappai; kappaj = c.kappaj; kappak = c.kappak; 
    kappai1 = c.kappai1; kappaj1 = c.kappaj1; kappak1 = c.kappak1;
    
    proche = c.proche; proche1 = c.proche1;
    
    xi = c.xi; yj = c.yj; zk = c.zk;
    
    S = c.S;
    
    
    fex = c.fex; fey = c.fey; fez = c.fez;
    
    
    for(int l=0;l<5;l++){ 
        
        lambda[l] = c. lambda[l]; 
        rp[l] = c.rp[l]; rm[l] = c.rm[l]; 
        
         fluxi[l] = c.fluxi[l]; fluxj[l] = c.fluxj[l]; fluxk[l] = c.fluxk[l];  
	       //flux_modif[l] = c.flux_modif[l];
				 dtfxi[l] = c.dtfxi[l]; dtfyj[l] = c.dtfyj[l]; dtfzk[l] = c.dtfzk[l];
				 
        Qci[l] = c.Qci[l]; Qcj[l] = c.Qcj[l]; Qck[l] = c.Qck[l];
        
        delw[l] = c.delw[l]; delwnu[l] = c.delwnu[l];
        
        ve[l] = c.ve[l];  
        
        cf2[l] = c.cf2[l]; cf3[l] = c.cf3[l]; cf4[l] = c.cf4[l]; cf5[l] = c.cf5[l]; cf6[l] = c.cf6[l]; 
        cf7[l] = c.cf7[l]; cf8[l] = c.cf8[l]; cf9[l] = c.cf9[l]; cf10[l] = c.cf10[l]; cf11[l] = c.cf11[l]; 
        
        psic0[l] = c.psic0[l]; psic1[l] = c.psic1[l]; psic2[l] = c.psic2[l]; 
        psic3[l] = c.psic3[l]; psic4[l] = c.psic4[l]; 
        
        psid0[l] = c.psid0[l]; psid1[l] = c.psid1[l]; psid2[l] = c.psid2[l]; 
        psid3[l] = c.psid3[l]; psid4[l] = c.psid4[l]; 
        
        
        
        
        psic0r[l] = c.psic0r[l]; psic1r[l] = c.psic1r[l]; psic2r[l] = c.psic2r[l]; 
        psic3r[l] = c.psic3r[l]; psic4r[l] = c.psic4r[l]; 
        
        psid0r[l] = c.psid0r[l]; psid1r[l] = c.psid1r[l]; psid2r[l] = c.psid2r[l]; 
        psid3r[l] = c.psid3r[l]; psid4r[l] = c.psid4r[l]; 
        
        psid[l] = c.psid[l]; 
        am[l] = c.am[l]; am1[l] = c.am1[l]; 
        
        for(int m=0;m<5;m++){ 
            vpr[l][m] = c.vpr[l][m]; 
        }
        
    } 
    
    ordre = c.ordre; 
    
    for(int l=0; l<ordre;l++){ 
        co[l] = c.co[l]; 
    } 
    
    return *this; //renvoye le pointeur de l'object
}


//Destructeur
Cellule::~Cellule(){
}

bool Cellule :: is_in_cell(double x0,double y0, double z0)
{
    bool k = false;
    
    if( (abs(x0-x) <= dx*(1./2.)) && (abs(y0-y) <= dy*(1./2.)) && (abs(z0-z) <= dz*(1./2.)) )
	{
        k= true;
	}
    return k;
}

void Cellule :: Affiche (){
    
	  
   cout<< " x = "<< x<< " y = "<<y<< " z = "<<z<<endl;

		//if(std::abs(alpha-1.)<eps){
    //cout<< " ki = "<< kappai<< " kj = "<< kappaj<< " kk = "<< kappak<<endl;
    cout<<"alpha ="<<alpha<<endl;
// 		cout<<"pression ="<<p<<endl;
// 	  cout<<"pression parois x= "<<phi_x<<endl;
// 		cout<<"pression parois y= "<<phi_y<<endl;
// 		cout<<"pression parois z= "<<phi_z<<endl;
// 		cout<<"pression en x= "<<pdtx<<endl;
// 		cout<<"pression en y= "<<pdty<<endl;
// 		cout<<"pression en z= "<<pdtz<<endl;
// 		cout<<"u= "<<u<<endl;
// 		cout<<"v= "<<v<<endl;
// 		cout<<"w= "<<w<<endl;
// 		cout<<"E= "<<rhoE<<endl;
// 		
// 		
// 		cout<<"flux en i "<<endl;
// 		for(int i=0; i<5; i++)
// 			cout<<fluxi[i]<<endl;
// 		cout<<"flux en  j "<<endl;
// 		for(int i=0; i<5; i++)
// 			cout<<fluxj[i]<<endl;
// 		cout<<"flux sigma k "<<endl;
// 		for(int i=0; i<5; i++)
// 			cout<<fluxk[i]<<endl;
// 		
// 	  cout<<"flux sigma en i "<<endl;
// 	 for(int i=0; i<5; i++)
// 		 cout<<dtfxi[i]<<endl;
// 	 cout<<"flux en sigma j "<<endl;
// 	 for(int i=0; i<5; i++)
// 		 cout<<dtfyj[i]<<endl;
// 	 cout<<"flux en sigma k "<<endl;
// 	 for(int i=0; i<5; i++)
// 		 cout<<dtfzk[i]<<endl;
		//}
    
}


//Definition des methodes de la classe Grille 

//Constructeur

Grille::Grille(): grille(Nx+2*marge, vector< vector<Cellule> >(Ny+2*marge, vector<Cellule>(Nz+2*marge)) ){
	//Grille::Grille(){ 
	
// 		grille.resize(Nx+2*marge);
// 		for (int i = 0; i < Nx+2*marge; i++) {
// 			grille[i].resize(Ny+2*marge);
// 			
// 			for (int j = 0; j <Ny+2*marge; j++)
// 				grille[i][j].resize(Nz+2*marge);
// 		}
		
    x = X0; y = Y0; z = Z0;
    
    dx = deltax; dy = deltay; dz = deltaz;
    
    for(int i=0;i<Nx+2*marge;i++){
        for(int j=0;j<Ny+2*marge;j++){ 
            for(int k=0;k<Nz+2*marge;k++){ 
//                 Cellule c(x+dx/2.+(i-marge)*dx,y+dy/2.+(j-marge)*dy, z+dz/2.+(k-marge)*dz,dx, dy, dz); 
//                grille[i][j][k] = c;
							grille[i][j][k] = Cellule(x+dx/2.+(i-marge)*dx,y+dy/2.+(j-marge)*dy, z+dz/2.+(k-marge)*dz,dx, dy, dz);
            }
        } 
    }
} 

// Grille::Grille(int Nx0, int Ny0, int Nz0,
//                double dx0, double x0, double dy0,
//                double y0, double dz0, double z0):grille(Nx0+2*marge, vector< vector<Cellule> > (Ny0+2*marge, vector<Cellule>(Nz0+2*marge))){ 
//     
//     
//     x = x0; y = y0; z = z0;
//     
//     dx = dx0; dy = dy0; dz = dz0;
//     
//     for(int i=0;i<Nx0+2*marge;i++){
//         for(int j=0;j<Ny0+2*marge;j++){ 
//             for(int k=0;k<Nz0+2*marge;k++){ 
//                 Cellule c(x+dx/2.+(i-marge)*dx,y+dy/2.+(j-marge)*dy, z+dz/2.+(k-marge)*dz, dx, dy, dz); 
//                 grille[i][j][k] = c; 
//             }
//         } 
//     }
// }

//Destructeur
Grille::~Grille(){
}


void Grille:: affiche()
{
    int s=0;
		double vol=0.;
    for(int i=marge;i<Nx+marge;i++){
        for(int j=marge;j<Ny+marge;j++){ 
            for(int k=marge;k<Nz+marge;k++){ 
                s++;
								vol +=(dx*dy*dz)*grille[i][j][k].alpha;
                //cout<<"cellule i="<<i-marge<<"j= "<<j-marge<<endl;
									grille[i][j][k].Affiche();
            }
        }
    }
   cout<<"volume solide := "<<vol<<endl;
}

void Grille:: affiche(string r)
{
    int s=0;
		double vol=0.;
		for(int i=marge;i<Nx+marge;i++){
		  for(int j=marge;j<Ny+marge;j++){
			Cellule cb = grille[i][j][marge];
            for(int k=marge;k<Nz+marge;k++){ 
			  s++;
			  Cellule c = grille[i][j][k];
			  if(abs(c.w)>eps){
				cout << r << " " << c.x << " " << c.y << " " << c.z << " w=" << c.w << endl;
				//getchar();
			  }
			}
		  }
		}
}

//Accss a une cellule i 
Cellule Grille::cellule(int i,int j, int k){ 
    return grille[i][j][k]; 
}

Cellule Grille::in_cell(Point_3 p){
  int i,j,k;
  i = (int) (floor(CGAL::to_double((p.operator[](0)-x)/dx))+marge);
  j = (int) (floor(CGAL::to_double((p.operator[](1)-y)/dy))+marge);
  k = (int) (floor(CGAL::to_double((p.operator[](2)-z)/dz))+marge);
  if(i<0){
	i = 0;
  }
  if(i>Nx+2*marge-1){
	i = Nx+2*marge-1;
  }
  if(j<0){
	j = 0;
  }
  if(j>Ny+2*marge-1){
	j = Ny+2*marge-1;
  }
  if(k<0){
	k = 0;
  }
  if(k>Nz+2*marge-1){
	k = Nz+2*marge-1;
  }
  return cellule(i,j,k);
}


//Sous-programme de definition des conditions initiales 
void Grille::init(){
    
    Cellule c; 
    for(int i=0;i<Nx+2*marge;i++){
        for(int j=0;j<Ny+2*marge;j++){
            for(int k=0;k<Nz+2*marge;k++){
                c = grille[i][j][k];
                c.dx = deltax; c.dy = deltay; c.dz = deltaz;
                c.x = x+c.dx/2.+(i-marge)*c.dx;
                c.y = y+c.dy/2.+(j-marge)*c.dy;
                c.z = z+c.dz/2.+(k-marge)*c.dz;
                c.rho = Rho(c.x, c.y, c.z); 
                c.u = U(c.x,c.y, c.z);
                c.v = V(c.x,c.y, c.z);
                c.w = W(c.x,c.y, c.z);
                c.p   = P(c.x,c.y, c.z);
                c.impx = c.rho*c.u; c.impy = c.rho*c.v; c.impz = c.rho*c.w;
                c.rhoE= c.rho*c.u*c.u/2. + c.rho*c.v*c.v/2. + c.rho*c.w*c.w/2. + c.p/(gam-1.); 
                c.kappai = c.kappaj = c.kappak = c.alpha = 0.;
                grille[i][j][k] = c;
            }
        }
    } 
        
}



//Sous-programme de calcul du pas de temps fluide
double Grille::pas_temps(double t, double T){ 
    
    double dt = 10000.;
    //Cellule c;
    //Restriction CFL sur le fluide
    for(int i=marge;i<Nx+marge;i++){
        for(int j=marge;j<Ny+marge;j++){
            for(int k=marge;k<Nz+marge;k++){
              Cellule c = grille[i][j][k]; 
                double c2 = gam*c.p/c.rho;
                double dt1 = cfl*min(c.dx/(sqrt(c2)+abs(c.u)),min(c.dy/(sqrt(c2)+abs(c.v)), c.dz/(sqrt(c2)+abs(c.w))));
                dt = min(dt,dt1); 
            }
        } 
    }
    dt = min(dt,T-t); 
    return dt; 
}

//Calcul des quantites conservatives totales 
double Grille::Masse(){ 
    double m = 0.; 
    //Cellule c; 
    for(int i=marge;i<Nx+marge;i++){
        for(int j=marge;j<Ny+marge;j++){
            for(int k=marge;k<Nz+marge;k++){
							Cellule c = grille[i][j][k]; 
                m += c.rho*c.dx*c.dy*c.dz*(1.-c.alpha);
            }
        }
    } 
    return m; 
} 

double Grille::Impulsionx(){ 
    double impx = 0.; 
    //Cellule c ;
    for(int i=marge;i<Nx+marge;i++){
        for(int j=marge;j<Ny+marge;j++){
            for(int k=marge;k<Nz+marge;k++){
							Cellule c = grille[i][j][k]; 
                impx += c.impx*c.dx*c.dy*c.dz*(1.-c.alpha);
            }
        }
    } 
    return impx; 
}

double Grille::Impulsiony(){ 
    double impy = 0.; 
    //Cellule c ;  
    for(int i=marge;i<Nx+marge;i++){
        for(int j=marge;j<Ny+marge;j++){
            for(int k=marge;k<Nz+marge;k++){
							Cellule c = grille[i][j][k];  
                impy += c.impy*c.dx*c.dy*c.dz*(1.-c.alpha);
            }
        }
    } 
    return impy; 
}

double Grille::Impulsionz(){ 
    double impz = 0.; 
    //Cellule c ;
    for(int i=marge;i<Nx+marge;i++){
        for(int j=marge;j<Ny+marge;j++){
            for(int k=marge;k<Nz+marge;k++){
							Cellule c = grille[i][j][k];  
                impz += c.impz*c.dx*c.dy*c.dz*(1.-c.alpha);
            }
        }
    } 
    return impz; 
}

double Grille::Energie(){ 
    double E = 0.; 
   // Cellule c ;
    for(int i=marge;i<Nx+marge;i++){ 
        for(int j=marge;j<Ny+marge;j++){ 
            for(int k=marge;k<Nz+marge;k++){
							Cellule  c = grille[i][j][k]; 
                E += c.rhoE*c.dx*c.dy*c.dz*(1.-c.alpha); 
            }
        } 
    }
    return E; 
} 

//Fonction qui melange les cellules a densite ou pression negative 
void Grille::melange(const double dt){ 
	for(int i=marge;i<Nx+marge;i++){
		for(int j=marge;j<Ny+marge;j++){
			for(int k=marge;k<Nz+marge;k++){
				Cellule c = grille[i][j][k];
				if(c.rho<0. && abs(c.alpha-1.)>eps){
					cout << "densite negative en : " << c.x << " " << c.y << " " << c.z <<" rho " << c.rho << endl;
					c.rho = c.rho0;
					c.u = c.impx0/c.rho;
					c.v = c.impy0/c.rho;
					c.w = c.impz0/c.rho;
					c.p = (gam-1)*(c.rhoE0-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2.-c.rho*c.w*c.w/2.);
					c.impx = c.impx0;
					c.impy = c.impy0;
					c.impz = c.impz0;
					c.rhoE = c.rhoE0;
				}
				if(c.p<0. && abs(c.alpha-1.)>eps){
					cout << "pression negative en : " << c.x << " " << c.y << " " << c.z << " p " << c.p << endl;
					c.rho = c.rho0;
					c.u = c.impx0/c.rho;
					c.v = c.impy0/c.rho;
					c.w = c.impz0/c.rho;
					c.p = (gam-1)*(c.rhoE0-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2.-c.rho*c.w*c.w/2.);
					c.impx = c.impx0;
					c.impy = c.impy0;
					c.impz = c.impz0;
					c.rhoE = c.rhoE0;
					
				}
				grille[i][j][k] = c;
				
			}
		}
	}
}


//Resolution des equations pour le fluide dans la direction x

void Grille::solve_fluidx(const double dt){ 
    
   const double sigma = dt/dx; 
    double dw1 =0., dw2=0., dw3=0., dw4=0., dw5 = 0.;
    
    //Calcul des variables au temps t+dt
    //Cellule c, ci; 
    for(int i=1;i<Nx+2*marge-1;i++){ 
        for(int j=1;j<Ny+2*marge-1;j++){ 
            for(int k=1;k<Nz+2*marge-1;k++){
                
							Cellule c = grille[i][j][k]; 
							Cellule ci = grille[i-1][j][k];    //Cellule  i-1
                
                //Stockage de la pression utilise pendant le pas de temps
                c.pdtx = dt*c.p;
                
                
                dw1 = -sigma*(c.fluxi[0]-ci.fluxi[0]); 
                dw2 = -sigma*(c.fluxi[1]-ci.fluxi[1]); 
                dw3 = -sigma*(c.fluxi[2]-ci.fluxi[2]); 
                dw4 = -sigma*(c.fluxi[3]-ci.fluxi[3]);
                dw5 = -sigma*(c.fluxi[4]-ci.fluxi[4]);
                
                for(int l=0; l<5; l++){
                    c.dtfxi[l] = sigma*c.fluxi[l];
                }
                c.rho  += dw1; c.impx += dw2; c.impy += dw3; c.impz += dw4; c.rhoE += dw5;
                
                c.u = c.impx/c.rho; c.v = c.impy/c.rho; c.w = c.impz/c.rho;
                c.p = (gam-1.)*(c.rhoE-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2.-c.rho*c.w*c.w/2.); 
               // c.pdtx= c.p;
                grille[i][j][k] = c; 
            }
        }
    }
    melange(dt);   
} 


//Resolution des equations pour le fluide dans la direction y

void Grille::solve_fluidy(const double dt){ 
    const double sigma = dt/dy; 
    double dw1 =0., dw2=0., dw3=0., dw4=0., dw5 = 0.;
    
    //Cellule c, cj;
    //Calcul des variables au temps t+dt 
    for(int i=1;i<Nx+2*marge-1;i++){ 
        for(int j=1;j<Ny+2*marge-1;j++){ 
            for(int k=1;k<Nz+2*marge-1;k++){ 
                
							Cellule c = grille[i][j][k]; 
							Cellule cj = grille[i][j-1][k];    //Cellule j-1
                
                
                //Stockage de la pression utilise pendant le pas de temps
                c.pdty = dt*c.p;
                
                dw1 = -sigma*(c.fluxj[0]-cj.fluxj[0]); 
                dw2 = -sigma*(c.fluxj[1]-cj.fluxj[1]); 
                dw3 = -sigma*(c.fluxj[2]-cj.fluxj[2]); 
                dw4 = -sigma*(c.fluxj[3]-cj.fluxj[3]);
                dw5 = -sigma*(c.fluxj[4]-cj.fluxj[4]);
                
                for(int l=0;l<5;l++){
                    c.dtfyj[l] = sigma*c.fluxj[l];
                }
                
                c.rho  += dw1; c.impx += dw2; c.impy += dw3; c.impz += dw4; c.rhoE += dw5;
                
                c.u = c.impx/c.rho; c.v = c.impy/c.rho; c.w = c.impz/c.rho;
                c.p = (gam-1.)*(c.rhoE-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2.-c.rho*c.w*c.w/2.);
                //c.pdty= c.p;
                grille[i][j][k] = c;
            }
        }
    }
    melange(dt);    
} 


void Grille::solve_fluidz(const double dt){
    
   const double sigma = dt/dz;
    double dw1 =0., dw2=0., dw3=0., dw4=0., dw5 = 0.;
    
    //Cellule c, ck;
    //Calcul des variables au temps t+dt 
    
    for(int i=1;i<Nx+2*marge-1;i++){
        for(int j=1;j<Ny+2*marge-1;j++){
            for(int k=1;k<Nz+2*marge-1;k++){
                
							Cellule c = grille[i][j][k]; 
							Cellule ck = grille[i][j][k-1];    //Cellule k-1
                
                
                //Stockage de la pression utilise pendant le pas de temps
                c.pdtz = dt*c.p;
                
                dw1 = -sigma*(c.fluxk[0]-ck.fluxk[0]); 
                dw2 = -sigma*(c.fluxk[1]-ck.fluxk[1]); 
                dw3 = -sigma*(c.fluxk[2]-ck.fluxk[2]); 
                dw4 = -sigma*(c.fluxk[3]-ck.fluxk[3]);
                dw5 = -sigma*(c.fluxk[4]-ck.fluxk[4]);
                
                for(int l=0;l<5;l++){
                    c.dtfzk[l] = sigma*c.fluxk[l];
                }
                
                c.rho  += dw1; c.impx += dw2; c.impy += dw3; c.impz += dw4; c.rhoE += dw5;
                c.u = c.impx/c.rho; c.v = c.impy/c.rho; c.w = c.impz/c.rho;
                c.p = (gam-1.)*(c.rhoE-c.rho*c.u*c.u/2.-c.rho*c.v*c.v/2.-c.rho*c.w*c.w/2.);
                //c.pdtz= c.p;
                grille[i][j][k] = c;
            }
        }
    }
    melange(dt);  
}


//Correction d'entropie selon x
void Grille::corentx(double sigma){
    
    //Initialisation des variables
    //Cellule cel;
	for(int i=0;i<Nx+2*marge;i++){
	    for(int j=0;j<Ny+2*marge;j++){
            for(int k=0;k<Nz+2*marge;k++){
							  Cellule cel = grille[i][j][k];
                cel.S = log(cel.p) - gam*log(cel.rho);
                cel.ve[0] = (1.-gam)/cel.p*cel.rhoE-(cel.S-gam-1.);
                cel.ve[1] = (gam-1.)/cel.p*cel.impx;
                cel.ve[2] = (gam-1.)/cel.p*cel.impy;
                cel.ve[3] = (gam-1.)/cel.p*cel.impz;
                cel.ve[4] = (1.-gam)*cel.rho/cel.p;
                cel.fex = -cel.impx*cel.S;
                cel.fey = -cel.impy*cel.S;
                cel.fez = -cel.impz*cel.S;
                cel.Qci[0] = cel.Qci[1] = cel.Qci[2] = cel.Qci[3] = cel.Qci[4] = 0.;
                grille[i][j][k] = cel;
            }
        }
    }
    //Calcul du correcteur d'entropie selon x
    
    double df0 = 0., df1=0., df2=0., df3=0., df4=0.;
    double F0 = 0., F1=0., F2=0., F3=0., F4=0.;
    
    //Cellule c, cd;
	for(int i=0;i<Nx+2*marge-1;i++){
	    for(int j=0;j<Ny+2*marge;j++){ 
            for(int k=0;k<Nz+2*marge;k++){
							Cellule c = grille[i][j][k];
							Cellule cd = grille[i+1][j][k];
                double alpha = 0.;
                //Calcul de pe
                double pe = (cd.ve[0]-c.ve[0])*(cd.rho-c.rho);
                pe += (cd.ve[1]-c.ve[1])*(cd.impx-c.impx);
                pe += (cd.ve[2]-c.ve[2])*(cd.impy-c.impy);
                pe += (cd.ve[3]-c.ve[3])*(cd.impz-c.impz);
                pe += (cd.ve[4]-c.ve[4])*(cd.rhoE-c.rhoE);
                
                //Calcul des differences du flux dans les deux cellules adjacentes
                df0 = (cd.impx-c.impx);
                df1 = (cd.rho*cd.u*cd.u+cd.p)-(c.rho*c.u*c.u+c.p);
                df2 = (cd.rho*cd.u*cd.v)-(c.rho*c.u*c.v);
                df3 = (cd.rho*cd.u*cd.w)-(c.rho*c.u*c.w);
                df4 = (cd.rhoE*cd.u+cd.p*cd.u)-(c.rhoE*c.u+c.p*c.u);
                
                //Calcul du flux centre
                F0 = 1./2.*(cd.impx+c.impx);
                F1 = 1./2.*((cd.rho*cd.u*cd.u+cd.p)+(c.rho*c.u*c.u+c.p));
                F2 = 1./2.*((cd.rho*cd.u*cd.v)+(c.rho*c.u*c.v));
                F3 = 1./2.*((cd.rho*cd.u*cd.w)+(c.rho*c.u*c.w));
                F4 = 1./2.*((cd.rhoE*cd.u+cd.p*cd.u)+(c.rhoE*c.u+c.p*c.u));
                
                //Calcul de qef
                double qef = cd.fex - c.fex;
                qef -= 0.5*(cd.ve[0]+c.ve[0])*df0;
                qef -= 0.5*(cd.ve[1]+c.ve[1])*df1;
                qef -= 0.5*(cd.ve[2]+c.ve[2])*df2;
                qef -= 0.5*(cd.ve[3]+c.ve[3])*df3;
                qef -= 0.5*(cd.ve[4]+c.ve[4])*df4;
                
                //Calcul de q-q*
                double qmqet = qef;
                qmqet += (cd.ve[0]-c.ve[0])*(c.fluxi[0]-F0);
                qmqet += (cd.ve[1]-c.ve[1])*(c.fluxi[1]-F1);
                qmqet += (cd.ve[2]-c.ve[2])*(c.fluxi[2]-F2);
                qmqet += (cd.ve[3]-c.ve[3])*(c.fluxi[3]-F3);
                qmqet += (cd.ve[4]-c.ve[4])*(c.fluxi[4]-F4);
                qmqet *= -2.*sigma;
                
                //Calcul de alpha
                //if(qmqet<0. && pe>eps){
                if(pe>eps){
                    alpha = 2.*max(qef,0.)/pe;
                    //alpha = max(-qmqet,0.)/pe;
                }
                
                //Calcul du correcteur d'entropie a droite
                
                c.Qci[0] = alpha*(cd.rho-c.rho);
                c.Qci[1] = alpha*(cd.impx-c.impx);
                c.Qci[2] = alpha*(cd.impy-c.impy);
                c.Qci[3] = alpha*(cd.impz-c.impz);
                c.Qci[4] = alpha*(cd.rhoE-c.rhoE);
                
                for(int l=0;l<5;l++){
                    c.fluxi[l] -= c.Qci[l]; //modification flux
                }
                grille[i][j][k] = c;
            }
        }
    }
    
    //Modification des flux
    
}//Fin de la correction d'entropie selon x

//Correction d'entropie selon y
void Grille::corenty(double sigma){
    
    //Initialisation des variables
   // Cellule cel;
	for(int i=0;i<Nx+2*marge;i++){
        for(int j=0;j<Ny+2*marge;j++){ 
            for(int k=0;k<Nz+2*marge;k++){
							Cellule cel = grille[i][j][k];
                cel.S = log(cel.p) - gam*log(cel.rho);
                cel.ve[0] = (1.-gam)/cel.p*cel.rhoE-(cel.S-gam-1.);
                cel.ve[1] = (gam-1.)/cel.p*cel.impx;
                cel.ve[2] = (gam-1.)/cel.p*cel.impy;
                cel.ve[3] = (gam-1.)/cel.p*cel.impz;
                cel.ve[4] = (1.-gam)*cel.rho/cel.p;
                cel.fex = -cel.impx*cel.S;
                cel.fey = -cel.impy*cel.S;
                cel.fez = -cel.impz*cel.S;
                cel.Qcj[0] = cel.Qcj[1] = cel.Qcj[2] = cel.Qcj[3] = cel.Qcj[4] = 0.;
                grille[i][j][k] = cel;
            }
        }
    }
    //Calcul du correcteur d'entropie selon y
    
   // Cellule c, ch;
    double df0 = 0., df1=0., df2=0., df3=0., df4=0.;
    double F0 = 0., F1=0., F2=0., F3=0., F4=0.;
    
	for(int i=0;i<Nx+2*marge;i++){
        for(int j=0;j<Ny+2*marge-1;j++){ 
            for(int k=0;k<Nz+2*marge;k++){
							Cellule c = grille[i][j][k];
							Cellule ch = grille[i][j+1][k];
                double alpha = 0.;
                //Calcul de pe
                double pe = (ch.ve[0]-c.ve[0])*(ch.rho-c.rho);
                pe += (ch.ve[1]-c.ve[1])*(ch.impx-c.impx);
                pe += (ch.ve[2]-c.ve[2])*(ch.impy-c.impy);
                pe += (ch.ve[3]-c.ve[3])*(ch.impz-c.impz);
                pe += (ch.ve[4]-c.ve[4])*(ch.rhoE-c.rhoE);
                
                //Calcul des differences du flux dans les deux cellules adjacentes
                df0 = (ch.impy-c.impy);
                df1 = (ch.rho*ch.u*ch.v)-(c.rho*c.u*c.v);
                df2 = (ch.rho*ch.v*ch.v+ch.p)-(c.rho*c.v*c.v+c.p);
                df3 = (ch.rho*ch.v*ch.w)-(c.rho*c.v*c.w);
                df4 = (ch.rhoE*ch.v+ch.p*ch.v)-(c.rhoE*c.v+c.p*c.v);
                
                //Calcul du flux centre
                F0 = 1./2.*(ch.impy+c.impy);
                F1 = 1./2.*((ch.rho*ch.u*ch.v)+(c.rho*c.u*c.v));
                F2 = 1./2.*((ch.rho*ch.v*ch.v+ch.p)+(c.rho*c.v*c.v+c.p));
                F3 = 1./2.*((ch.rho*ch.v*ch.w)+(c.rho*c.v*c.w));
                F4 = 1./2.*((ch.rhoE*ch.v+ch.p*ch.v)+(c.rhoE*c.v+c.p*c.v));
                
                //Calcul de qef
                double qef = ch.fey - c.fey;
                qef -= 0.5*(ch.ve[0]+c.ve[0])*df0;
                qef -= 0.5*(ch.ve[1]+c.ve[1])*df1;
                qef -= 0.5*(ch.ve[2]+c.ve[2])*df2;
                qef -= 0.5*(ch.ve[3]+c.ve[3])*df3;
                qef -= 0.5*(ch.ve[4]+c.ve[4])*df4;
                
                //Calcul de q-q*
                double qmqet = qef;
                qmqet += (ch.ve[0]-c.ve[0])*(c.fluxj[0]-F0);
                qmqet += (ch.ve[1]-c.ve[1])*(c.fluxj[1]-F1);
                qmqet += (ch.ve[2]-c.ve[2])*(c.fluxj[2]-F2);
                qmqet += (ch.ve[3]-c.ve[3])*(c.fluxj[3]-F3);
                qmqet += (ch.ve[3]-c.ve[3])*(c.fluxj[4]-F4);
                qmqet *= -2.*sigma;
                //Calcul de alpha
                //if(qmqet<0. && pe>eps){
                if(pe>eps){
                    alpha = 2.*max(qef,0.)/pe;
                    //alpha = max(-qmqet,0.)/pe;
                }
                //Calcul du correcteur d'entropie en haut
                c.Qcj[0] = alpha*(ch.rho-c.rho);
                c.Qcj[1] = alpha*(ch.impx-c.impx);
                c.Qcj[2] = alpha*(ch.impy-c.impy);
                c.Qcj[3] = alpha*(ch.impz-c.impz);
                c.Qcj[4] = alpha*(ch.rhoE-c.rhoE);
                
                for(int l=0;l<5;l++){
                    c.fluxj[l] -= c.Qcj[l]; //modification flux
                }
                grille[i][j][k] = c;
            }
        }
    }
    
}  //Fin de la correction d'entropie selon y


//Correction d'entropie selon z
void Grille::corentz(double sigma){
    //Initialisation des variables
    //Cellule cel;
	for(int i=0;i<Nx+2*marge;i++){
        for(int j=0;j<Ny+2*marge;j++){ 
            for(int k=0;k<Nz+2*marge;k++){
							  Cellule cel = grille[i][j][k];
                cel.S = log(cel.p) - gam*log(cel.rho);
                cel.ve[0] = (1.-gam)/cel.p*cel.rhoE-(cel.S-gam-1.);
                cel.ve[1] = (gam-1.)/cel.p*cel.impx;
                cel.ve[2] = (gam-1.)/cel.p*cel.impy;
                cel.ve[3] = (gam-1.)/cel.p*cel.impz;
                cel.ve[4] = (1.-gam)*cel.rho/cel.p;
                cel.fex = -cel.impx*cel.S;
                cel.fey = -cel.impy*cel.S;
                cel.fez = -cel.impz*cel.S;
                cel.Qck[0] = cel.Qck[1] = cel.Qck[2] = cel.Qck[3] = cel.Qck[4] = 0.;
                grille[i][j][k] = cel;
            }
        }
    }
    //Calcul du correcteur d'entropie selon z
    
   // Cellule c, ch;
    double df0 = 0., df1=0., df2=0., df3=0., df4=0.;
    double F0 = 0., F1=0., F2=0., F3=0., F4=0.;
    
	for(int i=0;i<Nx+2*marge;i++){
        for(int j=0;j<Ny+2*marge;j++){ 
            for(int k=0;k<Nz+2*marge-1;k++){
							  Cellule c = grille[i][j][k];
							  Cellule ch = grille[i][j][k+1];
                double alpha = 0.;
                //Calcul de pe
                double pe = (ch.ve[0]-c.ve[0])*(ch.rho-c.rho);
                pe += (ch.ve[1]-c.ve[1])*(ch.impx-c.impx);
                pe += (ch.ve[2]-c.ve[2])*(ch.impy-c.impy);
                pe += (ch.ve[3]-c.ve[3])*(ch.impz-c.impz);
                pe += (ch.ve[4]-c.ve[4])*(ch.rhoE-c.rhoE);
                
                //Calcul des differences du flux dans les deux cellules adjacentes
                df0 = (ch.impz-c.impz);
                df1 = (ch.rho*ch.u*ch.w)-(c.rho*c.u*c.w);
                df2 = (ch.rho*ch.v*ch.w)-(c.rho*c.v*c.w);
                df3 = (ch.rho*ch.w*ch.w+ch.p)-(c.rho*c.w*c.w+c.p);
                df4 = (ch.rhoE*ch.w+ch.p*ch.w)-(c.rhoE*c.w+c.p*c.w);
                
                //Calcul du flux centre
                F0 = 1./2.*(ch.impz+c.impz);
                F1 = 1./2.*((ch.rho*ch.u*ch.w)+(c.rho*c.u*c.w));
                F3 = 1./2.*((ch.rho*ch.w*ch.w+ch.p)+(c.rho*c.w*c.w+c.p));
                F2 = 1./2.*((ch.rho*ch.v*ch.w)+(c.rho*c.v*c.w));
                F4 = 1./2.*((ch.rhoE*ch.w+ch.p*ch.w)+(c.rhoE*c.w+c.p*c.w));
                
                //Calcul de qef
                double qef = ch.fez - c.fez;
                qef -= 0.5*(ch.ve[0]+c.ve[0])*df0;
                qef -= 0.5*(ch.ve[1]+c.ve[1])*df1;
                qef -= 0.5*(ch.ve[2]+c.ve[2])*df2;
                qef -= 0.5*(ch.ve[3]+c.ve[3])*df3;
                qef -= 0.5*(ch.ve[4]+c.ve[4])*df4;
                
                //Calcul de q-q*
                double qmqet = qef;
                qmqet += (ch.ve[0]-c.ve[0])*(c.fluxk[0]-F0);
                qmqet += (ch.ve[1]-c.ve[1])*(c.fluxk[1]-F1);
                qmqet += (ch.ve[2]-c.ve[2])*(c.fluxk[2]-F2);
                qmqet += (ch.ve[3]-c.ve[3])*(c.fluxk[3]-F3);
                qmqet += (ch.ve[3]-c.ve[3])*(c.fluxk[4]-F4);
                qmqet *= -2.*sigma;
                //Calcul de alpha
                //if(qmqet<0. && pe>eps){
                if(pe>eps){
                    alpha = 2.*max(qef,0.)/pe;
                    //alpha = max(-qmqet,0.)/pe;
                }
                //Calcul du correcteur d'entropie en haut
                c.Qck[0] = alpha*(ch.rho-c.rho);
                c.Qck[1] = alpha*(ch.impx-c.impx);
                c.Qck[2] = alpha*(ch.impy-c.impy);
                c.Qck[3] = alpha*(ch.impz-c.impz);
                c.Qck[4] = alpha*(ch.rhoE-c.rhoE);
                
                for(int l=0;l<5;l++){
                    c.fluxk[l] -= c.Qck[l]; //Modification des flux
                }	
                grille[i][j][k] = c;
            }
        }
    }
    
    
    
}  //Fin de la correction d'entropie selon z


//Calcul du flux en x entre les cellules 
void Grille::fnumx(const double sigma, double t){ 
   // Cellule c, ci, cg, cd, cg2, cg3, cg4, cg5, cd2, cd3, cd4;
    //Initialisation du flux au flux centre
    for(int i=0; i<Nx+2*marge-1; i++){ 
        for(int j=0; j<Ny+2*marge-1; j++){
            for(int k=0; k<Nz+2*marge-1; k++){
							Cellule c = grille[i][j][k];      //Cellule de reference pour le calcul du flux  
							Cellule ci = grille[i+1][j][k];   //Cellule en i
                
                //Calcul d'indicateurs de l'ordre 
                for(int l=0;l< c.ordre;l++){ 
                    c.co[l]=1.; 
                } 
                for(int l=c.ordre;l<ordremax;l++){ 
                    c.co[l]=0.;
                } 
                
                //partie flux centre
                c.fluxi[0] = (c.impx+ci.impx)/2.; 
                c.fluxi[1] = (c.rho*c.u*c.u+c.p+ci.rho*ci.u*ci.u+ci.p)/2.;
                c.fluxi[2] = (c.rho*c.u*c.v+ci.rho*ci.u*ci.v)/2.;
                c.fluxi[3] = (c.rho*c.u*c.w+ci.rho*ci.u*ci.w)/2.;
                c.fluxi[4] = ((c.rhoE+c.p)*c.u+(ci.rhoE+ci.p)*ci.u)/2.; 
                grille[i][j][k] = c;
            }
        }
    } 
    
    //Boucle sur les cellules : calcul des variables preliminaires au calcul des flux (variables de limitation) 
    for(int i=0; i<Nx+2*marge-1; i++){
        for(int j=0; j<Ny+2*marge-1; j++){
            for(int k=0; k<Nz+2*marge-1; k++){
                //Definition des deux cellules encadrant le flux en i+1/2 : 
                //la cellule c a gauche (en i) et la cellule cd a droite (en i+1) 
                Cellule ci = grille[i+1][j][k]; 
								Cellule c = grille[i][j][k]; 
                
                //Calcul des variables de Roe 
                double roe = sqrt(ci.rho/c.rho); 
                double rhor = roe*c.rho; 
                double ur = (roe*ci.u+c.u)/(1.+roe); 
                double vr = (roe*ci.v+c.v)/(1.+roe); 
                double wr = (roe*ci.w+c.w)/(1.+roe); 
                double Hr = (roe*(ci.rho*ci.u*ci.u/2.+ ci.rho*ci.v*ci.v/2.+ ci.rho*ci.w*ci.w/2. 
                            + ci.p*gam/(gam-1.))/ci.rho
                            + (c.rho*c.u*c.u/2.+c.rho*c.v*c.v/2.+ ci.rho*ci.w*ci.w/2.
                            + c.p*gam/(gam-1.))/c.rho)/(1.+roe);
                double ur2 = ur*ur;
                double vr2 = vr*vr;
                double wr2 = wr*wr;
                double cr2 = (gam-1.)*(Hr-ur2/2.-vr2/2.-wr2/2.); 
                
                //Test sur la vitesse du son 
                if(cr2<=0. && abs(c.alpha-1.)>eps){
				  cout << "calcul des flux selon x" << endl;
				  cout << "i=" << i << " j=" << j << " k=" << k<< " vitesse du son negative : c2=" << cr2 << endl;
				  cout << "x=" << c.x << " y=" << c.y << " z=" << c.z<< " alpha=" << c.alpha << endl;
				  cout << "t=" << t << endl; 
				  cout << "c.p=" << c.p << endl; 
				  cout << "c.rho=" << c.rho << endl; 
				  cout << "c.u=" << c.u << endl; 
				  cout << "c.v=" << c.v << endl; 
				  cout << "c.w=" << c.w << endl; 
				  cout << "ci.p=" << ci.p << endl; 
				  cout << "ci.rho=" << ci.rho << endl; 
                    cout << "ci.u=" << ci.u << endl;
                    cout << "ci.v=" << ci.v << endl;
                    cout << "ci.w=" << ci.w << endl;
                    cout << "ur=" << ur << endl; 
                    cout << "ur2=" << ur2 << endl; 
                    cout << "vr=" << vr << endl; 
                    cout << "vr2=" << vr2 << endl; 
                    cout << "wr2=" << wr2 << endl; 
                    cout << "Hr=" << Hr << endl;
                    getchar();
                } 
                
                //Vitesse du son 
                double cr = sqrt(cr2);
                
                //Valeurs propres 
                c.lambda[0] = ur-cr; 
                c.lambda[1] = ur;
                c.lambda[2] = ur;
                c.lambda[3] = ur;
                c.lambda[4] = ur+cr; 
                
                //Calcul des differences entre Wd et Wg 
                double drho = ci.rho - c.rho; 
                double du = ci.u - c.u;
                double dv = ci.v - c.v;
                double dw = ci.w - c.w;
                double dp = ci.p - c.p; 
                
                //Calcul des deltaV (differences entre Wd et Wg dans la base des vecteurs propres du systeme) 
                double ros2c = rhor/cr/2.; 
                c.delw[0] = dp/cr2/2. - ros2c*du; 
                c.delw[1] = drho - dp/cr2;
                c.delw[2] = 2.*ros2c*dv;
                c.delw[3] = 2.*ros2c*dw;   // a verifie!!!!!
                c.delw[4] = dp/cr2/2. + ros2c*du; 
                
                //Calcul de la correction complete dans la base des vecteurs propres 
                double xnu[5]; 
                for(int l=0;l<5;l++){ 
                    xnu[l]  = sigma*abs(c.lambda[l]); 
                    c.delwnu[l] = abs(c.lambda[l])*(1.-xnu[l])*c.delw[l]; 
                    //Calcul des coefficients correctifs pour l'ordre superieur 
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
                
                
                for(int l=0;l<5;l++){ 
                    //Calcul des correctifs centres 
                    c.psic0[l] = (c.cf2[l]-2.*c.cf4[l]+6.*c.cf6[l]-20.*c.cf8[l]+70.*c.cf10[l])*c.delw[l]; 
                    c.psic1[l] = (c.cf4[l]-4.*c.cf6[l]+15.*c.cf8[l]-56.*c.cf10[l])*c.delw[l]; 
                    c.psic2[l] = (c.cf6[l]-6.*c.cf8[l]+28.*c.cf10[l])*c.delw[l]; 
                    c.psic3[l] = (c.cf8[l]-8.*c.cf10[l])*c.delw[l]; 
                    c.psic4[l] = (c.cf10[l])*c.delw[l]; 
                    //Calcul des correctifs decentres 
                    c.psid0[l] = (126.*c.cf11[l]-35.*c.cf9[l]+10.*c.cf7[l]-3.*c.cf5[l]+c.cf3[l])*c.delw[l]; 
                    c.psid1[l] = (84.*c.cf11[l]-21.*c.cf9[l]+5.*c.cf7[l]-c.cf5[l])*c.delw[l]; 
                    c.psid2[l] = (36.*c.cf11[l]-7.*c.cf9[l]+c.cf7[l])*c.delw[l]; 
                    c.psid3[l] = (9.*c.cf11[l]-c.cf9[l])*c.delw[l]; 
                    c.psid4[l] = (c.cf11[l])*c.delw[l]; 
                } 
                
                //Calcul des vecteurs propres a gauche 
                c.vpr[0][0] = 1.; 
                c.vpr[1][0] = ur-cr; 
                c.vpr[2][0] = vr;
                c.vpr[3][0] = wr;
                c.vpr[4][0] = Hr-ur*cr; 
                
                c.vpr[0][1] = 1.; 
                c.vpr[1][1] = ur;
                c.vpr[2][1] = vr;
                c.vpr[3][1] = wr;
                c.vpr[4][1] = ur2/2. + vr2/2. + wr2/2.;
                
                c.vpr[0][2] = 0.; 
                c.vpr[1][2] = 0.;
                c.vpr[2][2] = cr;
                c.vpr[3][2] = 0.;
                c.vpr[4][2] = vr*cr;
                
                c.vpr[0][3] = 0.; 
                c.vpr[1][3] = 0.;
                c.vpr[2][3] = 0.;
                c.vpr[3][3] = cr;
                c.vpr[4][3] = wr*cr;
                
                c.vpr[0][4] = 1.; 
                c.vpr[1][4] = ur+cr;
                c.vpr[2][4] = vr;
                c.vpr[3][4] = wr;
                c.vpr[4][4] = Hr+ur*cr; 
                
                //Calcul des corrections dans la base des vecteurs propres du systeme 
                for(int l=0;l<5;l++){ 
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
                for(int m=0;m<5;m++){ 
                    for(int l=0;l<5;l++){ 
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
                grille[i][j][k] = c;
            }
        }
    }    //Fin de la boucle de calcul des valeurs initiales 
    
    
    //Boucle de calcul des indicateurs de monotonicite
    
    for(int l=0;l<5;l++){ 
        for(int i=1;i<Nx+2*marge-1;i++){ 
            for(int j=1;j<Ny+2*marge-1;j++){
                for(int k=1;k<Nz+2*marge-1;k++){
									Cellule c = grille[i][j][k]; 
									Cellule cg = grille[i-1][j][k]; 
                    c.am[l] = c.lambda[l]*c.delw[l]-cg.lambda[l]*cg.delw[l]; 
                    grille[i][j][k] = c; 
                }
            }
        } 
        //Calcul de dj^m4 
        for(int i=1;i<Nx+2*marge-2;i++){ 
            for(int j=0;j<Ny+2*marge;j++){
                for(int k=0;k<Nz+2*marge;k++){   
									Cellule c = grille[i][j][k]; 
									Cellule cd = grille[i+1][j][k]; 
                    double z1 = 4.*c.am[l]-cd.am[l]; 
                    double z2 = 4.*cd.am[l]-c.am[l]; 
                    double z3 = c.am[l]; 
                    double z4 = cd.am[l]; 
                    c.am1[l] = (sign(z1)+sign(z2))/2.*abs((sign(z1)+sign(z3))/2.)*(sign(z1)
                               + sign(z4))/2.*min(abs(z1),min(abs(z2),min(abs(z3),abs(z4))));
                    grille[i][j][k] = c; 
                }
            }
        } 
    } 
    
    //Boucle de calcul de r+ et r- 
    for(int l=0;l<5;l++){ 
        for(int i=marge;i<Nx+2*marge-4;i++){ 
            for(int j=marge;j<Ny+2*marge-4;j++){ 
                for(int k=marge;k<Nz+2*marge-4;k++){
									Cellule c = grille[i][j][k]; 
									Cellule cd = grille[i+1][j][k]; 
									Cellule cg = grille[i-1][j][k]; 
                    c.rp[l] = sign(c.delw[l])*sign(cg.delw[l])*(abs(cg.delw[l])+eps)/(abs(c.delw[l])+eps); 
                    c.rm[l] = sign(c.delw[l])*sign(cd.delw[l])*(abs(cd.delw[l])+eps)/(abs(c.delw[l])+eps); 
                    //Corrections d'ordre superieur 
                    Cellule  cg2 = grille[i-2][j][k]; 
										Cellule cg3 = grille[i-3][j][k]; 
										Cellule cg4 = grille[i-4][j][k]; 
										Cellule cg5 = grille[i-5][j][k]; 
										Cellule cd2 = grille[i+2][j][k]; 
										Cellule cd3 = grille[i+3][j][k]; 
										Cellule cd4 = grille[i+4][j][k]; 
                    c.psid[l] = -c.psid0[l]+cg.psid0[l]+cd.psid1[l]-cg2.psid1[l]-cd2.psid2[l]+cg3.psid2[l]
                                + cd3.psid3[l]-cg4.psid3[l]-cd4.psid4[l]+cg5.psid4[l]; 
                    grille[i][j][k] = c; 
                }
            }
        } 
    } 
    
    //Boucle de calcul des flux 
    for(int i=marge-1;i<Nx+marge;i++){ 
        for(int j=marge-1;j<Ny+marge;j++){ 
            for(int k=marge-1;k<Nz+marge;k++){ 
                //Cellule de reference 
                Cellule  c = grille[i][j][k]; 
                //Cellules voisines 
                Cellule cg = grille[i-1][j][k]; 
								Cellule cg2 = grille[i-2][j][k]; 
								Cellule cg3 = grille[i-3][j][k]; 
								Cellule cg4 = grille[i-4][j][k]; 
								Cellule cd = grille[i+1][j][k]; 
								Cellule cd2 = grille[i+2][j][k]; 
								Cellule cd3 = grille[i+3][j][k]; 
								Cellule cd4 = grille[i+4][j][k]; 
                
                //Flux TVD 
                double tvd[5]; 
                double psict[5]; 
                
                //Initialisation 
                for(int l=0; l<5; l++){ 
                    tvd[l] = 0.;
                    //Partie centre
                    psict[l] = c.psic0r[l] + cg.psic1r[l] + cd.psic1r[l] + cg2.psic2r[l] + cd2.psic2r[l] 
                               + cg3.psic3r[l] + cd3.psic3r[l] + cg4.psic4r[l] + cd4.psic4r[l]; 
                } 
                
                //Limiteur 
                double psic; 
                for(int l=0; l<5; l++){ 
                    psic = c.psic0[l] + cg.psic1[l] + cd.psic1[l] + cg2.psic2[l] + cd2.psic2[l] 
                           + cg3.psic3[l] + cd3.psic3[l] + cg4.psic4[l] + cd4.psic4[l]; 
                    //Partie descentre
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
                    psi = (double) sign(c.delwnu[l])*psi/(abs(c.delwnu[l]+eps)); 
                    double psimax1 = 2.*r*(1.-xnume)/(xnu*(1.-xnu)); 
                    double psimax2 = 2./(1.-xnu); 
                    double psitvd = max(0.,min(psi,min(psimax1,psimax2))); 
                    
                    //Critere de monotonicite
                    if((c.delwnu[l] != 0.) && (abs(psi-psitvd)>eps)){ 
                        double dfo = psi*c.delwnu[l]/2.; 
                        double dabsf = psimax2*c.delwnu[l]/2.; 
                        double dful = psimax1*c.delwnu[l]/2.; 
                        double dfmd = dabsf/2.-c.am1[l]/2.; 
                        Cellule camont = grille[i-is][j][k];   //Cellule en amont ddescentre
                        double dflc = dful/2.+((1.-xnume)/xnu)*camont.am1[l]/2.; 
                        double dfmin = max(min(0.,min(dabsf,dfmd)),min(0.,min(dful,dflc))); 
                        double dfmax = min(max(0.,max(dabsf,dfmd)),max(0.,max(dful,dflc))); 
                        if((dfmin-dfo)*(dfmax-dfo)>0.){ 
                            psi = psitvd; 
                        } 
                    } 
                    
                    //A decommenter si on veut utiliser uniquement le TVD et pas le MP 
                    //psi = psitvd; 
                    
                    //A decommenter si on veut utiliser le schema sans TVD ni MP 
                    //psi = 0.; 
                    
                    double ctvd = psi*c.delwnu[l]/2.-abs(c.lambda[l])*c.delw[l]/2.; 
                    for(int m=0;m<5;m++){ 
                        tvd[m] += ctvd*c.vpr[m][l];
                    } 
                } 
                //Fin du calcul de la correction TVD 
                
                // Calcul final du flux a droite de la cellule c 
                for(int l=0;l<5;l++){ 
                    c.fluxi[l] += tvd[l]; 
                }
                
                grille[i][j][k] = c; 
            }
        }
    } //Fin de la boucle sur les cellules
    
    
    // Correction d'entropie
    
    corentx(sigma);
    
    
    if(BC_x_in ==  1 || BC_x_out ==  1){
        for(int j=0;j<Ny+2*marge;j++){
            for(int k=0;k<Nz+2*marge;k++){
                
                if(BC_x_in ==  1){
                    Cellule c = grille[marge-1][j][k];
                    Cellule cp = grille[marge][j][k];
                    double p0 = (gam-1.)*(cp.rhoE0-1./2.*(cp.impx0*cp.impx0 + cp.impy0*cp.impy0 
                                + cp.impz0*cp.impz0)/cp.rho0);
                    c.fluxi[0] = 0.;
                    c.fluxi[1] = p0;
                    c.fluxi[2] = 0.;
                    c.fluxi[3] = 0.;
                    c.fluxi[4] = 0.;
                    
                    grille[marge-1][j][k] = c;
                }
                
                if(BC_x_out ==  1){
                    Cellule c2 = grille[Nx+marge-1][j][k];	  
                    double p02 = (gam-1.)*(c2.rhoE0-1./2.*(c2.impx0*c2.impx0 + c2.impy0*c2.impy0 
                                + c2.impz0*c2.impz0)/c2.rho0);;
                    c2.fluxi[0] = 0.;
                    c2.fluxi[1] = p02;
                    c2.fluxi[2] = 0.;
                    c2.fluxi[3] = 0.;
                    c2.fluxi[4] = 0.;
                    
                    grille[Nx+marge-1][j][k] = c2;
                }
                
            }
        }
    }
    
}


//Calcul du flux en y entre les cellules 
void Grille::fnumy(const double sigma, double t){ 
    
    //Cellule c, cj, camont, cg, cd, cg2, cg3, cg4, cg5, cd2, cd3, cd4;
    //Initialisation du flux au flux centre
    for(int i=0; i<Nx+2*marge-1; i++){
        for(int j=0; j<Ny+2*marge-1; j++){ 
            for(int k=0; k<Nz+2*marge-1; k++){
							Cellule c = grille[i][j][k];      //Cellule de reference pour le calcul du flux  
							Cellule cj = grille[i][j+1][k];   //Cellule en j
                //Calcul d'indicateurs de l'ordre 
                for(int l=0;l< c.ordre;l++){ 
                    c.co[l]=1.; 
                } 
                for(int l=c.ordre;l<ordremax;l++){ 
                    c.co[l]=0.;
                } 
                
                //partie flux centre  
                c.fluxj[0] = (c.impy+cj.impy)/2.; 
                c.fluxj[1] = (c.rho*c.u*c.v+cj.rho*cj.u*cj.v)/2.;
                c.fluxj[2] = (c.rho*c.v*c.v +c.p +cj.rho*cj.v*cj.v + cj.p)/2.;
                c.fluxj[3] = (c.rho*c.v*c.w+cj.rho*cj.v*cj.w)/2.;
                c.fluxj[4] = ((c.rhoE+c.p)*c.v+(cj.rhoE+cj.p)*cj.v)/2.; 
                grille[i][j][k] = c;
            }
        }
    } 
    
    //Boucle sur les cellules : calcul des variables preliminaires au calcul des flux (variables de limitation) 
    
    for(int i=0; i<Nx+2*marge-1; i++){
        for(int j=0; j<Ny+2*marge-1; j++){
            for(int k=0; k<Nz+2*marge-1; k++){
                //Definition des deux cellules encadrant le flux en i+1/2 : la cellule c a gauche (en i) et la cellule cd a droite (en i+1) 
                Cellule cj = grille[i][j+1][k]; 
								Cellule c = grille[i][j][k]; 
                
                //Calcul des variables de Roe 
                double roe = sqrt(cj.rho/c.rho); 
                double rhor = roe*c.rho; 
                double ur = (roe*cj.u+c.u)/(1.+roe); 
                double vr = (roe*cj.v+c.v)/(1.+roe); 
                double wr = (roe*cj.w+c.w)/(1.+roe); 
                double Hr = (roe*(cj.rho*cj.u*cj.u/2.+ cj.rho*cj.v*cj.v/2.+ cj.rho*cj.w*cj.w/2. 
                            + cj.p*gam/(gam-1.))/cj.rho
                            + (c.rho*c.u*c.u/2.+c.rho*c.v*c.v/2.+ cj.rho*cj.w*cj.w/2.
                            + c.p*gam/(gam-1.))/c.rho)/(1.+roe);
                double ur2 = ur*ur; 
                double vr2 = vr*vr;
                double wr2 = wr*wr;
                double cr2 = (gam-1.)*(Hr-ur2/2.-vr2/2.-wr2/2.); 
                
                //Test sur la vitesse du son 
                //if((cr2<=0.) && (c.alpha<1.) && (cd.alpha<1.)){
                if(cr2<=0. && abs(c.alpha-1.)>eps){
                    cout << "calcul des flux selon y" << endl;
                    cout << "i=" << i << " j=" << j <<" k = "<<k<< " vitesse du son negative : c2=" << cr2 << endl;
                    cout << "x=" << c.x << " y=" << c.y << " z=" << c.z << " alpha=" << c.alpha << endl;
                    cout << "t=" << t << endl; 
                    cout << "c.p=" << c.p << endl; 
                    cout << "c.rho=" << c.rho << endl; 
                    cout << "c.u=" << c.u << endl; 
                    cout << "c.v=" << c.v << endl; 
                    cout << "c.w=" << c.w << endl; 
                    cout << "cj.p=" << cj.p << endl; 
                    cout << "cj.rho=" << cj.rho << endl; 
                    cout << "cj.u=" << cj.u << endl;
                    cout << "cj.v=" << cj.v << endl;
                    cout << "cj.w=" << cj.w << endl;
                    cout << "ur=" << ur << endl; 
                    cout << "ur2=" << ur2 << endl; 
                    cout << "vr=" << vr << endl; 
                    cout << "vr2=" << vr2 << endl; 
                    cout << "wr2=" << wr2 << endl; 
                    cout << "Hr=" << Hr << endl;
                    getchar(); 
                } 
                
                //Vitesse du son 
                double cr = sqrt(cr2); 
                
                //Valeurs propres 
                c.lambda[0] = vr-cr; 
                c.lambda[1] = vr;
                c.lambda[2] = vr;
                c.lambda[3] = vr;
                c.lambda[4] = vr+cr; 
                
                //Calcul des differences entre Wd et Wg 
                double drho = cj.rho - c.rho; 
                double du = cj.u - c.u;
                double dv = cj.v - c.v;
                double dw = cj.w - c.w;
                double dp = cj.p - c.p; 
                
                //Calcul des deltaV (differences entre Wd et Wg dans la base des vecteurs propres du systeme) 
                double ros2c = rhor/cr/2.; 
                c.delw[0] = dp/cr2/2. - ros2c*dv; 
                c.delw[1] = drho - dp/cr2;
                c.delw[2] = 2.*ros2c*du;
                c.delw[3] = 2.*ros2c*dw;   // a verifie!!!!!
                c.delw[4] = dp/cr2/2. + ros2c*dv; 
                
                //Calcul de la correction complete dans la base des vecteurs propres 
                double xnu[5]; 
                for(int l=0;l<5;l++){ 
                    xnu[l]  = sigma*abs(c.lambda[l]); 
                    c.delwnu[l] = abs(c.lambda[l])*(1.-xnu[l])*c.delw[l]; 
                    //Calcul des coefficients correctifs pour l'ordre superieur 
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
                
                
                for(int l=0;l<5;l++){ 
                    //Calcul des correctifs centres 
                    c.psic0[l] = (c.cf2[l]-2.*c.cf4[l]+6.*c.cf6[l]-20.*c.cf8[l]+70.*c.cf10[l])*c.delw[l]; 
                    c.psic1[l] = (c.cf4[l]-4.*c.cf6[l]+15.*c.cf8[l]-56.*c.cf10[l])*c.delw[l]; 
                    c.psic2[l] = (c.cf6[l]-6.*c.cf8[l]+28.*c.cf10[l])*c.delw[l]; 
                    c.psic3[l] = (c.cf8[l]-8.*c.cf10[l])*c.delw[l]; 
                    c.psic4[l] = (c.cf10[l])*c.delw[l]; 
                    //Calcul des correctifs decentres 
                    c.psid0[l] = (126.*c.cf11[l]-35.*c.cf9[l]+10.*c.cf7[l]-3.*c.cf5[l]+c.cf3[l])*c.delw[l]; 
                    c.psid1[l] = (84.*c.cf11[l]-21.*c.cf9[l]+5.*c.cf7[l]-c.cf5[l])*c.delw[l]; 
                    c.psid2[l] = (36.*c.cf11[l]-7.*c.cf9[l]+c.cf7[l])*c.delw[l]; 
                    c.psid3[l] = (9.*c.cf11[l]-c.cf9[l])*c.delw[l]; 
                    c.psid4[l] = (c.cf11[l])*c.delw[l]; 
                } 
                
                //Calcul des vecteurs propres a gauche 
                c.vpr[0][0] = 1.; 
                c.vpr[1][0] = ur; 
                c.vpr[2][0] = vr-cr;
                c.vpr[3][0] = wr;
                c.vpr[4][0] = Hr-vr*cr; 
                
                c.vpr[0][1] = 1.; 
                c.vpr[1][1] = ur;
                c.vpr[2][1] = vr;
                c.vpr[3][1] = wr;
                c.vpr[4][1] = ur2/2. + vr2/2. + wr2/2.;
                
                c.vpr[0][2] = 0.; 
                c.vpr[1][2] = 0.;
                c.vpr[2][2] = cr;
                c.vpr[3][2] = 0.;
                c.vpr[4][2] = ur*cr;
                
                c.vpr[0][3] = 0.; 
                c.vpr[1][3] = 0.;
                c.vpr[2][3] = 0.;
                c.vpr[3][3] = cr;
                c.vpr[4][3] = wr*cr;
                
                c.vpr[0][4] = 1.; 
                c.vpr[1][4] = ur;
                c.vpr[2][4] = vr+cr;
                c.vpr[3][4] = wr;
                c.vpr[4][4] = Hr+vr*cr; 
                
                //Calcul des corrections dans la base des vecteurs propres du systeme 
                for(int l=0;l<5;l++){ 
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
                for(int m=0;m<5;m++){ 
                    for(int l=0;l<5;l++){ 
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
                grille[i][j][k] = c; 
            }
        }
    }    //Fin de la boucle de calcul des valeurs initiales 
    
    
    //Boucle de calcul des indicateurs de monotonicite
    for(int l=0;l<5;l++){ 
        for(int i=1;i<Nx+2*marge-1;i++){
            for(int j=1;j<Ny+2*marge-1;j++){ 
                for(int k=1;k<Nz+2*marge-1;k++){
									Cellule c = grille[i][j][k]; 
									Cellule cg = grille[i][j-1][k]; 
                    c.am[l] = c.lambda[l]*c.delw[l]-cg.lambda[l]*cg.delw[l]; 
                    grille[i][j][k] = c; 
                }
            }
        } 
        //Calcul de dj^m4 
        for(int i=0;i<Nx+2*marge;i++){
            for(int j=1;j<Ny+2*marge-2;j++){ 
                for(int k=0;k<Nz+2*marge;k++){   
									Cellule c = grille[i][j][k]; 
									Cellule  cd = grille[i][j+1][k]; 
                    double z1 = 4.*c.am[l]-cd.am[l]; 
                    double z2 = 4.*cd.am[l]-c.am[l]; 
                    double z3 = c.am[l]; 
                    double z4 = cd.am[l]; 
                    c.am1[l] = (sign(z1)+sign(z2))/2.*abs((sign(z1)+sign(z3))/2.)*(sign(z1)
                               + sign(z4))/2.*min(abs(z1),min(abs(z2),min(abs(z3),abs(z4))));
                    grille[i][j][k] = c; 
                }
            }
        } 
    } 
    //Boucle de calcul de r+ et r- 
    for(int l=0;l<5;l++){ 
        for(int i=marge;i<Nx+2*marge-4;i++){ 
            for(int j=marge;j<Ny+2*marge-4;j++){ 
                for(int k=marge;k<Nz+2*marge-4;k++){
									Cellule c = grille[i][j][k]; 
									Cellule cd = grille[i][j+1][k]; 
									Cellule cg = grille[i][j-1][k]; 
                    c.rp[l] = sign(c.delw[l])*sign(cg.delw[l])*(abs(cg.delw[l])+eps)/(abs(c.delw[l])+eps); 
                    c.rm[l] = sign(c.delw[l])*sign(cd.delw[l])*(abs(cd.delw[l])+eps)/(abs(c.delw[l])+eps); 
                    //Corrections d'ordre superieur 
                    Cellule  cg2 = grille[i][j-2][k]; 
										Cellule cg3 = grille[i][j-3][k]; 
										Cellule cg4 = grille[i][j-4][k]; 
										Cellule cg5 = grille[i][j-5][k]; 
										Cellule cd2 = grille[i][j+2][k]; 
										Cellule cd3 = grille[i][j+3][k]; 
										Cellule cd4 = grille[i][j+4][k]; 
                    c.psid[l] = - c.psid0[l]+cg.psid0[l]+cd.psid1[l]-cg2.psid1[l]-cd2.psid2[l]
                    + cg3.psid2[l]+cd3.psid3[l]-cg4.psid3[l]-cd4.psid4[l]+cg5.psid4[l];
                    grille[i][j][k] = c;
                }
            }
        } 
    } 
		
    //Boucle de calcul des flux 
    for(int i=marge-1;i<Nx+marge;i++){ 
        for(int j=marge-1;j<Ny+marge;j++){ 
            for(int k=marge-1;k<Nz+marge;k++){ 
                //Cellule de reference 
                Cellule   c = grille[i][j][k]; 
                //Cellules voisines 
                Cellule  cg = grille[i][j-1][k]; 
								Cellule cg2 = grille[i][j-2][k]; 
								Cellule cg3 = grille[i][j-3][k]; 
								Cellule cg4 = grille[i][j-4][k]; 
								Cellule cd = grille[i][j+1][k]; 
								Cellule cd2 = grille[i][j+2][k]; 
								Cellule cd3 = grille[i][j+3][k]; 
								Cellule cd4 = grille[i][j+4][k]; 
                
                //Flux TVD 
                double tvd[5]; 
                double psict[5]; 
                
                //Initialisation 
                for(int l=0;l<5;l++){ 
                    tvd[l] = 0.; 
                    //Partie centre
                    psict[l] = c.psic0r[l] + cg.psic1r[l] + cd.psic1r[l] + cg2.psic2r[l] + cd2.psic2r[l]
                               + cg3.psic3r[l] + cd3.psic3r[l] + cg4.psic4r[l] + cd4.psic4r[l]; 
                } 
                
                //Limiteur 
                double psic; 
                for(int l=0;l<5;l++){ 
                    psic = c.psic0[l] + cg.psic1[l] + cd.psic1[l] + cg2.psic2[l] + cd2.psic2[l] 
                           + cg3.psic3[l] + cd3.psic3[l] + cg4.psic4[l] + cd4.psic4[l]; 
                    
                    //Partie descentre
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
                    psi = (double) sign(c.delwnu[l])*psi/(abs(c.delwnu[l]+eps)); 
                    double psimax1 = 2.*r*(1.-xnume)/(xnu*(1.-xnu)); 
                    double psimax2 = 2./(1.-xnu); 
                    double psitvd = max(0.,min(psi,min(psimax1,psimax2))); 
                    
                    //Critere de monotonicite
                    if((c.delwnu[l] != 0.) && (abs(psi-psitvd)>eps)){ 
                        double dfo = psi*c.delwnu[l]/2.; 
                        double dabsf = psimax2*c.delwnu[l]/2.; 
                        double dful = psimax1*c.delwnu[l]/2.; 
                        double dfmd = dabsf/2.-c.am1[l]/2.; 
                        Cellule camont = grille[i-is][j][k];   //Cellule en amont ddescentre
                        double dflc = dful/2.+((1.-xnume)/xnu)*camont.am1[l]/2.; 
                        double dfmin = max(min(0.,min(dabsf,dfmd)),min(0.,min(dful,dflc))); 
                        double dfmax = min(max(0.,max(dabsf,dfmd)),max(0.,max(dful,dflc))); 
                        if((dfmin-dfo)*(dfmax-dfo)>0.){ 
                            psi = psitvd; 
                        } 
                    } 
                    
                    //A decommenter si on veut utiliser uniquement le TVD et pas le MP 
                    //psi = psitvd; 
                    
                    //A decommenter si on veut utiliser le schema sans TVD ni MP 
                    //psi = 0.; 
                    
                    double ctvd = psi*c.delwnu[l]/2.-abs(c.lambda[l])*c.delw[l]/2.; 
                    for(int m=0;m<5;m++){ 
                        tvd[m] += ctvd*c.vpr[m][l];
                    } 
                } 
                //Fin du calcul de la correction TVD 
                
                // Calcul final du flux a droite de la cellule c 
                for(int l=0;l<5;l++){ 
                    c.fluxj[l] += tvd[l]; 
                }
                
                grille[i][j][k] = c; 
            }
        }
    } //Fin de la boucle sur les cellules

    // Correction d'entropie
    corenty(sigma);
    
    
    if(BC_y_in ==  1 || BC_y_out ==  1){
        for(int i=0;i<Nx+2*marge;i++){
            for(int k=0;k<Nz+2*marge;k++){
                
                if(BC_y_in ==  1){
                    Cellule c = grille[i][marge-1][k];
                    Cellule cp = grille[i][marge][k];
                    double p0 = (gam-1.)*(cp.rhoE0 - 1./2.*(cp.impx0*cp.impx0 + cp.impy0*cp.impy0 
                                + cp.impz0*cp.impz0)/cp.rho0);
                    // double p = c.p;
                    c.fluxj[0] = 0.;
                    c.fluxj[1] = 0.;
                    c.fluxj[2] = p0;
                    c.fluxj[3] = 0.;
                    c.fluxj[4] = 0.;
                    
                    grille[i][marge-1][k] = c;
                }
                
                if(BC_y_out ==  1){
                    Cellule c2 = grille[i][Ny+marge-1][k];
                    //Cellule cp = grille[Nx+marge-1][j];
                    double p02 = (gam-1.)*(c2.rhoE0-1./2.*(c2.impx0*c2.impx0 + c2.impy0*c2.impy0 
                                 + c2.impz0*c2.impz0)/c2.rho0);
                    //double p2 = c2.p;
                    c2.fluxj[0] = 0.;
                    c2.fluxj[1] = 0.;
                    c2.fluxj[2] = p02;
                    c2.fluxj[3] = 0.;
                    c2.fluxj[4] = 0.;
                    
                    grille[i][Ny+marge-1][k] = c2;
                }
                
            }
        }
    }
    
}



//Calcul du flux en z entre les cellules 
void Grille::fnumz(const double sigma, double t){ 
    //Cellule c, ck, cd, cg, cg2, cg3, cg4, cg5, cd2, cd3, cd4, camont;
    //Initialisation du flux au flux centre
    for(int i=0; i<Nx+2*marge-1; i++){
        for(int j=0; j<Ny+2*marge-1; j++){
            for(int k=0; k<Nz+2*marge-1; k++){ 
							Cellule  c = grille[i][j][k];      //Cellule de reference pour le calcul du flux  
							Cellule  ck = grille[i][j][k+1];   //Cellule en k
                //Calcul d'indicateurs de l'ordre 
                for(int l=0;l< c.ordre;l++){ 
                    c.co[l]=1.; 
                } 
                for(int l=c.ordre;l<ordremax;l++){ 
                    c.co[l]=0.;
                }  
                
                //partie flux centre  
                c.fluxk[0] = (c.impz+ck.impz)/2.; 
                c.fluxk[1] = (c.rho*c.u*c.w+ck.rho*ck.u*ck.w)/2.;
                c.fluxk[2] = (c.rho*c.v*c.w+ck.rho*ck.v*ck.w)/2.;
                c.fluxk[3] = (c.rho*c.w*c.w +c.p +ck.rho*ck.w*ck.w + ck.p)/2.;
                c.fluxk[4] = ((c.rhoE+c.p)*c.w+(ck.rhoE+ck.p)*ck.w)/2.; 
                grille[i][j][k] = c;
            }
        }
    } 
    
    //Boucle sur les cellules : calcul des variables preliminaires au calcul des flux (variables de limitation) 
    
    for(int i=0; i<Nx+2*marge-1; i++){
        for(int j=0; j<Ny+2*marge-1; j++){
            for(int k=0; k<Nz+2*marge-1; k++){
                //Definition des deux cellules encadrant le flux en i+1/2 : la cellule c a gauche (en i) et la cellule cd a droite (en i+1) 
                Cellule ck = grille[i][j][k+1]; 
								Cellule c = grille[i][j][k]; 
                
                //Calcul des variables de Roe 
                double roe = sqrt(ck.rho/c.rho); 
                double rhor = roe*c.rho; 
                double ur = (roe*ck.u+c.u)/(1.+roe); 
                double vr = (roe*ck.v+c.v)/(1.+roe); 
                double wr = (roe*ck.w+c.w)/(1.+roe); 
                double Hr = (roe*(ck.rho*ck.u*ck.u/2.+ ck.rho*ck.v*ck.v/2.+ ck.rho*ck.w*ck.w/2. 
                            + ck.p*gam/(gam-1.))/ck.rho+(c.rho*c.u*c.u/2.+c.rho*c.v*c.v/2.+ ck.rho*ck.w*ck.w/2.
                            + c.p*gam/(gam-1.))/c.rho)/(1.+roe);
                double ur2 = ur*ur; 
                double vr2 = vr*vr;
                double wr2 = wr*wr;
                double cr2 = (gam-1.)*(Hr-ur2/2.-vr2/2.-wr2/2.); 
                
                //Test sur la vitesse du son 
                //if((cr2<=0.) && (c.alpha<1.) && (cd.alpha<1.)){
                if(cr2<=0. && abs(c.alpha-1.)>eps){
                    cout << "calcul des flux selon y" << endl;
                    cout << "i=" << i << " j=" << j << " vitesse du son negative : c2=" << cr2 << endl;
                    cout << "x=" << c.x << " y=" << c.y << " z=" << c.z << " alpha=" << c.alpha << endl;
                    cout << "t=" << t << endl; 
                    cout << "c.p=" << c.p << endl; 
                    cout << "c.rho=" << c.rho << endl; 
                    cout << "c.u=" << c.u << endl; 
                    cout << "c.v=" << c.v << endl; 
                    cout << "c.w=" << c.w << endl; 
                    cout << "ck.p=" << ck.p << endl; 
                    cout << "ck.rho=" << ck.rho << endl; 
                    cout << "ck.u=" << ck.u << endl;
                    cout << "ck.v=" << ck.v << endl;
                    cout << "ck.w=" << ck.w << endl;
                    cout << "ur=" << ur << endl; 
                    cout << "ur2=" << ur2 << endl; 
                    cout << "vr=" << vr << endl; 
                    cout << "vr2=" << vr2 << endl; 
                    cout << "wr2=" << wr2 << endl; 
                    cout << "Hr=" << Hr << endl;
                    getchar(); 
                } 
                
                //Vitesse du son 
                double cr = sqrt(cr2); 
                
                //Valeurs propres 
                c.lambda[0] = wr-cr; 
                c.lambda[1] = wr;
                c.lambda[2] = wr;
                c.lambda[3] = wr;
                c.lambda[4] = wr+cr; 
                
                //Calcul des differences entre Wd et Wg 
                double drho = ck.rho - c.rho; 
                double du = ck.u - c.u;
                double dv = ck.v - c.v;
                double dw = ck.w - c.w;
                double dp = ck.p - c.p; 
                
                //Calcul des deltaV (differences entre Wd et Wg dans la base des vecteurs propres du systeme) 
                double ros2c = rhor/cr/2.; 
                c.delw[0] = dp/cr2/2. - ros2c*dw; 
                c.delw[1] = drho - dp/cr2;
                c.delw[2] = 2.*ros2c*du;
                c.delw[3] = 2.*ros2c*dv;   // a verifie!!!!!
                c.delw[4] = dp/cr2/2. + ros2c*dw; 
                
                //Calcul de la correction complete dans la base des vecteurs propres 
                double xnu[5]; 
                for(int l=0;l<5;l++){ 
                    xnu[l]  = sigma*abs(c.lambda[l]); 
                    c.delwnu[l] = abs(c.lambda[l])*(1.-xnu[l])*c.delw[l]; 
                    //Calcul des coefficients correctifs pour l'ordre superieur 
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
                
                
                for(int l=0;l<5;l++){ 
                    //Calcul des correctifs centres 
                    c.psic0[l] = (c.cf2[l]-2.*c.cf4[l]+6.*c.cf6[l]-20.*c.cf8[l]+70.*c.cf10[l])*c.delw[l]; 
                    c.psic1[l] = (c.cf4[l]-4.*c.cf6[l]+15.*c.cf8[l]-56.*c.cf10[l])*c.delw[l]; 
                    c.psic2[l] = (c.cf6[l]-6.*c.cf8[l]+28.*c.cf10[l])*c.delw[l]; 
                    c.psic3[l] = (c.cf8[l]-8.*c.cf10[l])*c.delw[l]; 
                    c.psic4[l] = (c.cf10[l])*c.delw[l]; 
                    //Calcul des correctifs decentres 
                    c.psid0[l] = (126.*c.cf11[l]-35.*c.cf9[l]+10.*c.cf7[l]-3.*c.cf5[l]+c.cf3[l])*c.delw[l]; 
                    c.psid1[l] = (84.*c.cf11[l]-21.*c.cf9[l]+5.*c.cf7[l]-c.cf5[l])*c.delw[l]; 
                    c.psid2[l] = (36.*c.cf11[l]-7.*c.cf9[l]+c.cf7[l])*c.delw[l]; 
                    c.psid3[l] = (9.*c.cf11[l]-c.cf9[l])*c.delw[l]; 
                    c.psid4[l] = (c.cf11[l])*c.delw[l]; 
                } 
                
                //Calcul des vecteurs propres a gauche 
                c.vpr[0][0] = 1.; 
                c.vpr[1][0] = ur; 
                c.vpr[2][0] = vr;
                c.vpr[3][0] = wr-cr;
                c.vpr[4][0] = Hr-wr*cr; 
                
                c.vpr[0][1] = 1.; 
                c.vpr[1][1] = ur;
                c.vpr[2][1] = vr;
                c.vpr[3][1] = wr;
                c.vpr[4][1] = ur2/2. + vr2/2. + wr2/2.;
                
                c.vpr[0][2] = 0.; 
                c.vpr[1][2] = 0.;
                c.vpr[2][2] = cr;
                c.vpr[3][2] = 0.;
                c.vpr[4][2] = ur*cr;
                
                c.vpr[0][3] = 0.; 
                c.vpr[1][3] = 0.;
                c.vpr[2][3] = 0;
                c.vpr[3][3] = cr;
                c.vpr[4][3] = vr*cr;
                
                c.vpr[0][4] = 1.; 
                c.vpr[1][4] = ur;
                c.vpr[2][4] = vr;
                c.vpr[3][4] = wr+cr;
                c.vpr[4][4] = Hr+wr*cr; 
                
                //Calcul des corrections dans la base des vecteurs propres du systeme 
                for(int l=0;l<5;l++){ 
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
                for(int m=0;m<5;m++){ 
                    for(int l=0;l<5;l++){ 
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
                grille[i][j][k] = c; 
            }
        }
    }    //Fin de la boucle de calcul des valeurs initiales 
    
    
    //Boucle de calcul des indicateurs de monotonicite
    
    for(int l=0;l<5;l++){ 
        for(int i=1;i<Nx+2*marge-1;i++){
            for(int j=1;j<Ny+2*marge-1;j++){
                for(int k=1;k<Nz+2*marge-1;k++){ 
									Cellule c = grille[i][j][k]; 
									Cellule cg = grille[i][j][k-1]; 
                    c.am[l] = c.lambda[l]*c.delw[l]-cg.lambda[l]*cg.delw[l]; 
                    grille[i][j][k] = c; 
                }
            }
        } 
        //Calcul de dj^m4 
        for(int i=0;i<Nx+2*marge;i++){
            for(int j=0;j<Ny+2*marge;j++){   
                for(int k=1;k<Nz+2*marge-2;k++){ 
									Cellule c = grille[i][j][k]; 
									Cellule cd = grille[i][j][k+1]; 
                    double z1 = 4.*c.am[l]-cd.am[l]; 
                    double z2 = 4.*cd.am[l]-c.am[l]; 
                    double z3 = c.am[l]; 
                    double z4 = cd.am[l]; 
                    c.am1[l] = (sign(z1)+sign(z2))/2.*abs((sign(z1)+sign(z3))/2.)*(sign(z1)
                               + sign(z4))/2.*min(abs(z1),min(abs(z2),min(abs(z3),abs(z4)))); 
                    grille[i][j][k] = c; 
                }
            }
        } 
    } 
    
    //Boucle de calcul de r+ et r- 
    for(int l=0;l<5;l++){ 
        for(int i=marge;i<Nx+2*marge-4;i++){ 
            for(int j=marge;j<Ny+2*marge-4;j++){
                for(int k=marge;k<Nz+2*marge-4;k++){ 
									Cellule c = grille[i][j][k]; 
									Cellule cd = grille[i][j][k+1]; 
									Cellule  cg = grille[i][j][k-1]; 
                    c.rp[l] = sign(c.delw[l])*sign(cg.delw[l])*(abs(cg.delw[l])+eps)/(abs(c.delw[l])+eps); 
                    c.rm[l] = sign(c.delw[l])*sign(cd.delw[l])*(abs(cd.delw[l])+eps)/(abs(c.delw[l])+eps); 
                    //Corrections d'ordre superieur 
                    Cellule cg2 = grille[i][j][k-2]; 
										Cellule cg3 = grille[i][j][k-3]; 
										Cellule cg4 = grille[i][j][k-4]; 
										Cellule cg5 = grille[i][j][k-5]; 
										Cellule cd2 = grille[i][j][k+2]; 
										Cellule cd3 = grille[i][j][k+3]; 
										Cellule cd4 = grille[i][j][k+4]; 
                    c.psid[l] = -c.psid0[l]+cg.psid0[l]+cd.psid1[l]-cg2.psid1[l]-cd2.psid2[l]
                    + cg3.psid2[l]+cd3.psid3[l]-cg4.psid3[l]-cd4.psid4[l]+cg5.psid4[l]; 
                    grille[i][j][k] = c; 
                }
            }
        } 
    } 
    
    //Boucle de calcul des flux 
    for(int i=marge-1;i<Nx+marge;i++){ 
        for(int j=marge-1;j<Ny+marge;j++){ 
            for(int k=marge-1;k<Nz+marge;k++){ 
                //Cellule de reference 
                Cellule c = grille[i][j][k]; 
                //Cellules voisines 
                Cellule cg = grille[i][j][k-1]; 
								Cellule cg2 = grille[i][j][k-2]; 
								Cellule cg3 = grille[i][j][k-3]; 
								Cellule cg4 = grille[i][j][k-4]; 
								Cellule cd = grille[i][j][k+1]; 
								Cellule cd2 = grille[i][j][k+2]; 
								Cellule cd3 = grille[i][j][k+3]; 
								Cellule cd4 = grille[i][j][k+4]; 
                
                //Flux TVD 
                double tvd[5]; 
                double psict[5]; 
                
                //Initialisation 
                for(int l=0;l<5;l++){ 
                    tvd[l] = 0.; 
                    //Partie centre
                    psict[l] = c.psic0r[l] + cg.psic1r[l] + cd.psic1r[l] + cg2.psic2r[l] + cd2.psic2r[l] 
                               + cg3.psic3r[l] + cd3.psic3r[l] + cg4.psic4r[l] + cd4.psic4r[l]; 
                } 
                
                //Limiteur 
                double psic; 
                for(int l=0;l<5;l++){ 
                    psic = c.psic0[l] + cg.psic1[l] + cd.psic1[l] + cg2.psic2[l] 
                           + cd2.psic2[l] + cg3.psic3[l] + cd3.psic3[l] + cg4.psic4[l] + cd4.psic4[l]; 
                    
                    //Partie descentre
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
                    psi = (double) sign(c.delwnu[l])*psi/(abs(c.delwnu[l]+eps)); 
                    double psimax1 = 2.*r*(1.-xnume)/(xnu*(1.-xnu)); 
                    double psimax2 = 2./(1.-xnu); 
                    double psitvd = max(0.,min(psi,min(psimax1,psimax2))); 
                    
                    //Critere de monotonicite
                    if((c.delwnu[l] != 0.) && (abs(psi-psitvd)>eps)){ 
                        double dfo = psi*c.delwnu[l]/2.; 
                        double dabsf = psimax2*c.delwnu[l]/2.; 
                        double dful = psimax1*c.delwnu[l]/2.; 
                        double dfmd = dabsf/2.-c.am1[l]/2.; 
												Cellule camont = grille[i-is][j][k];   //Cellule en amont ddescentre
                        double dflc = dful/2.+((1.-xnume)/xnu)*camont.am1[l]/2.; 
                        double dfmin = max(min(0.,min(dabsf,dfmd)),min(0.,min(dful,dflc))); 
                        double dfmax = min(max(0.,max(dabsf,dfmd)),max(0.,max(dful,dflc))); 
                        if((dfmin-dfo)*(dfmax-dfo)>0.){ 
                            psi = psitvd; 
                        } 
                    } 
                    
                    //A decommenter si on veut utiliser uniquement le TVD et pas le MP 
                    //psi = psitvd; 
                    
                    //A decommenter si on veut utiliser le schema sans TVD ni MP 
                    //psi = 0.; 
                    
                    double ctvd = psi*c.delwnu[l]/2.-abs(c.lambda[l])*c.delw[l]/2.; 
                    for(int m=0;m<5;m++){ 
                        tvd[m] += ctvd*c.vpr[m][l];
                    } 
                } 
                //Fin du calcul de la correction TVD 
                
                // Calcul final du flux a droite de la cellule c 
                for(int l=0;l<5;l++){ 
                    c.fluxk[l] += tvd[l]; 
                }
                
                grille[i][j][k] = c; 
            }
        }
    } //Fin de la boucle sur les cellules
    
    
    // Correction d'entropie
     corentz(sigma);
    
    if(BC_z_in ==  1 || BC_z_out ==  1){
        for(int i=0;i<Nx+2*marge;i++){
            for(int j=0;j<Ny+2*marge;j++){
                
                if(BC_z_in ==  1){
                    Cellule c = grille[i][j][marge-1];
                    Cellule cp = grille[i][j][marge];
                    double p0 = (gam-1.)*(cp.rhoE0 - 1./2.*(cp.impx0*cp.impx0 + cp.impy0*cp.impy0 
                                + cp.impz0*cp.impz0)/cp.rho0);
                    //double p = c.p;
                    c.fluxk[0] = 0.;
                    c.fluxk[1] = 0.;
                    c.fluxk[2] = 0.;
                    c.fluxk[3] = p0;
                    c.fluxk[4] = 0.;
                    
                    grille[i][j][marge-1] = c;
                }
                
                if(BC_z_out ==  1){
                    Cellule c2 = grille[i][j][Nz+marge-1]; 
                    //Cellule cp = grille[Nx+marge-1][j];
                    double p02 = (gam-1.)*(c2.rhoE0-1./2.*(c2.impx0*c2.impx0 + c2.impy0*c2.impy0 
                                + c2.impz0*c2.impz0)/c2.rho0);
                    //double p2 = c2.p;
                    c2.fluxk[0] = 0.;
                    c2.fluxk[1] = 0.;
                    c2.fluxk[2] = 0.;
                    c2.fluxk[3] = p02;
                    c2.fluxk[4] = 0.;
                    
                    grille[i][j][Nz+marge-1] = c2;
                }
                
            } 
        }
    }   
}



void Grille::Solve(const double dt, double t, int n){
    
    //Cellule c;   
    for(int i=0;i<Nx+2*marge;i++){
        for(int j=0;j<Ny+2*marge;j++){
            for(int k=0;k<Nz+2*marge;k++){
							Cellule  c = grille[i][j][k];
                c.rho0 = c.rho;
                c.impx0 = c.impx;
                c.impy0 = c.impy;
                c.impz0 = c.impz;
                c.rhoE0 = c.rhoE;
                grille[i][j][k] = c;
            }
        }
    }
   
    //alternance directionnelle a chaque pas de temps
	
    if(n%6==0){
        fnumx(dt/dx,t);     //Calcul des flux numeriques pour le fluide seul selon x
        solve_fluidx(dt);  //Resolution du fluide : 1er demi-pas de temps selon x
		BC();           //Imposition des conditions aux limites pour le fluide
        
        fnumy(dt/dy,t);     //Calcul des flux numeriques pour le fluide seul selon y
        solve_fluidy(dt);  //Resolution du fluide : 1er demi-pas de temps selon y
		BC();           //Imposition des conditions aux limites pour le fluide
        
        fnumz(dt/dz,t);     //Calcul des flux numeriques pour le fluide seul selon z
        solve_fluidz(dt);  //Resolution du fluide : 1er demi-pas de temps selon z
        BC();           //Imposition des conditions aux limites pour le fluide
        
    } 
    else if(n%6==2){
        fnumx(dt/dx,t);     //Calcul des flux numeriques pour le fluide seul selon x
        solve_fluidx(dt);  //Resolution du fluide : 1er demi-pas de temps selon x
		BC();           //Imposition des conditions aux limites pour le fluide
        
        fnumz(dt/dz,t);     //Calcul des flux numeriques pour le fluide seul selon z
        solve_fluidz(dt);  //Resolution du fluide : 1er demi-pas de temps selon z
		BC();           //Imposition des conditions aux limites pour le fluide
        
        fnumy(dt/dy,t);     //Calcul des flux numeriques pour le fluide seul selon y
        solve_fluidy(dt);  //Resolution du fluide : 1er demi-pas de temps selon y
		BC();           //Imposition des conditions aux limites pour le fluide
        
    } 
    
    else if(n%6==1){
        
        fnumy(dt/dy,t);     //Calcul des flux numeriques pour le fluide seul selon y
        solve_fluidy(dt);  //Resolution du fluide : 1er demi-pas de temps selon y
		BC();           //Imposition des conditions aux limites pour le fluide
        
        fnumx(dt/dx,t);     //Calcul des flux numeriques pour le fluide seul selon x
        solve_fluidx(dt);  //Resolution du fluide : 1er demi-pas de temps selon x
		BC();           //Imposition des conditions aux limites pour le fluide
        
        fnumz(dt/dz,t);     //Calcul des flux numeriques pour le fluide seul selon z
        solve_fluidz(dt);  //Resolution du fluide : 1er demi-pas de temps selon z
        BC();           //Imposition des conditions aux limites pour le fluide
        
    } 
    
    else if(n%6==3){
        
        fnumy(dt/dy,t);     //Calcul des flux numeriques pour le fluide seul selon y
        solve_fluidy(dt);  //Resolution du fluide : 1er demi-pas de temps selon y
		BC();           //Imposition des conditions aux limites pour le fluide
        
        fnumz(dt/dz,t);     //Calcul des flux numeriques pour le fluide seul selon z
        solve_fluidz(dt);  //Resolution du fluide : 1er demi-pas de temps selon z
		BC();           //Imposition des conditions aux limites pour le fluide
        
        fnumx(dt/dx,t);     //Calcul des flux numeriques pour le fluide seul selon x
        solve_fluidx(dt);  //Resolution du fluide : 1er demi-pas de temps selon x
		BC();           //Imposition des conditions aux limites pour le fluide
        
    }
    else if(n%6==4){
        
        fnumz(dt/dz,t);     //Calcul des flux numeriques pour le fluide seul selon z
        solve_fluidz(dt);  //Resolution du fluide : 1er demi-pas de temps selon z
		BC();           //Imposition des conditions aux limites pour le fluide
        
        fnumx(dt/dx,t);     //Calcul des flux numeriques pour le fluide seul selon x
        solve_fluidx(dt);  //Resolution du fluide : 1er demi-pas de temps selon x
		BC();           //Imposition des conditions aux limites pour le fluide
        
        fnumy(dt/dy,t);     //Calcul des flux numeriques pour le fluide seul selon y
        solve_fluidy(dt);  //Resolution du fluide : 1er demi-pas de temps selon y 
        BC();           //Imposition des conditions aux limites pour le fluide
        
        
    }
    else if(n%6==5){
        
        fnumz(dt/dz,t);     //Calcul des flux numeriques pour le fluide seul selon z
        solve_fluidz(dt);  //Resolution du fluide : 1er demi-pas de temps selon z
		BC();           //Imposition des conditions aux limites pour le fluide
        
        fnumy(dt/dy,t);     //Calcul des flux numeriques pour le fluide seul selon y
        solve_fluidy(dt);  //Resolution du fluide : 1er demi-pas de temps selon y
		BC();           //Imposition des conditions aux limites pour le fluide
        
        fnumx(dt/dx,t);     //Calcul des flux numeriques pour le fluide seul selon x
        solve_fluidx(dt);  //Resolution du fluide : 1er demi-pas de temps selon x
		BC();           //Imposition des conditions aux limites pour le fluide
        
    }
  
    
}



void Grille::BC(){ 
    
    // Inner Boundary Condition for x
    //Cellule c, cm, cp, cb;
    for(int i=0;i<marge;i++){
        for(int j=0;j<Ny+2*marge;j++){
            for(int k=0;k<Nz+2*marge;k++){
							Cellule c = grille[i][j][k];
							Cellule  cm = grille[2*marge-i-1][j][k];   //Cellule miroir
							Cellule cp = grille[Nx+i][j][k];          //Cellule periodique
                // cb = grille[marge][j][k];         //Cellule du bord
                
                if(BC_x_in ==  1){
                    // BC reflecting ("miroir");
                    c.rho = cm.rho;
                    c.u   = -cm.u;
                    c.v = cm.v;
                    c.w = cm.w;
                    c.p = cm.p;
                    c.impx = -cm.impx;
                    c.impy = cm.impy;
                    c.impz = cm.impz;
                    c.rhoE = cm.rhoE;
                }
                if(BC_x_in ==  3){
                    // BC outflow("transmisibles");
                    c.rho = cm.rho;
                    c.u   = cm.u;
                    c.v = cm.v;
                    c.w = cm.w;
                    c.p = cm.p;
                    c.impx = cm.impx;
                    c.impy = cm.impy;
                    c.impz = cm.impz;
                    c.rhoE = cm.rhoE;
                }
                else if(BC_x_in ==  2){
                    //BC periodic("periodique")
                    c.rho = cp.rho;
                    c.u   = cp.u;
                    c.v = cp.v;
                    c.w = cp.w;
                    c.p   = cp.p;
                    c.impx = cp.impx;
                    c.impy = cp.impy;
                    c.impz = cp.impz;
                    c.rhoE= cp.rhoE;
                }
                grille[i][j][k] = c;
            }
        }
    }
    
    
    // Outer Boundary Condition for x
    for(int i=Nx+marge;i<Nx+2*marge;i++){
        for(int j=0;j<Ny+2*marge;j++){
            for(int k=0;k<Nz+2*marge;k++){
							Cellule  c = grille[i][j][k];
							Cellule  cm = grille[2*Nx+2*marge-i-1][j][k];  //Cellule miroir
							Cellule  cp = grille[i-Nx][j][k];              //Cellule periodique
                //cb = grille[Nx+marge-1][j][k];        //Cellule du bord
                
                if(BC_x_out == 1){
                    //BC reflecting ("miroir");
                    c.rho = cm.rho;
                    c.u = -cm.u;
                    c.v = cm.v;
                    c.w = cm.w;
                    c.p   = cm.p;
                    c.impx = -cm.impx;
                    c.impy = cm.impy;
                    c.impz = cm.impz;
                    c.rhoE= cm.rhoE;
                }
                if(BC_x_out == 3){
                    //BC outflow("transmisibles");
                    c.rho = cm.rho;
                    c.u = cm.u;
                    c.v = cm.v;
                    c.w = cm.w;
                    c.p   = cm.p;
                    c.impx = cm.impx;
                    c.impy = cm.impy;
                    c.impz = cm.impz;
                    c.rhoE= cm.rhoE;
                }
                else if(BC_x_out == 2){
                    //BC periodic("periodique")
                    c.rho = cp.rho;
                    c.u = cp.u;
                    c.v = cp.v;
                    c.w = cp.w;
                    c.p   = cp.p;
                    c.impx = cp.impx;
                    c.impy = cp.impy;
                    c.impz = cp.impz;
                    c.rhoE= cp.rhoE;
                }
                grille[i][j][k] = c; 
            }
        }
    }
    
    
    
    // Inner Boundary Condition for y
    
    for(int i=0;i<Nx+2*marge;i++){
        for(int j=0;j<marge;j++){
            for(int k=0;k<Nz+2*marge;k++){
                
							Cellule c = grille[i][j][k];
							Cellule cm = grille[i][2*marge-j-1][k];             //Cellule miroir
							Cellule cp = grille[i][Ny+j][k];                   //Cellule periodique
                //cb = grille[i][marge][k];               //Cellule sur le bord
                
                if(BC_y_in == 1){
                    //BC reflecting ("miroir");
                    c.rho = cm.rho;
                    c.u   = cm.u;
                    c.v = -cm.v;
                    c.w   = cm.w;
                    c.p   = cm.p;
                    c.impx = cm.impx;
                    c.impy = -cm.impy;
                    c.impz = cm.impz;
                    c.rhoE= cm.rhoE;
                }
                if(BC_y_in == 3){
                    //BC outflow("transmisibles");
                    c.rho = cm.rho;
                    c.u   = cm.u;
                    c.v = cm.v;
                    c.w   = cm.w;
                    c.p   = cm.p;
                    c.impx = cm.impx;
                    c.impy = cm.impy;
                    c.impz = cm.impz;
                    c.rhoE= cm.rhoE;
                }
                else if(BC_y_in == 2){
                    //BC periodic("periodique")
                    c.rho = cp.rho;
                    c.u   = cp.u;
                    c.v = cp.v;
                    c.w = cp.w;
                    c.p   = cp.p;
                    c.impx = cp.impx;
                    c.impy = cp.impy;
                    c.impz = cp.impz;
                    c.rhoE= cp.rhoE;
                }
                grille[i][j][k] = c;
            }
        }
    }
    
    
    
    // Outer Boundary Condition for y
    
    for(int i=0;i<Nx+2*marge;i++){
        for(int j=Ny+marge;j<Ny+2*marge;j++){
            for(int k=0;k<Nz+2*marge;k++){
							Cellule  c = grille[i][j][k];
							Cellule  cm = grille[i][2*Ny+2*marge-j-1][k];      //Cellule miroir
							Cellule  cp = grille[i][j-Ny][k];                  //Cellule periodique
                // cb = grille[i][Ny+marge-1][k];        //Cellule sur le bord
                
                if(BC_y_out == 1){
                    //BC reflecting ("miroir");
                    c.rho = cm.rho;
                    c.u   = cm.u;
                    c.v = -cm.v;
                    c.w = cm.w;
                    c.p   = cm.p;
                    c.impx = cm.impx;
                    c.impy = -cm.impy;
                    c.impz = cm.impz;
                    c.rhoE= cm.rhoE;
                }
                if(BC_y_out == 3){
                    //BC outflow("transmisibles");
                    c.rho = cm.rho;
                    c.u   = cm.u;
                    c.v = cm.v;
                    c.w = cm.w;
                    c.p   = cm.p;
                    c.impx = cm.impx;
                    c.impy = cm.impy;
                    c.impz = cm.impz;
                    c.rhoE= cm.rhoE;
                }
                else if(BC_y_out == 2){
                    //BC periodic("periodique")
                    c.rho = cp.rho;
                    c.u = cp.u;
                    c.v = cp.v;
                    c.w = cp.w;
                    c.p   = cp.p;
                    c.impx = cp.impx;
                    c.impy = cp.impy;
                    c.impz = cp.impz;
                    c.rhoE= cp.rhoE;
                }
                
                grille[i][j][k] = c;
            }
        }
    }
    
    // Inner Boundary Condition for z
    
    for(int i=0;i<Nx+2*marge;i++){
        for(int j=0;j<Ny+2*marge;j++){
            for(int k=0;k<marge;k++){
							Cellule  c = grille[i][j][k];
							Cellule  cm = grille[i][j][2*marge-k-1];             //Cellule miroir
                Cellule  cp = grille[i][j][Nz+k];                   //Cellule periodique
                //cb = grille[i][j][marge];               //Cellule sur le bord
                
                if(BC_z_in == 1){
                    //BC reflecting ("miroir");
                    c.rho = cm.rho;
                    c.u = cm.u;
                    c.v = cm.v;
                    c.w = -cm.w;
                    c.p = cm.p;
                    c.impx = cm.impx;
                    c.impy = cm.impy;
                    c.impz = -cm.impz;
                    c.rhoE= cm.rhoE;
                }
                if(BC_z_in == 3){
                    //BC outflow("transmisibles");
                    c.rho = cm.rho;
                    c.u = cm.u;
                    c.v = cm.v;
                    c.w = cm.w;
                    c.p = cm.p;
                    c.impx = cm.impx;
                    c.impy = cm.impy;
                    c.impz = cm.impz;
                    c.rhoE= cm.rhoE;
                }
                else if(BC_z_in == 2){
                    //BC periodic("periodique")
                    c.rho = cp.rho;
                    c.u = cp.u;
                    c.v = cp.v;
                    c.w = cp.w;
                    c.p = cp.p;
                    c.impx = cp.impx;
                    c.impy = cp.impy;
                    c.impz = cp.impz;
                    c.rhoE= cp.rhoE;
                }
                grille[i][j][k] = c;
            }
        }
    }
    
    
    
    // Outer Boundary Condition for z
    
    for(int i=0;i<Nx+2*marge;i++){
        for(int j=0;j<Ny+2*marge;j++){
            for(int k=Nz+marge;k<Nz+2*marge;k++){
							Cellule  c = grille[i][j][k];
							Cellule  cm = grille[i][j][2*Nz+2*marge-k-1];      //Cellule miroir
							Cellule cp = grille[i][j][k-Nz];                  //Cellule periodique
                // cb = grille[i][j][Nz+marge-1];        //Cellule sur le bord
                
                if(BC_z_out == 1){
                    //BC reflecting ("miroir");
                    c.rho = cm.rho;
                    c.u = cm.u;
                    c.v = cm.v;
                    c.w = -cm.w;
                    c.p = cm.p;
                    c.impx = cm.impx;
                    c.impy = cm.impy;
                    c.impz = -cm.impz;
                    c.rhoE= cm.rhoE;
                }
                if(BC_z_out == 3){
                    //BC outflow("transmisibles");
                    c.rho = cm.rho;
                    c.u = cm.u;
                    c.v = cm.v;
                    c.w = cm.w;
                    c.p = cm.p;
                    c.impx = cm.impx;
                    c.impy = cm.impy;
                    c.impz = cm.impz;
                    c.rhoE= cm.rhoE;
                }
                else if(BC_z_out == 2){
                    //BC periodic("periodique")
                    c.rho = cp.rho;
                    c.u = cp.u;
                    c.v = cp.v;
                    c.w = cp.w;
                    c.p   = cp.p;
                    c.impx = cp.impx;
                    c.impy = cp.impy;
                    c.impz = cp.impz;
                    c.rhoE= cp.rhoE;
                }
                
                grille[i][j][k] = c;
            }
        }
    }
    
    
} 



void Grille::impression(int n){
    
/*	
	//Impression du fichier .dat
	const char* fichier;
	{
		std::ostringstream oss;
		oss << "resultats/xt" << n << ".dat";
		string s = oss.str();
		fichier = s.c_str();
	}
	std::ofstream xt(fichier,std::ios::out | std::ios::trunc);
	//Initialisation de la sortie tecplot
	xt << "TITLE=grandeurs obtenues" << endl;
	xt << "VARIABLES=X,Y,Z,'rho','u','v','w','p'" << endl;
	xt << "ZONE I=" << Nx <<  " ,J=" << Ny << ", K=" << Nz << ",F=POINT" << endl;
	for(int i=marge; i<Nx+marge; i++){ 
		for(int j=marge; j<Ny+marge; j++){ 
			for(int k=marge; k<Nz+marge; k++){ 
				Cellule c = grille[i][j][k]; 
				if(abs(c.alpha-1.)>eps){ 
					xt << c.x << " " << c.y << " "<< c.z  << " " << c.rho << " " << c.u << " " << c.v << " "<< c.w << " " << c.p <<endl;	
				} 
				else {
					xt << c.x << " " << c.y << " "<< c.z  << " " << 0 << " " << 0 << " " << 0 << " "<< 0 << " " << 0 <<endl;	
				}
				
			}
		} 
	}*/
	
    //Impression du fichier vtk
    //const char* fluidevtk;
    //{
  std::ostringstream oss;
  oss << "resultats/fluide" << n << ".vtk";
  string s = oss.str();
  //cout << s << endl;
  const char* const fluidevtk = s.c_str();
  //}
    
    
    //Ouverture des flux en donne en ecriture
    std::ofstream vtk(fluidevtk,ios::out);
    if(vtk)
    {
        // cout <<"ouverture de xt.vtk reussie" << endl;
    } else {
        cout <<"ouverture de fluide" << n << ".vtk rate" << endl;
    }
    //Initialisation du fichier vtk
    vtk << "# vtk DataFile Version 3.0" << endl;
    vtk << "#Simulation Euler" << endl;
    vtk << "ASCII" << endl;
    vtk<<"\n";
    vtk << "DATASET UNSTRUCTURED_GRID" << endl;
    vtk << "POINTS " << (Nx+1)*(Ny+1) *(Nz+1)<< " DOUBLE" << endl;
    
    for(int i=0; i<Nx+1; i++){
        for(int j=0; j<Ny+1; j++){ 
            for(int k=0; k<Nz+1; k++){ 
                vtk << i*dx << " " << j*dy << " " <<k*dz << " "<< endl;
            }
        }
    }
    vtk<<"\n";

	//Calcul du nombre de cellules vraiment fluides
	int Nfluides = 0;
	for(int i=marge; i<Nx+marge; i++){
        for(int j=marge; j<Ny+marge; j++){ 
            for(int k=marge; k<Nz+marge; k++){
			  Cellule c = grille[i][j][k]; 
			  if(abs(c.alpha-1.)>eps){
                Nfluides++;
			  }
			}
        }
    }
	
    vtk << "CELLS " << Nfluides << " " << 9*Nfluides<< endl;
    
    for(int i=marge; i<Nx+marge; i++){
        for(int j=marge; j<Ny+marge; j++){ 
            for(int k=marge; k<Nz+marge; k++){
			  Cellule c = grille[i][j][k]; 
			  if(abs(c.alpha-1.)>eps){
                vtk << 8 << " " << (k-marge)+(j-marge)*(Nz+1)+(i-marge)*(Nz+1)*(Ny+1) << " " << ((k-marge)+1)+(j-marge)*(Nz+1)+(i-marge)*(Nz+1)*(Ny+1) << " " << ((k-marge)+1)+((j-marge)+1)*(Nz+1) + (i-marge)*(Nz+1)*(Ny+1) << " "<< (k-marge)+((j-marge)+1)*(Nz+1) + (i-marge)*(Nz+1)*(Ny+1)<< " " <<  (k-marge)+(j-marge)*(Nz+1)+((i-marge)+1)*(Nz+1)*(Ny+1)<< " "<< ((k-marge)+1)+(j-marge)*(Nz+1) + ((i-marge)+1)*(Nz+1)*(Ny+1) << " " <<  ((k-marge)+1)+((j-marge)+1)*(Nz+1)+((i-marge)+1)*(Nz+1)*(Ny+1)<< " " << (k-marge)+((j-marge)+1)*(Nz+1)+((i-marge)+1)*(Nz+1)*(Ny+1)<< endl;
			  }
			}
        }
    }
    vtk<<"\n";
    vtk << "CELL_TYPES " <<Nfluides<<endl;
    for(int k=0; k<Nfluides; k++){ 
        vtk<<12<<endl;
    }
    
    vtk<<"\n";
    vtk << "CELL_DATA " << Nfluides << endl;
    //Pression
    vtk << "SCALARS pression double 1" << endl;
    vtk << "LOOKUP_TABLE default" << endl;
    for(int i=marge; i<Nx+marge; i++){
        for(int j=marge; j<Ny+marge; j++){ 
            for(int k=marge; k<Nz+marge; k++){ 
							Cellule c = grille[i][j][k]; 
							if(abs(c.alpha-1.)>eps){ 
								vtk << grille[i][j][k].p << endl;
							} else {
							  //vtk << 0. << endl;
							}
				
            }
        }
    }
    
    vtk<<"\n";
    //Densite
    vtk << "SCALARS densite double 1" << endl;
    vtk << "LOOKUP_TABLE default" << endl;
    for(int i=marge; i<Nx+marge; i++){
        for(int j=marge; j<Ny+marge; j++){ 
            for(int k=marge; k<Nz+marge; k++){ 
							Cellule c = grille[i][j][k]; 
							if(abs(c.alpha-1.)>eps){ 
								vtk << grille[i][j][k].rho << endl;
							} else {
							  //vtk << 0. << endl;
							}
						
            }
        }
    }
    
    vtk<<"\n";
    //Vitesse u
    vtk << "SCALARS u double 1" << endl;
    vtk << "LOOKUP_TABLE default" << endl;
    for(int i=marge; i<Nx+marge; i++){
        for(int j=marge; j<Ny+marge; j++){ 
            for(int k=marge; k<Nz+marge; k++){ 
							Cellule c = grille[i][j][k]; 
							if(abs(c.alpha-1.)>eps){ 
								vtk << grille[i][j][k].u << endl;
							} else {
							  //vtk << 0. << endl;
							}
                
            }
        }
    }
    vtk<<"\n";
    //Vitesse v
    vtk << "SCALARS v double 1" << endl;
    vtk << "LOOKUP_TABLE default" << endl;
    for(int i=marge; i<Nx+marge; i++){
        for(int j=marge; j<Ny+marge; j++){ 
            for(int k=marge; k<Nz+marge; k++){ 
							Cellule c = grille[i][j][k]; 
							if(abs(c.alpha-1.)>eps){ 
								vtk << grille[i][j][k].v << endl;
							} else {
							  //vtk << 0. << endl;
							}
                
            }
        }
    }
    vtk<<"\n";
    //Vitesse w
    vtk << "SCALARS w double 1" << endl;
    vtk << "LOOKUP_TABLE default" << endl;
    for(int i=marge; i<Nx+marge; i++){
        for(int j=marge; j<Ny+marge; j++){ 
            for(int k=marge; k<Nz+marge; k++){ 
							Cellule c = grille[i][j][k]; 
							if(abs(c.alpha-1.)>eps){ 
								vtk << grille[i][j][k].w << endl;
							} else {
							  //vtk << 0. << endl;
							}
                
            }
        }
    }
    
    
    //Impression du fichier vtk
    //   const char* fluidevtk2;
    //   {
    //     std::ostringstream oss;
    //     oss << "resultats/Test" << n << ".vtk";
    //     string s = oss.str();
    //     //cout << s << endl;
    //     fluidevtk2 = s.c_str();
    //   }
    //   //Ouverture des flux en donne en ecriture
    //   std::ofstream vtk2(fluidevtk2,ios::out);
    //   if(vtk2)
    //     {
    //      // cout <<"ouverture de xt.vtk reussie" << endl;
    //     } else {
    //     cout <<"ouverture de xt" << n << ".vtk rate" << endl;
    //   }
    //   //Initialisation du fichier vtk
    //   vtk2 << "# vtk DataFile Version 3.0" << endl;
    //   vtk2 << "#Simulation Euler" << endl;
    //   vtk2 << "ASCII" << endl;
    //   vtk2<<"DATASET RECTILINEAR_GRID"<<endl;
    //   vtk2<<"DIMENSIONS"<<" "<< Nx <<" "<<Ny<<" "<<Nz<<endl;
    // 
    //   vtk2<<"X_COORDINATES"<< " "<< Nx<<" " <<"float"<<endl;
    //   for(int i=marge;i<Nx+marge;i++){
    //        vtk2<< x+dx/2.+(i-marge)*dx << " ";
    // 
    //   }
    // vtk2<<endl;
    // 
    // vtk2<<"Y_COORDINATES"<< " "<< Ny<<" "<<"float"<<endl;
    // 
    //     for(int j=marge;j<Ny+marge;j++){ 
    //        vtk2<<y+dy/2.+(j-marge)*dy<< " ";
    //   }
    // vtk2<<endl;
    // 
    // vtk2<<"Z_COORDINATES"<< " "<< Nz<<" "<<"float"<<endl;
    // 	for(int k=marge;k<Nz+marge;k++){ 
    //         vtk2<<z+dz/2.+(k-marge)*dz<< " ";
    //   }
    // vtk2<<endl;
    // 
    // vtk2<<"POINT_DATA"<<" "<< Nx * Ny *Nz <<endl;
    // vtk2<<"SCALARS Densite double"<<endl;
    // vtk2<<"LOOKUP_TABLE default"<<endl;
    // for(int k=marge; k<Nz+marge; k++){ 
    //     for(int j=marge; j<Ny+marge; j++){ 
    // 	for(int i=marge; i<Nx+marge; i++){
    // 
    // 		vtk2 << grille[i][j][k].rho << endl;
    // 
    //       }
    //     }
    //   }
    // 
    // vtk2<<"SCALARS Pression double"<<endl;
    // vtk2<<"LOOKUP_TABLE default"<<endl;
    // for(int k=marge; k<Nz+marge; k++){ 
    //     for(int j=marge; j<Ny+marge; j++){ 
    // 	for(int i=marge; i<Nx+marge; i++){
    // 
    // 		vtk2 << grille[i][j][k].p << endl;
    // 
    //       }
    //     }
    //   }
    // 
    // vtk2<<"SCALARS U double"<<endl;
    // vtk2<<"LOOKUP_TABLE default"<<endl;
    // for(int k=marge; k<Nz+marge; k++){ 
    //     for(int j=marge; j<Ny+marge; j++){ 
    // 	for(int i=marge; i<Nx+marge; i++){
    // 
    // 		vtk2 << grille[i][j][k].u << endl;
    // 
    //       }
    //     }
    //   }
    // 
    // vtk2<<"SCALARS V double"<<endl;
    // vtk2<<"LOOKUP_TABLE default"<<endl;
    // for(int k=marge; k<Nz+marge; k++){ 
    //     for(int j=marge; j<Ny+marge; j++){ 
    // 	for(int i=marge; i<Nx+marge; i++){
    // 
    // 		vtk2 << grille[i][j][k].v << endl;
    // 
    //       }
    //     }
    //   }
    // 
    // vtk2<<"SCALARS W double"<<endl;
    // vtk2<<"LOOKUP_TABLE default"<<endl;
    // for(int k=marge; k<Nz+marge; k++){ 
    //     for(int j=marge; j<Ny+marge; j++){ 
    // 	for(int i=marge; i<Nx+marge; i++){
    // 
    // 		vtk2 << grille[i][j][k].w << endl;
    // 
    //       }
    //     }
    //   }
    
    //   int nb_points = (Nx+1)*(Ny+1) *(Nz+1);
    //   int values = 5; 
    //   int n_cells = (Nx)*(Ny)*(Nz);
    //   char name[30] = "Simulation Euler"; 
    //   double tab_points[nb_points][N_dim];
    //   double tab_values[n_cells][values];
    //   int tab_connectivity[n_cells][8];
    //   
    //   int b=0;
    //    for(int k=0; k<Nz+1; k++){ 
    //    for(int j=0; j<Ny+1; j++){ 
    //     for(int i=0; i<Nx+1; i++){
    //      tab_points[b][0] = i*dx; 
    //      tab_points[b][1] = j*dy;
    //      tab_points[b][2] = k*dz;
    //      b++;
    // 	}
    //     }
    //   }
    //   
    //   cout<<"tab conectivity : " <<endl;
    //    int l=0;
    //    for(int k=0; k<Nz; k++){ 
    //     for(int j=0; j<Ny; j++){ 
    //      for(int i=0; i<Nx; i++){
    //      tab_connectivity[l][0] = i+j*(Nx+1)+k*(Nx+1)*(Ny+1)+1;
    // 	tab_connectivity[l][1] = (i+1)+j*(Nx+1)+k*(Nx+1)*(Ny+1)+1;
    // 	tab_connectivity[l][2] = (i+1)+(j+1)*(Nx+1) + k*(Nx+1)*(Ny+1)+1;
    // 	tab_connectivity[l][3] = i+(j+1)*(Nx+1) + k*(Nx+1)*(Ny+1)+1;
    // 	tab_connectivity[l][4] = i+j*(Nx+1)+(k+1)*(Nx+1)*(Ny+1)+1;
    // 	tab_connectivity[l][5] = (i+1)+j*(Nx+1) + (k+1)*(Nx+1)*(Ny+1) +1;
    // 	tab_connectivity[l][6] = (i+1)+(j+1)*(Nx+1)+(k+1)*(Nx+1)*(Ny+1)+1;
    // 	tab_connectivity[l][7] = i+(j+1)*(Nx+1)+(k+1)*(Nx+1)*(Ny+1)+1;
    // 	l++;
    // 	  
    //         }
    //       }      
    //   }
    //   int r=0;
    //   for(int k=marge; k<Nz+marge; k++){ 
    //     for(int j=marge; j<Ny+marge; j++){ 
    // 	for(int i=marge; i<Nx+marge; i++){
    // 
    // 		tab_values[r][0] = grille[i][j][k].rho;
    // 		r++;
    // 
    //       }
    //     }
    //   }
    // 
    //   //Impression du fichier inp
    //    const char* fluideinp;
    //    {
    //      std::ostringstream oss;
    //      oss << "resultats/fluide." << n<< ".inp";
    //      string s = oss.str();
    //      cout << s << endl;
    //      fluideinp = s.c_str();
    //    }
    //   
    //   std::ofstream out(fluideinp,ios::out);
    //     if(out)
    //       {
    //       // cout <<"ouverture de fluide"<< n << ".inp reussie" << endl;
    //       } else {
    //       cout <<"ouverture de fluide" << n << ".inp rate" << endl;
    //     }
    //     
    //     //Write the mesh info
    //      out << nb_points << " "
    //        << n_cells << " "
    //        << 1
    //        << " 0 0\n";
    //      
    //      
    //      // Write the coordinates
    //       
    //       
    //      int k=1; // 1-based node number for UCD
    //      
    //      for(int l=0; l<nb_points; l++){ 
    //          	 
    //         out << k++ << "\t";
    //      	out<< tab_points[l][0]<<" " << tab_points[l][1]<<" " << tab_points[l][2]<<"\n";
    //    
    //        }
    // 
    //         int m=1; // 1-based node number for UCD
    //         
    //         for(int l=0; l<n_cells; l++){ 
    //             	 
    //             out << m++  <<" "<< 0<<" "<<"hex" << "\t";
    //         	out<< tab_connectivity[l][0]<<" " << tab_connectivity[l][1]<<" " << tab_connectivity[l][2]<<" " << tab_connectivity[l][3]<<" " 
    //         	   << tab_connectivity[l][4]<<" " << tab_connectivity[l][5]<<" " << tab_connectivity[l][6]<<" " << tab_connectivity[l][7]<<"\n";
    // 
    //           }
    //           
    
    // out<<1<<""<< 1<<endl;
    // out<<"pressure,"<<endl;
    // 
    // int j=1; // 1-based node number for UCD
    //         
    //         for(int l=0; l<n_cells; l++){ 
    //        out << j++<< " "<< tab_values[l][0] <<"\n";
    //}
    
}
#endif
