#include "parametres.hpp"

/*!
 *  \file parametres.cpp
 * \brief Initialisation des variables fluide: densit&eacute;, pression, vitesse
 */

// En cas de reprise definition des conditions initiales pour le fluide
double rhoi[Nx][Ny][Nz]; //!< En cas de reprise (\a rep=true) d&eacute;finition des conditions initiales pour le fluide: densite
double ui[Nx][Ny][Nz]; //!< En cas de reprise (\a rep=true) d&eacute;finition des conditions initiales pour le fluide:vitesse selon x
double vi[Nx][Ny][Nz]; //!< En cas de reprise (\a rep=true) d&eacute;finition des conditions initiales pour le fluide: vitesse selon y
double wi[Nx][Ny][Nz]; //!< En cas de reprise (\a rep=true) d&eacute;finition des conditions initiales pour le fluide: vitesse selon z
double pi[Nx][Ny][Nz]; //!< En cas de reprise (\a rep=true) d&eacute;finition des conditions initiales pour le fluide: pression

/**
 \fn void reprise()
 \brief Reprise d'une simulation (param&egrave;tre \a rep=true): initialisation des variables du probl&egrave;me &agrave; partir des fichier "resultats/fluide" << \a numrep << ".vtk" et "resultats/solide" << \a numrep << ".vtk". 
 \note Fonction utilis&eacute;e principalement pour le debugage. <b> Il  n'est pas n&eacute;cessaire de la re-coder! </b>
 \return void
 */
void reprise(){
  std::ostringstream oss;
  oss << "resultats/fluide" << numrep << ".vtk";
  string s = oss.str();
  const char* nom = s.c_str();
  std::ifstream init(nom,std::ios::in);
  string dump;
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  char points[7];
  int Np;
  char Double[7];
  init >> points >> Np >> Double;
  //cout << Np << endl;
  double indexp[Np][3];
  for(int i=0;i<Np;i++){
    init>>indexp[i][0]>>indexp[i][1]>>indexp[i][2];
    //cout << indexp[i][0] << indexp[i][1] << indexp[i][2] << endl;
    if(!init){
      cout << "probleme " << i << endl;
      getchar();
    }
  }
  //getline(init,dump);
  char cells[6];
  int N;
  int Ninfo;
  init>>cells>>N>>Ninfo;
  //cout << N << endl;
  //On enregistre le numero du point de coin de (i,j,k) inferieur pour chaque cellule
  //Utilise fondamentalement la construction de la sortie fluide dans impression()
  //Recuperation des indices ou alpha est non nul
  int index[N][3];
  int indice;
  double x,y,z;
  int trash;
  init >> trash >> indice >> trash >> trash >> trash >> trash >> trash >> trash >> trash;
  //cout << trash;
  x = indexp[indice][0];
  y = indexp[indice][1];
  z = indexp[indice][2];
  int l=0;
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      for(int k=0;k<Nz;k++){
	if(abs(i*deltax-x)<0.1*deltax && abs(j*deltay-y)<0.1*deltay && abs(k*deltaz-z)<0.1*deltaz){
	  index[l][0] = i;
	  index[l][1] = j;
	  index[l][2] = k;
	  l++;
	  if(l<N){
	    init >> trash >> indice >> trash >> trash >> trash >> trash >> trash >> trash >> trash;
	    x = indexp[indice][0];
	    y = indexp[indice][1];
	    z = indexp[indice][2];
	  }
	}
      }
    }
  }
  //On jette le type de cellules
  for(int i=0;i<N+2;i++){
    getline(init,dump);
  }
  
  //Recuperation de la pression
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  //cout << dump << endl;
  for(int l=0;l<N;l++){
    double p;
    init >> p;
    int i = index[l][0];
    int j = index[l][1];
    int k = index[l][2];
    pi[i][j][k] = p;
    //Entre deux points de alpha<1, on remplit avec la meme valeur
    int i1 = Nx-1;
    int j1 = Ny-1;
    int k1 = Nz-1;
    if(l<N-1){
      i1 = index[l+1][0];
      j1 = index[l+1][1];
      k1 = index[l+1][2];
    }
    if(k1+Nz*j1+Nz*Ny*i1-(k+Nz*j+Nz*Ny*i)>1){
      for(int it=k+Nz*j+Nz*Ny*i+1;it<k1+Nz*j1+Nz*Ny*i1;it++){
	int ktemp = it%Nz;
	int jtemp = ((it-ktemp)/Nz)%Ny;
	int itemp = (((it-ktemp)/Nz)-jtemp)/Ny;
	pi[itemp][jtemp][ktemp] = p;
      }
    }
    //cout << pi[i][j][k];
  }
  //Recuperation de la densite
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  //cout << dump << endl;
  for(int l=0;l<N;l++){
    double rho;
    init >> rho;
    int i = index[l][0];
    int j = index[l][1];
    int k = index[l][2];
    rhoi[i][j][k] = rho;
    //Entre deux points de alpha<1, on remplit avec la meme valeur
    int i1 = Nx-1;
    int j1 = Ny-1;
    int k1 = Nz-1;
    if(l<N-1){
      i1 = index[l+1][0];
      j1 = index[l+1][1];
      k1 = index[l+1][2];
    }
    if(k1+Nz*j1+Nz*Ny*i1-(k+Nz*j+Nz*Ny*i)>1){
      for(int it=k+Nz*j+Nz*Ny*i+1;it<k1+Nz*j1+Nz*Ny*i1;it++){
	int ktemp = it%Nz;
	int jtemp = ((it-ktemp)/Nz)%Ny;
	int itemp = (((it-ktemp)/Nz)-jtemp)/Ny;
	rhoi[itemp][jtemp][ktemp] = rho;
      }
    }
    //cout << rhoi[i][j][k];
  }
  //Recuperation de u
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  //cout << dump << endl;
  for(int l=0;l<N;l++){
    double u;
    init >> u;
    int i = index[l][0];
    int j = index[l][1];
    int k = index[l][2];
    ui[i][j][k] = u;
    //Entre deux points de alpha<1, on remplit avec la meme valeur
    int i1 = Nx-1;
    int j1 = Ny-1;
    int k1 = Nz-1;
    if(l<N-1){
      i1 = index[l+1][0];
      j1 = index[l+1][1];
      k1 = index[l+1][2];
    }
    if(k1+Nz*j1+Nz*Ny*i1-(k+Nz*j+Nz*Ny*i)>1){
      for(int it=k+Nz*j+Nz*Ny*i+1;it<k1+Nz*j1+Nz*Ny*i1;it++){
	int ktemp = it%Nz;
	int jtemp = ((it-ktemp)/Nz)%Ny;
	int itemp = (((it-ktemp)/Nz)-jtemp)/Ny;
	ui[itemp][jtemp][ktemp] = u;
      }
    }
    //cout << ui[i][j][k];
  }
  //Recuperation de v
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  //cout << dump << endl;
  for(int l=0;l<N;l++){
    double v;
    init >> v;
    int i = index[l][0];
    int j = index[l][1];
    int k = index[l][2];
    vi[i][j][k] = v;
    //Entre deux points de alpha<1, on remplit avec la meme valeur
    int i1 = Nx-1;
    int j1 = Ny-1;
    int k1 = Nz-1;
    if(l<N-1){
      i1 = index[l+1][0];
      j1 = index[l+1][1];
      k1 = index[l+1][2];
    }
    if(k1+Nz*j1+Nz*Ny*i1-(k+Nz*j+Nz*Ny*i)>1){
      for(int it=k+Nz*j+Nz*Ny*i+1;it<k1+Nz*j1+Nz*Ny*i1;it++){
	int ktemp = it%Nz;
	int jtemp = ((it-ktemp)/Nz)%Ny;
	int itemp = (((it-ktemp)/Nz)-jtemp)/Ny;
	vi[itemp][jtemp][ktemp] = v;
      }
    }
    //cout << vi[i][j][k];
  }
  //Recuperation de w
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  //cout << dump << endl;
  for(int l=0;l<N;l++){
    double w;
    init >> w;
    int i = index[l][0];
    int j = index[l][1];
    int k = index[l][2];
    wi[i][j][k] = w;
    //Entre deux points de alpha<1, on remplit avec la meme valeur
    int i1 = Nx-1;
    int j1 = Ny-1;
    int k1 = Nz-1;
    if(l<N-1){
      i1 = index[l+1][0];
      j1 = index[l+1][1];
      k1 = index[l+1][2];
    }
    if(k1+Nz*j1+Nz*Ny*i1-(k+Nz*j+Nz*Ny*i)>1){
      for(int it=k+Nz*j+Nz*Ny*i+1;it<k1+Nz*j1+Nz*Ny*i1;it++){
	int ktemp = it%Nz;
	int jtemp = ((it-ktemp)/Nz)%Ny;
	int itemp = (((it-ktemp)/Nz)-jtemp)/Ny;
	wi[itemp][jtemp][ktemp] = w;
      }
    }
  }
}


/**
 \fn double Rho(double x,double y, double z)
 \brief Initialisation du probl&egrave;me: valeur initiale de la densit&eacute; du fluide:
 * -rho
\param (x,y,z) position du centre de la cellule
 \return double
 */

double Rho(double x,double y, double z){ 
  double rho;
  if(x<1.5)  {rho = 1.4;}
  else {rho = 1.4;}
  //En cas de reprise
  if(rep){
    int i = (int) ((x)/deltax);
    int j = (int) ((y)/deltay);
    int k = (int) ((z)/deltaz);
    i = min(max(i,0),Nx-1);
    j = min(max(j,0),Ny-1);
    k = min(max(k,0),Nz-1);
    rho = rhoi[i][j][k];
  }
  return rho; 
} 

/**
 \fn double U(double x,double y, double z)
 \brief Initialisation du probl&egrave;me: valeur de la vitesse du fluide suivant x :
 * -u
 \param (x,y,z) position du centre de la cellule
 \return double
 */
double U(double x,double y, double z){ 
  double u = 1.;
  //En cas de reprise
  if(rep){
    int i = (int) ((x)/deltax);
    int j = (int) ((y)/deltay);
    int k = (int) ((z)/deltaz);
    i = min(max(i,0),Nx-1);
    j = min(max(j,0),Ny-1);
    k = min(max(k,0),Nz-1);
    u = ui[i][j][k];
  }
  return u; 
 }
/**
 \fn double V(double x,double y, double z)
 \brief Initialisation du probl&egrave;me: valeur de la vitesse du fluide suivant y: 
 * -v
 \param (x,y,z)  position du centre de la cellule
 \return double
 */
 double V(double x,double y, double z){ 
  double v = 0.; 
  //En cas de reprise
  if(rep){
    int i = (int) ((x)/deltax);
    int j = (int) ((y)/deltay);
    int k = (int) ((z)/deltaz);
    i = min(max(i,0),Nx-1);
    j = min(max(j,0),Ny-1);
    k = min(max(k,0),Nz-1);
    v = vi[i][j][k];
  }
  return v; 
 } 
/**
 \fn double W(double x,double y, double z)
 \brief Initialisation du probl&egrave;me: valeur de la vitesse du fluide suivant z :
 * -w
 \param (x,y,z)  position du centre de la cellule
 \return double
 */
double W(double x,double y, double z){ 
   double w = 0.; 
   //En cas de reprise
   if(rep){
    int i = (int) ((x)/deltax);
    int j = (int) ((y)/deltay);
    int k = (int) ((z)/deltaz);
    i = min(max(i,0),Nx-1);
    j = min(max(j,0),Ny-1);
    k = min(max(k,0),Nz-1);
    w = wi[i][j][k];
   }
   return w; 
 } 
/**
 \fn double P(double x,double y, double z,double dx, double dy, double dz)
 \brief Initialisation du probl&egrave;me: valeur de la pression du fluide:
 * -p
 \param (x,y,z) position du centre de la cellule 
 \param (dx, dy,dz) taille de la cellule
 \return double
 */
 double P(double x,double y, double z, double dx, double dy, double dz){ 

  double p;
  if(x<0.16){p = 1.;} 
  else {p = 1.;}
  //En cas de reprise
  if(rep){
    int i = (int) ((x)/deltax);
    int j = (int) ((y)/deltay);
    int k = (int) ((z)/deltaz);
    i = min(max(i,0),Nx-1);
    j = min(max(j,0),Ny-1);
    k = min(max(k,0),Nz-1);
    p = pi[i][j][k];
  }
  return p; 
}
