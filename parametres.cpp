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

#include "parametres.hpp"
//!\file
//!\authors Laurent Monasse and Maria Adela Puscas

//Definition of the fluid initial conditions
//In case of restart of the simulation
double rhoi[Nx][Ny][Nz];
double ui[Nx][Ny][Nz];
double vi[Nx][Ny][Nz];
double wi[Nx][Ny][Nz];
double pi[Nx][Ny][Nz];

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
  double indexp[Np][3];
  for(int i=0;i<Np;i++){
    init>>indexp[i][0]>>indexp[i][1]>>indexp[i][2];
    if(!init){
      cout << "Problem " << i << endl;
      getchar();
    }
  }
  char cells[6];
  int N;
  int Ninfo;
  init>>cells>>N>>Ninfo;
  
  //Recovery of indices where alpha is not zero
  int index[N][3];
  int indice;
  double x,y,z;
  int trash;
  init >> trash >> indice >> trash >> trash >> trash >> trash >> trash >> trash >> trash;
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
  //Dump the type of cells
  for(int i=0;i<N+2;i++){
    getline(init,dump);
  }
	
  //Recover pressure
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  for(int l=0;l<N;l++){
    double p;
    init >> p;
    int i = index[l][0];
    int j = index[l][1];
    int k = index[l][2];
    pi[i][j][k] = p;
    //Between two points with alpha<1, fill with the same value
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
  }
  //Recovery of density
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  for(int l=0;l<N;l++){
    double rho;
    init >> rho;
    int i = index[l][0];
    int j = index[l][1];
    int k = index[l][2];
    rhoi[i][j][k] = rho;
    //Between two points with alpha<1, fill with the same value
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
  }
  //Recovery of u
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  for(int l=0;l<N;l++){
    double u;
    init >> u;
    int i = index[l][0];
    int j = index[l][1];
    int k = index[l][2];
    ui[i][j][k] = u;
    //Between two points with alpha<1, fill with the same value
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
  }
  //Recovery of v
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  for(int l=0;l<N;l++){
    double v;
    init >> v;
    int i = index[l][0];
    int j = index[l][1];
    int k = index[l][2];
    vi[i][j][k] = v;
    //Between two points with alpha<1, fill with the same value
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
  }
  //Recovery of w
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  getline(init,dump);
  for(int l=0;l<N;l++){
    double w;
    init >> w;
    int i = index[l][0];
    int j = index[l][1];
    int k = index[l][2];
    wi[i][j][k] = w;
    //Between two points with alpha<1, fill with the same value
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

double Rho(double x,double y, double z){ 
  double rho;
  //if(x<0.1 && y<0.2 && y>0.1 && z>0.083)  {rho = 1.1845467;}
  //{rho = 1.4;}
  //else {rho = 1.1845467;}
  rho = 1.4;
	
// if(x>0.2 && x<0.6 &&  y>0.41 && y<0.46){rho=0.;}
// else {rho=1.4;}
  //In case of restart
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
double U(double x,double y, double z){ 
  double u = 1.;
  //In case of restart
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
double V(double x,double y, double z){ 
  double v = 1.; 
  //In case of restart
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
double W(double x,double y, double z){ 
  double w = 1.; 
  //In case of restart
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
double P(double x,double y, double z, double dx, double dy, double dz){ 
	
  double p;
  //if(x<0.1 && y<0.2 && y>0.1 && z>0.083){p = 101325.;} 
  //{p = 1.;} 
  //else {p = 101325.;}
  p=1.;
  // if(x>0.2 && x<0.6 &&  y>0.41 && y<0.46){p=0.;}
// else {p=1.;}
  //In case of restart
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
