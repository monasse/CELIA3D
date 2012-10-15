#include "parametres.hpp"

double Rho(double x,double y, double z){ 
  double rho;
  if(x<0.16)
  	{rho = 1.4;} 
   else {rho = 1.4;}

  return rho; 
} 
double U(double x,double y, double z){ 
   double u;
  if(x<0.16)
	{u = 1.;} 
   else {u = 1.;}

  return u; 
 }
 double V(double x,double y, double z){ 
  double v = 0.; 

 return v; 
 } 
double W(double x,double y, double z){ 
   double w = 0.; 

 return w; 
 } 
 double P(double x,double y, double z, double dx, double dy, double dz){ 

  double p; 
  if(x<0.16)
	{p = 1.;}
 else {p = 1.;}

  return p; 
}
