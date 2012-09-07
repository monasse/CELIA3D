#ifndef ABS_H
#define ABS_H

// inline double abs(double x){
//   double a = x;
//   if(a<0.){ a = -a;}
//   return a;
// }


//Definition de la fonction sign 
inline double sign(double z){ 
  double s = 0.; 
  if(z>0.){s = 1.;} else {s = -1.;} 
  return s; 
} 

#endif
