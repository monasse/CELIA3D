#ifndef ABS_H
#define ABS_H

// inline double abs(double x){
//   double a = x;
//   if(a<0.){ a = -a;}
//   return a;
// }


//Definition de la fonction sign 
inline double sign(double x){ 

  double s = (x > 0) ? 1 : ((x < 0) ? -1 : 0);
	
  return s; 
} 

#endif
