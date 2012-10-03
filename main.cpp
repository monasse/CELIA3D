#include <iostream>
#include <ctime>

//#include "fluide.hpp"
#include "fluide.cpp"
#include "couplage.cpp"
//#include "parametres.hpp"
#include "parametres.cpp"
#include "solide.cpp"


using namespace std;          // espace de nom standard


int main(){
  char tempsF[]="resultats/tempsF.dat";
	
	//Ouverture des flux en donne en ecriture
	std::ofstream out(tempsF,ios::out);
	if(out)
	{
		// cout <<"ouverture de xt.vtk reussie" << endl;
	} else {
		cout <<"ouverture de .dat rate" << endl;
	}
	
	Grille Fluide;
	
	
	double t=0.;

	Solide S;
	S.init("maillage.dat"); //Initialisation du solide a partir du fichier "maillage.dat"


	Fluide.init();
	Fluide.parois(S);
	Fluide.BC();
	//S.Affiche();
	

	int iter=0;
	
	clock_t start,end;
	
	start =clock();

	int kimp = 0; //Numero de suivi de l'impression
	double next_timp = dtimp; //Instant de la prochaine impression
	Fluide.impression(kimp);
	S.impression(kimp);
	kimp++;
	
	
	for (int n=0; (t<T) && n<Nmax; n++){
		
	  if(t>next_timp){
		Fluide.impression(kimp);
		S.impression(kimp);
		kimp++;
		next_timp += dtimp;
	  }
	  cout<<"Energie: "<< Fluide.Energie()+S.Energie() << " Solide:" << S.Energie() <<"  "<<"Masse : "<<"  "<< Fluide.Masse() <<endl;
	  double dt = min(Fluide.pas_temps(t, T),S.pas_temps(t,T));
		//Fluide.affiche("avant Solve");
		Fluide.Solve(dt, t, n);
		//Fluide.affiche("Solve");
		S.solve(dt);
		Fluide.modif_fnum(dt);
		//Fluide.affiche("modif_fnum");
		Fluide.mixage();
		//Fluide.affiche("mixage");
		Fluide.fill_cel(S);
		//Fluide.affiche("fill_cell");
		Fluide.BC();
		//Fluide.affiche("BC");
		out<< n << " temps actuel "<<t<<" dt "<<dt<<"\n";
		cout<<"iteration="<<n<< " dt="<<dt<<" t="<<t<<endl;
		t+= dt;
		iter++;
		
	}
	end=clock();
	
	Fluide.impression(kimp);
	S.impression(kimp);
	
	out<< "Temps final  "<< t<<endl;
	out<<"nb iter= "<< iter<<endl;    
	out <<"Temps de calcul " <<(double) (end-start)/CLOCKS_PER_SEC << endl;     
	cout<<"nb iter= "<< iter<<endl;    
	cout <<"Temps de calcul " <<(double) (end-start)/CLOCKS_PER_SEC << endl;    
	
	return 0;
	
	
}



// //test 
// 
// //Point_3 tr(Triangle_3 Tn, Triangle_3 Tn1, Point_3 Xn)
// 
// Point_3 A(1,0,0),B(2,0,0), C(1,1,0), E(3,0,0), F(4,0,0), G(3,1,0), Xn(0.,0.,0.), Z(1.,0.,0.), K(1.2,0.2,0.), R(1.,0.2,0.); 
// Triangle_3 Tn(A,B,C);
// Triangle_3 Tn1(E,F,G);
// Triangle_3 T(Z,K,R);
// 
// //Point_3 M = tr(Tn1,Tn,tr(Tn,Tn1,Xn));
// Triangle_3 Q = tr(Tn1,Tn,tr(Tn,Tn1,T));
// //cout<<"Point M :"<<M<<endl;
// //cout<<"Triangle Q :"<<Q<<endl;
// 
// Point_3 A1(3,0,1),B1(4,0,1), C1(3,1,1), X1(3.3,0.2,1.); 
// Triangle_3 Tri(A1,B1,C1);
// //cout<<"Point 3d :"<<X1<<endl;
// 
// //Point_2 tr(Triangle_3 Tn1, Point_3 Xn)
// Point_2 P2d = tr(Tri,X1);
// //cout<<"Point 2d :"<<P2d<<endl;
// 
// //Point_3 tr(Triangle_3 Tn1, Point_2 Xn){
// 	Point_3 P3d=tr(Tri,P2d);
// 	//cout<<"Point 3d :"<<P3d<<endl;
// 	
// 	Triangles tn, tn1;
// 	tn.push_back(T); 
// 	
// 	Point_3 N1(3.2,0.2,0.);
// 	Point_3 N2(3.4,0.2,0.);
// 	Point_3 N3(3.3,0.4,0.);
// 	tn1.push_back(Triangle_3(N1,N2,N3));
// 	
// 	
// 	
// 	Triangle_3 Tntest = tr(Tn,Tn1,T);
// 	cout<<"Triangle Test n->n+1 :"<<Tntest<<endl;
// 	
// 	Triangle_2 Intest=tr(Tn1,Tntest);
// 	cout<<"Triangle n->n+1 -> 3d- 2d:"<<Intest<<endl;
// 	
// 	//Triangle_2 tr(Triangle_3 Tn1, Triangle_3 T);
// 	Triangle_2 Test2d=tr(Tn1,tn1[0]);
// 	cout<<"Triangle 3d- 2d:"<<Test2d<<endl;
// 	
// 	Triangle_3 Test3d =tr(Tn1, Test2d);
// 	cout<<"Triangle 3d- 2d -3d :"<<Test3d<<endl;
// 	
// 	//void sous_maillage_faceTn_faceTn1(Triangle_3 Tn, Triangles tn, Triangle_3 Tn1, Triangles tn1,Triangles& T3d_n,Triangles& T3d_n1)
// 	Triangles Qn,Qn1;
// 	
// 	sous_maillage_faceTn_faceTn1(Tn,tn, Tn1,tn1,Qn,Qn1);
// 	
// 	for(int i=0; i<Qn.size(); i++){ 
// 		
// 		cout<<"triangle 3d n "<<Qn[i]<<endl;
// 	}		
// 	for(int i=0; i<Qn1.size(); i++){ 
// 		
// 		cout<<"triangle 3d n+1 "<<Qn1[i]<<endl;
// 	}
// 	//fin test 

