#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#ifndef PARAMETRES_H
#define PARAMETRES_H

using std::string;

//Parametres pour le fluide
const double gam = 1.4;            //constante des gaz parfaits 
const double eps = 0.00000000001;      //constante proche de 0 pour le controle 
const double epsa = 0.5;           //fraction de cellule coupee
const int ordremax = 11;           //ordre maximal du schema
const int marge = 6;              //marge de cellules a appliquer au debut et a la fin du tableau des cellules  

const int N_dim=3;

const double X0 = 0;              //pozition de l'origine
const double Y0 = 0;
const double Z0 = 0;

const int Nx =14;                 //nombre de cellules de fluide selon x
const int Ny =14;                 //nombre de cellules de fluide selon y
const int Nz =14;                 //nombre de cellules de fluide selon z

const double domainex = 2.;            //Largeur du domaine fluide selon x
const double domainey = 2.;          //Largeur du domaine fluide selon y
const double domainez = 2.;          //Largeur du domaine fluide selon z

const double deltax = domainex/Nx;      //Pas d'espace pour le fluide selon x
const double deltay = domainey/Ny;       //Pas d'espace pour le fluide selon y
const double deltaz = domainez/Nz;       //Pas d'espace pour le fluide selon z

//Parametres solides
const double rhos = 1.; //Densite du solide 
const double nu = 0.; //Coefficient de Poisson du materiau
const double E = 10000.; //Module d'Young du materiau

//Parametres temporels
const double T = 2.;             //temps total de simulation
const double cfl = 0.5;            //valeur de la cfl fluide
const double cfls = 0.1;           //Valeur de la cfl solide
const int nimp = 1000;                //Nombre d'impressions
const double dtimp = T/nimp;        //Pas de temps entre deux impressions
const int Nmax = 1;           //nombre maximal d'iterations en temps

//Conditions aux limites
//Type de CL :  1 = reflecting ("miroir"); 2 = periodic("periodique"); 3= outflow("transmisibles");

const int BC_x_in =  3;                 // Inner Boundary Condition for x
const int BC_x_out = 3;                 // Outer Boundary Condition for x
const int BC_y_in =  3;                 // Inner Boundary Condition for y
const int BC_y_out = 3;                 // Outer Boundary Condition for y
const int BC_z_in =  3;                 // Inner Boundary Condition for z
const int BC_z_out = 3;                 // Outer Boundary Condition for z

double Rho(double x = 0.,double y = 0., double z = 0.);

double U(double x = 0.,double y = 0., double z = 0.);

double V(double x = 0.,double y = 0., double z = 0.); 

double W(double x = 0.,double y = 0., double z = 0.);

double P(double x = 0.,double y = 0., double z = 0., double dx = deltax, double dy = deltay, double dz = deltaz); 

#endif
