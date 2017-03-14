#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sstream>
#ifndef PARAMETRES_H
#define PARAMETRES_H

/*!
 \file
\authors Laurent Monasse and Maria Adela Puscas
\brief Parameters of the problem. The parameters specific to coupling are preceded by a "warning" sign
 */

using std::string;

//! \brief Flag for a potential restart from a previous recovery file: false if the simulation is not restarted, true if it is. \a rep indicates whether initial conditions (false) 
//! are used or a certain time output is used to restart the simulation (true). In the latter case, the index of the recovery file is given by \a numrep. 
int rep = false;
int numrep = 18;  //!<  Index of the recovery file


/*! 
 * \warning  <b> Specific coupling parameter ! </b>
 */
//! \brief  Flag for the type of scheme used: true for an explicit coupling and false for a semi-implicit coupling. 

bool explicite = true;
/*! 
 * \warning  <b> Specific coupling parameter ! </b>
 */

bool flag_2d = false;
bool couplage1d=false;
const bool exact_swap=false;

//Fluid parameters
const double gam = 1.4;                   //!<Perfect gas constant 
const double eps =  0.00000000000001;     //!<Numerical stabilization 
const double epsa = 0.5;                  //!<Limit of the size of small cut-cells
const int ordremax = 11;                  //!<Maximal order of the scheme
const int marge = 6;                      //!<Margin of cells on each side of the grid  
const double eps_vide =  0.0000000001;     
const int N_dim=3;

const double X0 = 0;              //!<Position of the origin of the grid
const double Y0 = 0;
const double Z0 = 0;

const int Nx =18;                 //!<Number of fluid cells in x
const int Ny =18;                 //!<Number of fluid cells in y
const int Nz =18;                 //!<Number of fluid cells in z

const double domainex = 0.3;            //!<Fluid domain size in x
const double domainey = 0.3;          //!<Fluid domain size in y
const double domainez = 0.3;          //!<Fluid domain size in z

const double deltax = domainex/Nx;      //!<Spatial discretization step in x
const double deltay = domainey/Ny;       //!<Spatial discretization step in y
const double deltaz = domainez/Nz;       //!<Spatial discretization step in z

//Solid parameters
const double rhos = 2698.9; //!<Solid density 
const double nu = 0.; //!<Poisson's ratio
const double E = 5.; //!<Young's modulus
const double k_max = 0.01; //!<Elongation at break

//Parametres temprels
const double T = 0.1;             //!<Total simulation time
const double cfl = 0.5;            //!<Fluid CFL condition
const double cfls = 0.5;           //!<Solid CFL condition
const int nimp = 10;                //!<Number of outputs
const double dtimp = T/nimp;        //!<Time-step between two consecutive outputs
const int Nmax = 1000000;           //!<Maximal number of time iterations

//!Boundary conditions
//!Types of BC:  1 = reflecting; 2 = periodic; 3= outflow; 

const int BC_x_in =  2;                 //!< Inner Boundary Condition for x
const int BC_x_out = 2;                 //!< Outer Boundary Condition for x
const int BC_y_in =  2;                 //!< Inner Boundary Condition for y
const int BC_y_out = 2;                 //!< Outer Boundary Condition for y
const int BC_z_in =  2;                 //!< Inner Boundary Condition for z
const int BC_z_out = 2;                 //!< Outer Boundary Condition for z

double Rho(double x = 0.,double y = 0., double z = 0.);

double U(double x = 0.,double y = 0., double z = 0.);

double V(double x = 0.,double y = 0., double z = 0.); 

double W(double x = 0.,double y = 0., double z = 0.);

double P(double x = 0.,double y = 0., double z = 0., double dx = deltax, double dy = deltay, double dz = deltaz); 

#endif
