#include "molecule.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include "../Eigen/Dense"
#include "../Eigen/Eigenvalues"
#include "../Eigen/Core"

using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;


 
void Molecule::print_(){
 cout << "\nPrinting geometry\n" << endl;
 for (int i = 0; i < natoms ; i++){
    printf("%10.10f  %10.10f %10.10f %10.10f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);}

 cout << "\nPrinting Hessian Matrix\n" << endl;
 for (int i = 0; i < 3*natoms ; i++){
    for (int j = 0; j < 3*natoms ; j++){
         printf("%10.10f \t",H[i][j]);
    }
    printf("\n");
 }
 
 cout << "\nPrinting Mass weighted Hessian Matrix inside print\n" << endl;
 cout << mw_H << endl;

}


Molecule::Molecule(const char *file_1, const char *file_2)
{
   std:: ifstream is (file_1);
   assert(is.good());
 
   is >> natoms;
   zvals = new double [natoms]; 
   mass = new double [3];
   //Reading coordinates (in bohr) 
   geom = new double* [natoms];
   for(int i=0; i<natoms; i++)
       geom[i] = new double[3];


   for(unsigned int i=0; i<natoms; i++){
      is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];
    }
    //  if (zvals[i] == 8.0){
         mass[0] = 15.99491461957;
         mass[1] = 1.00782503223;
         mass[2] = mass[1];
    cout << mass[0] << mass[1] << mass[2] << endl;
   is.close();


   hessn_mat = new double* [natoms*natoms*3];
   for (int i = 0; i <= natoms*natoms*3; i++){
       hessn_mat[i] = new double [3]; }

   H = new double* [natoms*3];
   for (int i = 0; i < natoms*3; i++){
       H[i] = new double [natoms*3]; }

   mw_H = Matrix(natoms*3,natoms*3);

   std:: ifstream line (file_2);
   assert(line.good());
   line >> num_atoms;

   if (num_atoms == natoms){      
   cout << "num_atoms = " << num_atoms << endl;}

   for (int i = 0; i < natoms*natoms*3; i++){
          line >>  hessn_mat[i][0] >> hessn_mat[i][1] >> hessn_mat[i][2];
   }
   
    for (int i=0; i< 3*natoms*natoms ; i++ ){
       for (int j =0; j < 3; j++){
          H[i/3][3*(i%3)+j] = hessn_mat[i][j]; }
    }
 
    for (int i =0; i<3*natoms; i++){
        for (int j =0; j<3*natoms; j++){
        mw_H(i,j) = H[i][j];
        mw_H(i,j) /= sqrt(mass[i/3]);
        mw_H(i,j) /= sqrt(mass[j/3]);
        }
    }       
        Eigen::SelfAdjointEigenSolver<Matrix> solver(mw_H);
        Matrix evecs = solver.eigenvectors();
        Matrix evals = solver.eigenvalues();

        cout << "\nEigenvalues of mass weighted hessian \n" << evals << endl;



        double conv = 4.3597447222071e-14; //E_h (kg cm^2 s^-2)
               conv /= pow(5.29177210903e-9,2); //Bohr^2 to cm^2
               conv /= 1.66053906660e-27; //amu to kg
               conv = sqrt(conv);
               conv *= 3.33e-11;
               conv /= (2*acos(-1.0));
        cout << "Printing vibrational frequencies in cm^-1" << endl;
        cout << conv*sqrt(-evals(0)) << endl;
        cout << conv*sqrt(evals(1)) << endl;
        cout << conv*sqrt(evals(2)) << endl;
        cout << conv*sqrt(evals(3)) << endl;
        cout << conv*sqrt(evals(4)) << endl;
        cout << conv*sqrt(evals(5)) << endl;
        cout << conv*sqrt(evals(6)) << endl;
        cout << conv*sqrt(evals(7)) << endl;
        cout << conv*sqrt(evals(8)) << endl;

  
    line.close();


} 
