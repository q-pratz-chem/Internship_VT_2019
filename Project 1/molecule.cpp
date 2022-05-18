#include "molecule.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include "../Eigen/Dense"
#include "../Eigen/Eigenvalues"
#include "../Eigen/Core"

using namespace std;


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
//typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

void Molecule::print_geom()
{
   for(int i=0; i < natom; i++){
   
   printf("%d %10.6f %10.6f %10.6f \n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
}
}

//Function for calculating interatomic distances r_ij
void Molecule::inter_dist()
{
   double power(double, int);  
   double x;
   double y;
   double z;
   double sum;
   for (unsigned int i=0; i<natom; i++){
 
      for( unsigned int j=0; j<natom; j++){
         x = (geom[j][0] -geom[i][0]);
         y = (geom[j][1] -geom[i][1]);
         z = (geom[j][2] -geom[i][2]);
         
         sum= power(x,2)+ power(y,2)+ power(z,2);
         r[i][j]= sqrt(sum);
      }
   }

   for(int i=0; i<natom; i++){
      for(int j=0; j<natom; j++){
         printf("%d %d %8.5f  ", i, j, r[i][j]);
      }
      printf("\n");
   }


}

double power(double a, int n)
{ 
   int i; 
   double val=a; 
 
   for(i=1; i<n; i++){
      val *=a;
   }
   return val;
}


//Function for calculation of vectors
double Molecule::uvector(int dim, int p, int q){
      return -(geom[p][dim] - geom[q][dim])/r[p][q] ;
}

double Molecule::phi( int i, int j, int k){
      return  acos(uvector(0,j,k)* uvector(0,j,i) + uvector(1,j,k)*uvector(1,j,i) + uvector(2,j,k)*uvector(2,j,i));
}

void Molecule::bond_angle()
{
    for(unsigned int i=0; i<natom; i++){
       for(unsigned int j=0; j<i; j++){
          for(unsigned int k=0; k<j; k++){
             if(r[k][j] < 4.0 && r[j][i] < 4.0){

    printf(" %2d %2d %2d  %10.6f \n", i,j,k, (phi(i,j,k)*(180.0/acos(-1.0))));
              }
           }
        }
     }
}


//Function for calculating outofplane angle theta_ikjl
void Molecule::outofplane_angle(){
    double exx;
    double eyy;
    double ezz;
    double theta;
    double eijk_x ; 
    double eijk_y ;
    double eijk_z ;
    double ejkl_x ;
    double ejkl_y ;
    double ejkl_z ;

    //Calculation for cross product
    for(unsigned int i=0; i<natom; i++){
       for(unsigned int k=0; k<natom; k++){
          for(unsigned int j=0; j<natom; j++){
             for (unsigned int l=0; l<j; l++){
                 sin_theta[j][k][l] = sin(phi(j,k,l)); 

                 ejkl_x = uvector(1, k,j)*uvector(2,k,l) - uvector(2,k,j)*uvector(1,k,l);     
                 ejkl_y = uvector(2, k,j)*uvector(0,k,l) - uvector(0,k,j)*uvector(2,k,l);     
                 ejkl_z = uvector(0, k,j)*uvector(1,k,l) - uvector(1,k,j)*uvector(0,k,l);     

                 exx = ejkl_x * uvector(0,k,i);
                 eyy = ejkl_y * uvector(1,k,i);
                 ezz = ejkl_z * uvector(2,k,i);
                 
                 theta = ((exx+eyy+ezz)/ sin_theta[j][k][l]);

                  if(theta< -1.0){
                     angle_outofplane[l]= asin(-1.0);
                     }
                  else if (theta > 1.0){
                     angle_outofplane[l]= asin(1.0);
                     }
                  else angle_outofplane[l]= asin(theta);

                  if (i!=j && i!=k && i!=l && j!=k  && j!=l && k!=l && r[i][k] < 4.0 && r[j][k] < 4.0 && r[j][l] < 4.0){
                  printf(" %d-%d-%d-%d  %10.6f \n", i,j,k,l, angle_outofplane[l]*(180.0/acos(-1.0)));  }   
             
           }
          }
       }
     }


}


void Molecule::dihedral_angle()
{
  double exx = 1.0;  
  double eyy = 1.0;
  double ezz = 1.0;
  double rslt_x = 0.0;
  double rslt_y = 0.0;
  double rslt_z = 0.0;
  double norm = 0.0;
  double sign ;
  double dot_prdct = 0.0;
  double tau = 1.0;
  double eijk_x = 0.0 ; 
  double eijk_y = 0.0 ;
  double eijk_z = 0.0 ;
  double ejkl_x = 0.0 ;
  double ejkl_y = 0.0 ;
  double ejkl_z = 0.0 ;



 for(unsigned int i=0; i<natom; i++){
     for(unsigned int j=0; j<i; j++){
        for(unsigned int k=0; k<j; k++){
           for (unsigned int l=0; l<k; l++){
             eijk_x = (uvector(1,i,j) * uvector(2,j,k) - uvector(2,i,j) * uvector(1,j,k)) ; 
             eijk_y = (uvector(2,i,j) * uvector(0,j,k) - uvector(0,i,j) * uvector(2,j,k)) ; 
             eijk_z = (uvector(0,i,j) * uvector(1,j,k) - uvector(1,i,j) * uvector(0,j,k)) ; 
   
             ejkl_x = (uvector(1,j,k) * uvector(2,k,l) - uvector(2,j,k) * uvector(1,k,l)) ; 
             ejkl_y = (uvector(2,j,k) * uvector(0,k,l) - uvector(0,j,k) * uvector(2,k,l)) ; 
             ejkl_z = (uvector(0,j,k) * uvector(1,k,l) - uvector(1,j,k) * uvector(0,k,l)) ; 
               
             exx = eijk_x * ejkl_x ; 
             eyy = eijk_y * ejkl_y ; 
             ezz = eijk_z * ejkl_z ; 
              
             sin_theta[j][k][l] = sin(phi(j,k,l)); 
             if ( sin_theta[i][j][k] !=0 && sin_theta[j][k][l] !=0) {
                tau = (exx+eyy+ezz)/ (sin_theta[i][j][k]*sin_theta[j][k][l]);}


              if(tau< -1.0){
                 tau = acos(-1.0); }
              else if (tau > 1.0){
                 tau = acos(1.0); }
              else tau = acos(tau);

             //Computing sign of the torsion
              rslt_x = eijk_y * ejkl_z - eijk_z * ejkl_y;
              rslt_y = eijk_z * ejkl_x - eijk_x * ejkl_z;
              rslt_z = eijk_x * ejkl_y - eijk_y * ejkl_x;

              norm = pow(rslt_x, 2) + pow(rslt_y, 2) + pow(rslt_z ,2);

              rslt_x /= norm;
              rslt_y /= norm;
              rslt_z /= norm;
              sign = 1.0;
              dot_prdct= rslt_x*uvector(0,j,k) + rslt_y*uvector(1,j,k) + rslt_z*uvector(2,j,k);

              if (dot_prdct < 0.0) {sign=-1.0;}
               tau  = tau*sign;

              if ( r[i][j] < 4.0 && r[j][k] < 4.0 && r[k][l] < 4.0){ 
               printf(" %d-%d-%d-%d  %10.6f \n", i,j,k,l, (tau*(180.0/acos(-1.0))));  }  
              }
           }
        }
     }
}

void Molecule::com(){
  double sum_mx=0;
  double sum_my=0;
  double sum_mz=0;
  double sum_m=0;

  for (unsigned int i=0; i<natom; i++)
  {
     sum_mx += mass[i]*geom[i][0];  
     sum_my += mass[i]*geom[i][1];
     sum_mz += mass[i]*geom[i][2];  
     sum_m += mass[i];
  }

    x_com = sum_mx/sum_m; 
    y_com = sum_my/sum_m; 
    z_com = sum_mz/sum_m; 

printf("X-component  %10.8f \n", x_com);
printf("Y-component  %10.8f \n", y_com);
printf("Z-component  %10.8f \n", z_com);


translate(x_com, y_com, z_com);

}

void Molecule::translate(double xcm, double ycm, double zcm){
   
   for (int i = 0; i < natom; i++){
        geom[i][0] = geom[i][0] - xcm ;
        geom[i][1] = geom[i][1] - ycm ;
        geom[i][2] = geom[i][2] - zcm ;

   }

}

void Molecule::moment_inertia()
{
Matrix MoI(3,3);
MoI(0,0) = 0.0 ; 
MoI(1,1) = 0.0 ;
MoI(2,2) = 0.0 ;
MoI(1,0) = 0.0 ;
MoI(2,0) = 0.0 ;
MoI(2,1) = 0.0 ;




   for (int i=0; i<natom; i++){
      MoI(0,0) += mass[i] * (pow(geom[i][1], 2) + pow(geom[i][2], 2));
      MoI(1,1) += mass[i] * (pow(geom[i][0], 2) + pow(geom[i][2], 2));
      MoI(2,2) += mass[i] * (pow(geom[i][0], 2) + pow(geom[i][1], 2));
      MoI(1,0) -= mass[i] * geom[i][1] * geom[i][0] ;
      MoI(2,0) -= mass[i] * geom[i][2] * geom[i][0] ;
      MoI(2,1) -= mass[i] * geom[i][2] * geom[i][1] ;
   }
    MoI(0,1) = MoI(1,0) ; 
    MoI(0,2) = MoI(2,0) ;
    MoI(1,2) = MoI(2,1) ;
 
    cout << "Moment of Inertia Tensor in  (amu bohr^2) \n" << MoI << endl;

    Eigen::SelfAdjointEigenSolver<Matrix> solver(MoI);
    Matrix evecs = solver.eigenvectors();
    Matrix evals = solver.eigenvalues();   

    cout << "\nPrinciple Moment of Inertia \n" << evals << endl;
  
    double conv_factr = 0.529177210903 * 0.529177210903; //amu to angstrom
    cout << "\nPrinciple Moment of Inertia in (amu A^2) \n" << evals*conv_factr << endl; 
   
    double factr_conv = 1.66053906660 * pow(10, -24); //amu to gm
           factr_conv *= pow(5.29177210903 * pow(10, -9),2); //amu to cm
    cout << "\nPrinciple Moment of Inertia in (gm cm^2) \n" << evals*factr_conv << endl; 


   //Classification of rotor
   if (evals(0) < 1e-4) 
     cout << "\nMolecule is linear\n" << endl;
   else if ((fabs(evals(0)-evals(1)) < 1e-4) &&( fabs(evals(1)-evals(2)) < 1e-4))
     cout << "\nMolecule is spherical top\n" << endl;
   else if ((fabs(evals(0)-evals(1)) < 1e-4) &&( fabs(evals(1)-evals(2)) > 1e-4))
     cout << "\nMolecule is oblate symmetric top\n" << endl;
   else if ((fabs(evals(0)-evals(1)) > 1e-4) &&( fabs(evals(1)-evals(2)) < 1e-4))
     cout << "\nMolecule is prolate symmetric top\n" << endl;
   else cout << "\nMolecule is  an asymmetric top\n" << endl;

   //Rotational constants
    double pi_ = acos(-1.0);
    double const_  = 6.626176e-34 ;//planck's constant in Kg m^2 s^-1
           const_ /= 8.0 * pi_ * pi_ ;
           const_ /= (1.66053906660e-27* 0.529177210903e-10 * 0.529177210903e-10);
    double rot_const = const_ * 1e-6;
    cout << "\nRotational constants in (MHz) are \n" << endl;
    cout << "A = " << rot_const/ evals(0) << "\n" << "B = " << rot_const / evals(1) << "\n" << "C = " << rot_const / evals(2) << endl;



}  


Molecule:: Molecule(const char *filename)
{ 

   std:: ifstream is (filename);
   assert(is.good());

   is >> natom;
   zvals = new int[natom]; 
   mass = new double[natom];
   
   //Memory allocation for geometry  
   geom = new double* [natom];
   for(int i=0; i<natom; i++)
      geom[i]= new double[3];
    
   //Reading the geom from file
   for(unsigned int i=0; i<natom; i++){
      is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];
      if (zvals[i] == 1)
         mass[i] =1.007825; 
      else if (zvals[i] == 6)
         mass[i] =12.000000;
      else if (zvals[i] == 8)
         mass[i] = 15.994914;
   }
  
   //Memory allocation for interatomic distance r[][]
   r = new double* [natom];
   for (int i=0; i<natom; i++){
      r[i] = new double[natom];
   }

   //Memory allocation for vectors of angle
  ex = new double* [natom];
  ey = new double* [natom];
  ez = new double* [natom];
  for (int i=0; i<natom; i++){
      ex[i] = new double [natom];
      ey[i] = new double [natom];
      ez[i] = new double [natom];
  } 

 
  //Array allocation for bond angle
  sin_theta = new double** [natom];
  for (int i= 0; i< natom; i++){
     sin_theta[i] = new double* [natom];
        for (int j= 0; j< natom; j++){
            sin_theta[i][j] = new double [natom];
        }
  }

 
  angle_outofplane = new double[natom];

  is.close();
}

Molecule:: ~Molecule()
{
   for(int i=0; i<natom; i++){
      delete[] geom[i];
      delete[] r[i];
      delete[] ex[i];
      delete[] ey[i];
      delete[] ez[i];
 //     delete[] moi[i];
   }
      delete[] geom; 
      delete[] r;
      delete[] zvals;
      delete[] mass;
   
      delete[] ex;
      delete[] ey;
      delete[] ez;
      delete[] angle_outofplane;
 //     delete[] moi;
}

