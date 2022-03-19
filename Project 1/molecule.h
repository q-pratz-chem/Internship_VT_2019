#include <string>

using namespace std;

class Molecule
{
   public:
      int natom;
      int *zvals;
      double **geom; 
      double **r;     
      double **ex; 
      double **ey; 
      double **ez;
      double ***sin_theta;
      double *angle_outofplane;
      double *mass;
      double x_com ;
      double y_com ;
      double z_com ;

      void print_geom();
      void inter_dist();      
      void bond_angle();   
          double uvector(int dim, int p, int q);
          double phi(int i, int j, int k);
      void outofplane_angle();   
      void dihedral_angle();
      void com();
      void translate(double xcm, double ycm, double zcm);
      void moment_inertia();

      Molecule(const char *filename);
      ~Molecule();
};
