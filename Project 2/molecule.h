#include <string>
#include "../Eigen/Dense"
#include "../Eigen/Eigenvalues"
#include "../Eigen/Core"

using namespace std;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
class  Molecule

{
  public:
     double *zvals;
     int natoms;
     int num_atoms;
     double **geom;
     double *mass;
     double **H; //hessian_matrix
     double **hessn_mat; //hessian_matrix
     Matrix mw_H; //mass weighted Hessian matrix


     void print_();

     Molecule(const char *file_1, const char *file_2);
//     ~Molecule();


};
