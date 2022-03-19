#include <iostream>
#include "molecule.h"
#include <fstream>
#include <cstdio>

using namespace std;

int main(int argc, char *argv[])
{
  Molecule mol("geom_water.dat", "hessian_mat.dat");

  cout << "Printing the geometry" << endl;
  mol.print_();


 
return 0;

}
