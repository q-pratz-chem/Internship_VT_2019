#include <iostream>
#include "molecule.h"

using namespace std;

int main(int argc, char *argv[])
{
   
   Molecule mol("geom_ch3cho.txt");
 
   cout << "\nnatoms = " << mol.natom << endl;
   cout << "\nPrinting the geometry...\n\n";
   mol.print_geom();

   cout << "\nPrinting the inter-atomic distance...\n\n";
   mol.inter_dist();
   cout << endl;

   cout << "\nPrinting the bond angles...\n\n";
   mol.bond_angle();
   cout << endl;

   cout << "\nPrinting the out of plane  angles...\n\n";
   mol.outofplane_angle();
   cout << endl;

  cout << "\nPrinting the dihedral  angles...\n\n";
  mol.dihedral_angle();
  cout << endl;

  cout << "\nPrinting the center of mass...\n\n";
  mol.com();
  cout << endl;


  cout << "\nPrinting the moment of inertia...\n\n";
  mol.moment_inertia();
  cout << endl;
   return 0;

}
