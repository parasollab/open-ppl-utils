
/* main.c */

#include<iostream.h>
#include<stdlib.h> // for atoi

extern "C" void __utils__setup_potential(int *mode);
extern "C" void __utils__potential(int *dim, double *phi,
               double *psi, double *omega, double *potential);


main() {
  int mode = 0;
  __utils__setup_potential(&mode);
  int numofResidue = 28;
  double *phi = new double[numofResidue-1];
      double *psi = new double[numofResidue-1];
      double *omega = new double[numofResidue-1];

  double pot;
  __utils__potential(&numofResidue, phi, psi, omega, &pot);



}
