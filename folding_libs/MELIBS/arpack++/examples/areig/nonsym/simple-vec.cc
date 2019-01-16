/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE Simple.cc
   Simple example program that illustrates how to solve a real
   nonsymmetric standard eigenvalue problem in regular mode
   using the AREig function.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "lnmatrxc.h"
#include "areig.h"
#include <math.h>
#include "arcomp.h"
#include "arlsmat.h"
#include "arlnsmat.h"
#include "arlssym.h"
#include "arlgsym.h"
#include "arlsnsym.h"
#include "arlgnsym.h"
#include "arlscomp.h"
#include "arlgcomp.h"

//    double M[7][7] = 
//    {
//      { 16,  -3,   0,   0,   0,   0,   0 },
//      {-13,  16,  -3,   0,   0,   0,   0 },
//      {  0, -13,  16,  -3,   0,   0,   0 },
//      {  0,   0, -13,  16,  -3,   0,   0 },
//      {  0,   0,   0, -13,  16,  -3,   0 },
//      {  0,   0,   0,   0, -13,  16,  -3 },
//      {  0,   0,   0,   0,   0, -13,  16 }
//    };

void CSCMatrix(double M[], int n, int &nnz, double *&A, 
	       int *&irow, int *&pcol )
{
  int nEdgeSize = 0;
  int nVertexSize = n;
  int nNumIndex = 0;
  int nRowIndex = 0;
  int nColumnIndex = 0;

  for (int i = 0; i < n; i ++)
    for (int j = 0; j < n; j ++)
      {
	if (M[i*n+j] != 0)
	  nEdgeSize ++;
      } 
  
  nnz = nEdgeSize;
  irow = new int[nEdgeSize];
  pcol = new int[nVertexSize +1];
  A = new double [nEdgeSize];

  for (nColumnIndex = 0; nColumnIndex < n; nColumnIndex ++)
    {
      pcol[nColumnIndex] = nNumIndex;
    
      for (nRowIndex  = 0; nRowIndex < n; nRowIndex ++)
	{
	  if (M[nRowIndex*n+nColumnIndex] != 0)
	    {
	      A[nNumIndex] = M[nRowIndex*n + nColumnIndex];
	      irow[nNumIndex] = nRowIndex;
	      nNumIndex ++;
	    }
	} 
    }
  pcol[nColumnIndex] = nNumIndex;
}
main(int nargc, char *argv[])
{

  // Defining variables needed to store A in CSC format.

  int     n;     // Dimension of matrix.
  int     nnz;   // Number of nonzero elements in A.
  int*    irow;  // Row index of all nonzero elements of A.
  int*    pcol;  // Pointer to the beginning of each column (in irow and A).
  double* A;     // Nonzero elements of A.
  int nev;
  ifstream fin;
  
  // Creating a double precision matrix.
  char* Filename;
  if (nargc <3)
    {
      printf("Usage: %s [Datafile] [Dimension] [nev]",argv[0]);
      exit(0);
    }
  fin.open(argv[1]);
  n = atoi(argv[2]);
  nev = atoi(argv[3]);

  double *M = new double[n*n];
  for (int i = 0; i < n; i ++)
    for (int j = 0; j < n; j ++)
      {
	fin >> M[i*n+j] ;
      }
  fin.close();

//      double M[7][7] = 
//      {
//        {1,0,0,0,0,0,0},
//        {0,2,0,0,0,0,0},
//        {0,0,3,0,0,0,0},
//        {0,0,0,4,0,0,0},
//        {0,0,0,0,5,0,0},
//        {0,0,0,0,0,6,0},
//        {0,0,0,0,0,0,7}
//      };

//    double M[7][7] = 
//    {
//      { 16,  -3,   0,   0,   0,   0,   0 },
//      {-13,  16,  -3,   0,   0,   0,   0 },
//      {  0, -13,  16,  -3,   0,   0,   0 },
//      {  0,   0, -13,  16,  -3,   0,   0 },
//      {  0,   0,   0, -13,  16,  -3,   0 },
//      {  0,   0,   0,   0, -13,  16,  -3 },
//      {  0,   0,   0,   0,   0, -13,  16 }
//    };
  CSCMatrix(M, n, nnz, A, irow, pcol);
  //  StiffnessMatrix(n, 10.0, nnz, A, irow, pcol);

  cout <<"n:" << n << endl;
  cout <<"nnz:" << nnz << endl;
  cout << "pcol: " << endl;
  for (int i = 0; i <= n; i ++)
    { cout << " " << pcol[i];}
  cout << endl << "irow[]:\n";
  for (int i = 0; i < nnz; i ++)
    { cout << " " << irow[i] << endl;}
  cout << "A[]:\n";
  for (int i = 0; i < nnz; i ++)
    { cout << " " << A[i] << endl;}

  // Defining AREig output variables.

  int     nconv;                       // Number of converged eigenvalues.
  double* EigValR = new double[nev];   // Real part of the eigenvalues.
  double* EigValI = new double[nev];   // Imaginary part of the eigenvalues.
  double* EigVec  = new double[(nev+1)*n];  // Eigenvectors.

//    nconv = AREig(EigValR, EigValI, EigVec, n, nnz, A, irow, pcol,nev);

    char* which = "LM";
    int ncv = 0;
    double tol = 0.0;
    int maxit = 0;
    double* resid = NULL;
    bool AutoShift = true;

// Creating a matrix in ARPACK++ format.
    ARluNonSymMatrix<double, double> matrix(n, nnz, A, irow, pcol);

// Defining the eigenvalue problem.

    ARluNonSymStdEig<double> prob(nev, matrix, which, ncv, tol,
                                   maxit, resid, AutoShift);

// Finding eigenvalues and the corresponding eigenvectors.

    prob.FindEigenvectors();

    for (int i = 0; i < prob.ConvergedEigenvalues(); i ++)
      {
	cout << " Eigenvalue[" << (i+1) << "] = ";
	cout << prob.Eigenvalue(i) << endl;
	cout << "Eigenvector[" << (i+1) << "]: \n";
	for (int j = 0; j < prob.GetN(); j++){
	  cout << prob.Eigenvector(i,j) << endl;
	}
      }

    cout << endl;
  delete [] irow;
  delete [] pcol;
  delete [] M;
  delete [] A;
  delete [] EigValR;
  delete [] EigValI;
  delete [] EigVec;

} // main


//  int main()
//  {

//    // Defining variables needed to store A in CSC format.

//    int     n;     // Dimension of matrix.
//    int     nnz;   // Number of nonzero elements in A.
//    int*    irow;  // Row index of all nonzero elements of A.
//    int*    pcol;  // Pointer to the beginning of each column (in irow and A).
//    double* A;     // Nonzero elements of A.

//    // Creating a double precision matrix.

//    n = 10;
//  //    double M[4][4] = 
//  //    {
//  //      {-1, 0, 2, 5},
//  //      {0, 3, 0, 0},
//  //      {0, 0, -4, 0},
//  //      {1, 0, 0, 2}
//  //    };
//  //    {
//  //      {3, 0, 0, -2, 0},
//  //      {0, 1, 0, 4, 0},
//  //      {3, 1, 2, 4, 5},
//  //      {2, 0, 3, 5, 5},
//  //      {4, 2, 3, 5, 3}
//  //    };
//    double M[3][3] =
//    {
//      {2, 3, 7},
//      {0, -1, -2},
//      {0, 0, 2}
//    };

//    StiffnessMatrix(n, 10.0, nnz, A, irow, pcol);

//    //CSCMatrix(M, n, nnz, A, irow, pcol );
//     cout <<"n:" << n << endl;
//     cout <<"nnz:" << nnz << endl;
//    cout << "pcol: " << endl;
//    for (int i = 0; i <= n; i ++)
//      { cout << " " << pcol[i];}
//    cout << endl << "irow[]:\n";
//    for (int i = 0; i < nnz; i ++)
//      { cout << " " << irow[i] << endl;}
//    cout << endl << "A[]:\n";
//    for (int i = 0; i < nnz; i ++)
//      { cout << " " << A[i] << endl;}


//    // Defining AREig output variables.

//    int     nconv;                       // Number of converged eigenvalues.
//    double* EigValR = new double[201];   // Real part of the eigenvalues.
//    double* EigValI = new double[201];   // Imaginary part of the eigenvalues.
//    double* EigVec  = new double[1201];  // Eigenvectors.

//    // Finding the five eigenvalues with largest magnitude
//    // and the related eigenvectors.
//    double nzval[] = { -1.0, 1.0, 3.0, 2.0, -4.0, 5.0, 2.0};
//    int irow1[] = {1, 4, 2, 1, 3, 1, 4 };
//    int pcol1[] = {0, 2, 3, 5, 7};
//    nconv = AREig(EigValR, EigValI, EigVec, n, nnz, A, irow, pcol, 1);
//    //  nconv = AREig(EigValR, EigValI, EigVec, 4, 7, nzval, irow1, pcol1, 2);
//    //  nconv = AREig(EigValR, EigValI, 4, 7, nzval, irow1, pcol1, 2);
//    //ARluNonSymMatrix<double> A2(4,7, nzval, irow1, pcol1);

  
//    // Printing eigenvalues.

//    cout << "Eigenvalues:" << endl;
//    for (int i=0; i<nconv; i++) {
//      cout << "  lambda[" << (i+1) << "]: " << EigValR[i];
//      if (EigValI[i]>=0.0) {
//        cout << " + " << EigValI[i] << " I" << endl;
//      }
//      else {
//        cout << " - " << fabs(EigValI[i]) << " I" << endl;
//      }
//    }
//    cout << endl;
  
  

//  } // main


