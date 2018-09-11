#include "Vector.hh"
#include "Matrix.hh"

namespace LinearAlgebra
{
// Find the min value of vector or the psoition. Algorythm scales as O(N)
double MinValue(const Vector<double> &vector);

double MinValueLocation(const Vector<double> &vector);

// Returns a vecot with all the elemts replaced with their absolute value
Vector<double> Absolute(const Vector<double> &vector);

// Return a vector which each element risen to the given power
Vector<double> Pow(const Vector<double> &vector, double power);

// Return a vector which the exponent of each element
Vector<double> Exp(const Vector<double> &vector);

// Returns the same vector with everything above index removed
Vector<double> TruncateLower(const Vector<double> &vector, int index);

// Returns the same vector with everything below index removed
Vector<double> TruncateUpper(const Vector<double> &vector, int index);

// Returns the transpose of the matrix.
Matrix<double> Transpose(const Matrix<double> &matrix);

// Returns the inverse of the matrix
Matrix<double> Inverse(const Matrix<double> &matrix);

// Return a lower triangula matrix that when multipled with the transpose, gives the
// matrix back. ( may only work with symetric matricies)
Matrix<double> CholeskyDecomposition(const Matrix<double> &matrix);

// Decomposes the matrix using the LU method
void LUPDecomposition(Matrix<double> &luMatrix, Vector<int> &rowPermutaions, int &rowChanges);

// Returns the eigenvalues and eigenvectors of the matrix. This method is only suitable for
// symetric matrix
void EigenSolver(const Matrix<double> &matrix, Matrix<double> &eigenVectors,
				 Vector<double> &eigenValues);

// Solves a generalised eiganvalue problem of the form Ax = Bex. Both A and B must be symetric
// matrices for this method to work
void GeneralisedEigenSolver(const Matrix<double> &matrixLeft, const Matrix<double> &matrixRight,
                            Matrix<double> &eigenVectors, Vector<double> &eigenValues);

}