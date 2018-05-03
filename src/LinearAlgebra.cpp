#include <cmath>
#include <limits>
#include <Eigenvalues>

#include "LinearAlgebra.hh"
#include "Vector.hh"
#include "Matrix.hh"

double LinearAlgebra::MinValue(const Vector<double> &vector)
{
	double small = vector(0);
	for (int i = 1; i < vector.Length(); ++i)
	{
		if (vector(i) < small) {small = vector(i);}
	}
	return small;
}

double LinearAlgebra::MinValueLocation(const Vector<double> &vector)
{
	double small = vector(0);
	int minPosition = 0;
	for (int i = 1; i < vector.Length(); ++i)
	{
		if (vector(i) < small)
		{
			small = vector(i);
			minPosition = i;
		}
	}
	return minPosition;
}

Vector<double> LinearAlgebra::Absolute(const Vector<double> &vector)
{
	Vector<double> absoulteVector(vector);
	for (int i = 0; i < absoulteVector.Length(); i++)
	{
		if (absoulteVector[i] < 0) { absoulteVector[i] *= -1; }
	}
	return absoulteVector;
}

Vector<double> LinearAlgebra::Pow(const Vector<double> &vector, double power)
{
	Vector<double> powerVector(vector.Length());
	for (int i = 0; i < powerVector.Length(); i++)
	{
		powerVector[i] =  std::pow(vector(i), power);
	}
	return powerVector;
}


Vector<double> LinearAlgebra::Exp(const Vector<double> &vector)
{
	Vector<double> powerVector(vector.Length());
	for (int i = 0; i < powerVector.Length(); i++)
	{
		powerVector[i] =  std::exp(vector(i));
	}
	return powerVector;
}

Vector<double> LinearAlgebra::TruncateLower(const Vector<double> &vector, int index)
{
	Vector<double> truncated(index);
	for (int i = 0; i < index ; i++)
	{
		truncated[i] = vector(i);
	}
	return truncated;
}

Vector<double> LinearAlgebra::TruncateUpper(const Vector<double> &vector, int index)
{
	Vector<double> truncated(vector.Length() - index);
	for (int i = 0; i < vector.Length() - index; i++)
	{
		truncated[i] = vector(i + index);
	}
	return truncated;
}

Matrix<double> LinearAlgebra::Transpose(const Matrix<double> &matrix)
{
	Matrix<double> transpose(matrix.GetRows(), matrix.GetColumns());
	for (int i = 0; i < matrix.GetRows(); i++)
	{
		for (int j = 0; j < matrix.GetColumns(); j++)
		{
			transpose[i][j] = matrix[j][i];
		}
	}
	return transpose;
}

Matrix<double> LinearAlgebra::Inverse(const Matrix<double> &matrix)
{
	int numberRows = matrix.GetRows();
	int numberColumns = matrix.GetColumns();
	Matrix<double> luMatrix = matrix.DeepCopy();


	Matrix<double> inverse(numberRows, numberColumns);
	Vector<int> rowPermutaions;
	int rowChanges;

	LUPDecomposition(luMatrix, rowPermutaions, rowChanges);

	for (int j = 0; j < numberColumns; j++)
	{
		for (int i = 0; i < numberRows; i++)
		{
			if (rowPermutaions[i] == j)
			{
				inverse[i][j] = 1.0;
			} else
			{
				inverse[i][j] = 0.0;
			}
			for (int k = 0; k < i; k++)
			{
				inverse[i][j] -= luMatrix[i][k] * inverse[k][j];
			}
		}

		for (int i = numberRows - 1; i >= 0; i--)
		{
			for (int k = i + 1; k < numberRows; k++)
			{
				inverse[i][j] -= luMatrix[i][k] * inverse[k][j];
			}
			inverse[i][j] = inverse[i][j] / luMatrix[i][i];
		}
	}

	return inverse;
}

Matrix<double> LinearAlgebra::CholeskyDecomposition(const Matrix<double> &matrix)
{
	int numberRows = matrix.GetRows();
	int numberColumns = matrix.GetColumns();
	Matrix<double> lowerMatrix(numberRows, numberColumns);

	for (int i = 0; i < numberRows; i++)
	{
		for (int k = 0; k < i + 1; k++)
		{

			if (k == i)
			{
				double sum(0);
				for (int j = 0; j < k; j++)
				{
					sum += lowerMatrix[k][j] * lowerMatrix[k][j];
				}
				lowerMatrix[k][k] = std::sqrt(matrix[k][k] - sum);
			} else
			{
				double sum(0);
				for (int j = 0; j < k; j++)
				{
					sum += lowerMatrix[i][j] * lowerMatrix[k][j];
				}
				lowerMatrix[i][k] = (1 / (lowerMatrix[k][k])) * (matrix[i][k] - sum);
			}
		}
	}
	return lowerMatrix;
}


void LinearAlgebra::LUPDecomposition(Matrix<double> &luMatrix, Vector<int> &rowPermutaions, int &rowChanges)
{
	rowChanges = 1;
	int numberRows = luMatrix.GetRows();
	int numberColumns = luMatrix.GetColumns();
	rowPermutaions = Vector<int>(numberRows);

	for (int i = 0; i < numberRows; i++)
	{
		rowPermutaions[i] = i;
	}

	for (int i = 0; i < numberRows; i++)
	{
		double maxA = 0.0;
		int imax = i;

		for (int k = i; k < numberColumns; k++)
		{
			if (std::abs(luMatrix[k][i]) > maxA)
			{
				maxA = std::abs(luMatrix[k][i]);
				imax = k;
			}
		}

		if (imax != i)
		{
			int j = rowPermutaions[i];
			rowPermutaions[i] = rowPermutaions[imax];
			rowPermutaions[imax] = j;

			for (int k = 0; k < numberColumns; k++)
			{
				double dummy = luMatrix[i][k];
				luMatrix[i][k] = luMatrix[imax][k];
				luMatrix[imax][k] = dummy;
			}
			rowChanges = -rowChanges;
		}

		for (int j = i + 1; j < numberColumns; j++)
		{
			luMatrix[j][i] /= luMatrix[i][i];
			for (int k = i + 1; k < numberColumns; k++)
			{
				luMatrix[j][k] -= luMatrix[j][i] * luMatrix[i][k];
			}
		}
	}
}

void LinearAlgebra::EigenSolver(const Matrix<double> &matrix, Matrix<double> &eigenVectors,
                                 Vector<double> &eigenValues)
{
	// Need to first get data from Matrix class and map to an Eigen matrix class
	Eigen::MatrixXd eMatrix(matrix.GetRows(), matrix.GetColumns());
	for (int i = 0; i < matrix.GetColumns(); i++)
	{
		eMatrix.col(i) = Eigen::Map<Eigen::VectorXd> (matrix[i], matrix.GetRows());
	}

	// Now set up the solver and find eigen vector/values
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver;
	eigenSolver.compute(eMatrix);

	// Now transfer data overt to eigenVector and eiganValues
	for (int i = 0; i < eigenVectors.GetRows(); i++)
	{
		eigenValues[i] = eigenSolver.eigenvalues()(i);
		for (int j = 0; j < eigenVectors.GetColumns(); j++)
		{
			eigenVectors[i][j] = eigenSolver.eigenvectors()(i,j);
		}
	}
}

void LinearAlgebra::GeneralisedEigenSolver(const Matrix<double> &matrixLeft, const Matrix<double> &matrixRight,
                            Matrix<double> &eigenVectors, Vector<double> &eigenValues)
{
	// Need to first get data from Matrix class and map to an Eigen matrix class
	Eigen::MatrixXd eMatrixL(matrixLeft.GetRows(), matrixLeft.GetColumns());
	Eigen::MatrixXd eMatrixR(matrixRight.GetRows(), matrixRight.GetColumns());
	for (int i = 0; i < matrixLeft.GetColumns(); i++)
	{
		eMatrixL.col(i) = Eigen::Map<Eigen::VectorXd> (matrixLeft[i], matrixLeft.GetRows());
		eMatrixR.col(i) = Eigen::Map<Eigen::VectorXd> (matrixRight[i], matrixRight.GetRows());
	}

	// Now set up the solver and find eigen vector/values
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver;
	eigenSolver.compute(eMatrixL, eMatrixR);

	// Now transfer data overt to eigenVector and eiganValues
	for (int i = 0; i < eigenVectors.GetColumns(); i++)
	{
		eigenValues[i] = eigenSolver.eigenvalues()(i);	
		for (int j = 0; j < eigenVectors.GetRows(); j++)
		{
			eigenVectors[j][i] = eigenSolver.eigenvectors()(j,i);
		}
	}
}
/*

void LinearAlgebra::GeneralisedEigenSolver(const Matrix<double> &matrixLeft, const Matrix<double> &matrixRight,
        Matrix<double> &eigenVectors, Vector<double> &eigenvalues)
{
	// First we need to decompose the right matrix to a lower and upper traingula
	Matrix<double> lowerMatrixRight = CholeskyDecomposition(matrixRight);
	// Now we form the new symetric matrix which has the same eiganvalues
	Matrix<double> symetricMatrix = Inverse(lowerMatrixRight) * matrixLeft * Inverse(Transpose(lowerMatrixRight));
	// Now find the eiganvalues and vectors of the matrix.
	EigenSolver(symetricMatrix, eigenVectors, eigenvalues);
	// Now find correct eiganvalues
	std::cout << "first" << std::endl;
	eigenVectors.Print();
	eigenVectors = Inverse(Transpose(lowerMatrixRight)) * eigenVectors;
	std::cout << "Second" << std::endl;
	eigenVectors.Print();

}



void LinearAlgebra::EigenSolver(const Matrix<double> &matrix, Matrix<double> &eigenVectors, Vector<double> &eigenValues,
                                int nrot, int maxSweeps)
{
	Vector<double> b(matrix.GetRows());
	Vector<double> z(matrix.GetRows());
	Matrix<double> A = matrix.DeepCopy();

	double small = std::numeric_limits<double>::epsilon();

	for (int i = 0; i < A.GetRows(); i++)
	{
		for (int j = 0; j < A.GetColumns(); j++)
		{
			eigenVectors[i][j] = 0.0;
			eigenVectors[i][i] = 1.0;
		}
	}
	for (int i = 0; i < A.GetRows(); i++)
	{
		b[i] = A[i][i];
		eigenValues[i] = A[i][i];
	}

	for (int i = 1; i < maxSweeps; i++)
	{
		double sum(0.0);
		for (int j = 0; j < A.GetRows() - 1; j++)
		{
			for (int k = j + 1; k < A.GetRows(); k++)
			{
				sum += std::abs(A[j][k]);
			}
		}
		if (sum == 0.0)
		{
			return;
		}
		double thresh;
		if (i < 4)
		{
			thresh = 0.2 * sum / (A.GetRows() * A.GetRows());
		} else
		{
			thresh = 0.0;
		}
		for (int j = 0; j < A.GetRows() - 1; j++)
		{
			double g;
			double h;
			double theta;
			double t;
			double c;
			double s;
			double tau;
			for (int k = j + 1; k < A.GetRows(); k++)
			{
				g = 100 * std::abs(A[j][k]);
//				if (i > 4 && g <= small * std::abs(eigenValues[j]) && g <= small * std::abs(eigenValues[k]))
				if (i > 4 && g <= small)
				{
					std::cout << g << std::endl;
					A[j][k] = 0.0;
				} else if (std::abs(A[j][k]) > thresh)
				{
					h =  eigenValues[k] - eigenValues[j];
					if (g <= small * std::abs(h))
					{
						t = (A[j][k]) / h;
					} else
					{
						theta = 0.5 * h / (A[j][k]);
						t = 1.0 / (std::abs(theta) + std::sqrt(1 + (theta * theta)));
						if (theta < 0.0)
						{
							t = -t;
						}
					}
					c = 1 / std::sqrt(1 + t * t);
					s = t * c;
					tau = s / (1 + c);
					h = t * A[j][k];
					z[j] -= h;
					z[k] += h;
					eigenValues[j] -= h;
					eigenValues[k] += h;
					A[j][k] = 0;
					for (int l = 0; l < j - 1; l++)
					{
						g = A[l][j];
						h = A[l][k];
						A[l][j] = g - s * (h + g * tau);
						A[l][k] = h + s * (g - h * tau);
					}
					for (int l = j + 1; l < k; l++)
					{
						g = A[j][l];
						h = A[l][k];
						A[j][l] = g - s * (h + g * tau);
						A[l][k] = h + s * (g - h * tau);
					}
					for (int l = k + 1; l < A.GetRows(); l++)
					{
						g = A[j][l];
						h = A[k][l];
						A[j][l] = g - s * (h + g * tau);
						A[k][l] = h + s * (g - h * tau);
					}
					for (int l = 0; l < A.GetRows(); l++)
					{
						g = eigenVectors[l][j];
						h = eigenVectors[l][k];
						eigenVectors[l][j] = g - s * (h + g * tau);
						eigenVectors[l][k] = h + s * (g - h * tau);
					}
					++nrot;
				}
			}
		}
		for (int j = 0; j < A.GetRows(); j++)
		{
			b[j] += z[j];
			eigenValues[j] = b[j];
			z[j] = 0;
		}
	}
}
*/