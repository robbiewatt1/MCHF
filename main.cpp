#include "Functions.hh"
#include "Orbital.hh"
#include "GaussianOrbital.hh"
#include "Constants.hh"
#include "Array3D.hh"
#include "Numerics.hh"
#include "HydrogenicOrbital.hh"
#include "LinearAlgebra.hh"
#include "Matrix.hh"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <climits>
#include <string>
#include "H5Cpp.h"
#include <cmath>
int main()
{
	Vector<GaussianOrbital> orbitV(4);

	orbitV[0] = GaussianOrbital(0,0,0,0.0,0.0,0.0);
//	orbitV[1] = GaussianOrbital(1,0,0,0.0,0.0,0.0);
//	orbitV[2] = GaussianOrbital(0,1,0,0.0,0.0,0.0);
//	orbitV[3] = GaussianOrbital(0,0,1,0.0,0.0,0.0);
//	orbitV[4] = GaussianOrbital(0,1,1,0.0,0.0,0.0);
//	orbitV[5] = GaussianOrbital(1,0,1,0.0,0.0,0.0);
//	orbitV[6] = GaussianOrbital(1,1,0,0.0,0.0,0.0);
	orbitV[1] = GaussianOrbital(0,0,2,0.0,0.0,0.0);
	orbitV[2] = GaussianOrbital(0,2,0,0.0,0.0,0.0);
	orbitV[3] = GaussianOrbital(2,0,0,0.0,0.0,0.0);

	Matrix<double> SMatrix(4,4);
	Matrix<double> EMatrix(4,4);
	Matrix<double> eigenVectors(4,4);
	Vector<double> eigenvalues(4);

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			EMatrix[i][j] = orbitV[i].KineticOverlap(orbitV[j]) + orbitV[i].NuclearOverlap(orbitV[j],1,0,0,0);
			SMatrix[i][j] = orbitV[i].Overlap(orbitV[j]) * orbitV[j].GetNormaliseConstant() * orbitV[i].GetNormaliseConstant();
		}
	}


//	LinearAlgebra::EigenSolver(EMatrix,eigenVectors,eigenvalues);
	LinearAlgebra::GeneralisedEigenSolver(EMatrix,SMatrix,eigenVectors,eigenvalues);
//	SMatrix.Print();
//	EMatrix.Print();
//	eigenVectors.Print();
//	eigenvalues.Print();
	return 0;
}
