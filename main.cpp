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
	Vector<GaussianOrbital> orbitV(11);

	orbitV[0] = GaussianOrbital(0,0,0,0.0,0.0,0.0);
	orbitV[1] = GaussianOrbital(1,0,0,0.0,0.0,0.0);
	orbitV[2] = GaussianOrbital(0,1,0,0.0,0.0,0.0);
	orbitV[3] = GaussianOrbital(0,0,1,0.0,0.0,0.0);
	orbitV[4] = GaussianOrbital(0,1,1,0.0,0.0,0.0);
	orbitV[5] = GaussianOrbital(1,0,1,0.0,0.0,0.0);
	orbitV[6] = GaussianOrbital(1,1,0,0.0,0.0,0.0);
	orbitV[7] = GaussianOrbital(0,0,2,0.0,0.0,0.0);
	orbitV[8] = GaussianOrbital(0,2,0,0.0,0.0,0.0);
	orbitV[9] = GaussianOrbital(2,0,0,0.0,0.0,0.0);
	orbitV[10] = GaussianOrbital(3,3,3,0.0,0.0,0.0);	

	Matrix<double> SMatrix(10,10);
	Matrix<double> EMatrix(10,10);
	Matrix<double> eigenVectors(10,10);
	Vector<double> eigenvalues(10);

	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < 10; j++)
		{
			//EMatrix[i][j] = 0.333333 * orbitV[i].KineticOverlap(orbitV[j]) + orbitV[i].NuclearOverlap(orbitV[j],1,0,0,0);
			EMatrix[i][j] = orbitV[i].KineticOverlap(orbitV[j]);
			SMatrix[i][j] = orbitV[i].Overlap(orbitV[j]) * orbitV[j].GetNormaliseConstant() * orbitV[i].GetNormaliseConstant();
		}
	}
	EMatrix.Print();
	LinearAlgebra::GeneralisedEigenSolver(EMatrix,SMatrix,eigenVectors,eigenvalues);
	eigenVectors.Print();
	eigenvalues.Print();

	std::cout << orbitV[10].KineticOverlap(orbitV[10]) << std::endl;
	return 0;
}
