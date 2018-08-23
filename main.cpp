#include <cmath>
#include <iostream>
#include <fstream>

#include "Molecule.hh"
#include "Vector.hh"
#include "STOnGOrbit.hh"
#include "Numerics.hh"
#include "GaussianOrbital.hh"
#include "Functions.hh"
#include "BoysFunction.hh"

int main(int argc, char* argv[])
{
	if(argc != 2)
	{
		std::cerr << "ERROR: Incorrect use of " << argv[0] << std::endl;
		std::cerr << "Correct use: " << argv[0] << " filename" << std::endl;
		return 1;
	}
	
	std::ifstream inputFile(argv[1]);
	if(!inputFile)
	{
		std::cerr << "Could not find file " << argv[1] << std::endl;
		return 1;
	}


	Vector<Vector<double>> positions;
	Vector<double> charges;
	double Z, posx, posy, posz;
	while(inputFile >> Z >> posx >> posy >> posz)
	{
		charges.Append(Z);

		Vector<double> position(3);
		position[0] = posx;
		position[1] = posy;
		position[2] = posz;
		positions.Append(position);
	}

	int maxL = 9;
	int maxU = 1;
	BoysFunction boyFn = BoysFunction(2 * maxL, maxU, 10*maxU);
	std::ofstream outfile("./OutputData/H2+OscStr.dat");
	for(int i = 60; i >= 1; i--)
	{
		srand((unsigned int) time(0));
		positions[0][0] = (double)i / 10.0;
		Molecule mol = Molecule(charges, positions, maxL, boyFn);
		mol.CalculateEnergy();
		double me  = mol.OscilatorStrength(0,1);
		double me2 = mol.OscilatorStrength(0,2);
		double me3 = mol.OscilatorStrength(0,3);
		mol.GetEnergyLevels().Print();

		outfile << positions[0][0] << "\t" << me << std::endl;
	}
	return 0;
}
