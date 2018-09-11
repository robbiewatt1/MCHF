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

	int maxL = 2;
	int maxU = std::pow(2,1);
	BoysFunction boyFn = BoysFunction(4 * maxL, maxU, 16*maxU);
	std::ofstream outfile("./OutputData/CarbonEnergy.dat");

//	GaussianOrbital test = GaussianOrbital(0, 0, 0, 1.0 , positions[0]);
//	GaussianOrbital test2 = GaussianOrbital(0, 0, 1, 1.0 , positions[0]);

//	std::cout << test2.ElectronRepulsion(test2, test, test, boyFn) << std::endl;
//	std::cout << test.ElectronRepulsion(test, test2, test2, boyFn) << std::endl;
//	std::cout << test2.ElectronRepulsion(test2, test2, test, boyFn) << std::endl;
//	std::cout << test2.ElectronRepulsion(test, test2, test2, boyFn) << std::endl;
	Molecule mol = Molecule(1, charges, positions, maxL, boyFn, "./OrbitalData/STO6/");
	mol.CalculateEnergy();

/*
	for(int i = 1000; i >= 1; i--)
	{
		std::cout << i << std::endl;
		positions[0][0] = (double)i / 20.0;
		Molecule mol = Molecule(charges, positions, maxL, boyFn, "./OrbitalData/STO6/");
		mol.CalculateEnergy();
		double me  = mol.OscilatorStrength(0,1);

		outfile << positions[0][0] << "\t" << mol.GetEnergyLevels()[0] << "\t" << mol.GetEnergyLevels()[1] << "\t" << me << std::endl;
	}
*/
	return 0;
}
