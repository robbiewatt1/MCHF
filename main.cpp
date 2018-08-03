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
	int maxU = 1000000;
	BoysFunction boyFn = BoysFunction(2 * maxL, maxU, 100*maxU);
	std::ofstream outfile("./OutputData/test.dat");
	for(int i = 100; i >= 1; i--)
	{
		positions[0][0] = (double)i / 6.0;
		Molecule mol = Molecule(charges, positions, maxL, boyFn);
		mol.CalculateEnergy();
		mol.GetEnergyLevels().Print();
		outfile << positions[0][0] << "\t" << mol.GetEnergyLevels()[0] << "\t" << mol.GetEnergyLevels()[1] << "\t" << std::endl;
			//	<<  mol.GetEnergyLevels()[2] << "\t" << mol.GetEnergyLevels()[3] << "\t" << mol.GetEnergyLevels()[4]  
				//<< "\t" << mol.GetEnergyLevels()[5] << "\t" << mol.GetEnergyLevels()[6] 
				//<< "\t" << mol.GetEnergyLevels()[7] << "\n";
	//	std::cout <<positions[0][0] << " " << mol.GetEnergyLevels()[0] << " " << mol.GetEnergyLevels()[1] << " " << mol.GetEnergyLevels()[2] << mol.GetEnergyLevels()[3] << " " << mol.GetEnergyLevels()[4] << " " << mol.GetEnergyLevels()[5] << std::endl;

	}



/*
    int maxL = 2;
    int maxU = 1000000;
    BoysFunction boyFn = BoysFunction(2 * maxL, maxU, 100*maxU);

    std::ofstream file("./test");
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 1000; ++j)
		{
			file << boyFn.Interpolate(i,j) << "\t";
		}
		file << "\n";
	}
*/	

	return 0;
}
