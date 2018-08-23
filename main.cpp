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

	int maxL = 7;
	int maxU = 100;
	BoysFunction boyFn = BoysFunction(2 * maxL, maxU, 100*maxU);
	std::ofstream outfile("./OutputData/testMe.dat");
	for(int i = 10; i >= 1; i--)
	{
		srand((unsigned int) time(0));
		positions[0][0] = (double)i / 10.0;
		Molecule mol = Molecule(charges, positions, maxL, boyFn);
		mol.CalculateEnergy();
		Vector<double> me = mol.MatrixElement(0,1);
		Vector<double> me2 = mol.MatrixElement(0,2);
		Vector<double> me3 = mol.MatrixElement(0,3);
		std::cout << mol.GetEnergyLevels()[0] << std::endl;
		std::cout << (me[0] * me[0] + me[1] * me[1] + me[2] * me[2]) + (me2[0] * me2[0] + me2[1] * me2[1] + me2[2] * me2[2]) + (me3[0] * me3[0] + me3[1] * me3[1] + me3[2] * me3[2]) << std::endl;

		outfile << positions[0][0] << "\t" << me[0] << "\t" << me[1] << "\t" << me[2] << std::endl;
		//outfile << positions[0][0] << "\t" << mol.GetEnergyLevels()[0] << "\t" << mol.GetEnergyLevels()[1] << "\t" << std::endl;
			//	<<  mol.GetEnergyLevels()[2] << "\t" << mol.GetEnergyLevels()[3] << "\t" << mol.GetEnergyLevels()[4]  
				//<< "\t" << mol.GetEnergyLevels()[5] << "\t" << mol.GetEnergyLevels()[6] 
				//<< "\t" << mol.GetEnergyLevels()[7] << "\n";
	//	std::cout <<positions[0][0] << " " << mol.GetEnergyLevels()[0] << " " << mol.GetEnergyLevels()[1] << " " << mol.GetEnergyLevels()[2] << mol.GetEnergyLevels()[3] << " " << mol.GetEnergyLevels()[4] << " " << mol.GetEnergyLevels()[5] << std::endl;

	}
	return 0;
}
