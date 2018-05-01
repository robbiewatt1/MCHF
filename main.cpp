#include <cmath>
#include <iostream>
#include <fstream>

#include "Molecule.hh"
#include "Vector.hh"
#include "STOnGOrbit.hh"
#include "Numerics.hh"

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
/*
	std::ofstream outfile("./OutputData/HeH+2.dat");
	Molecule mol = Molecule(charges, positions, 6);
	mol.CalculateEnergy();
	outfile << positions[0][0] << "\t" << mol.GetEnergyLevels()[0] << "\t" << mol.GetEnergyLevels()[1] << "\t" 
			<<  mol.GetEnergyLevels()[2] << "\t" << mol.GetEnergyLevels()[3] << "\t" << mol.GetEnergyLevels()[4]  
			<< "\t" << mol.GetEnergyLevels()[5] << "\t" << mol.GetEnergyLevels()[6] << "\t" << mol.GetEnergyLevels()[7] << "\n";

*/

	Vector<double> xaxis(101), yaxis(101), zaxis(101);
	for(int i = 0; i < 101; i++)
	{
		xaxis[i] = -10 + 20.0*(double)i/100;
		yaxis[i] = -10 + 20.0*(double)i/100;
		zaxis[i] = -10 + 20.0*(double)i/100;
	}
	Molecule mol = Molecule(charges, positions, 2);
	mol.CalculateEnergy();
	
	mol.SetXAxis(xaxis);
	mol.SetYAxis(yaxis);
	mol.SetZAxis(zaxis); 

	mol.CaculateMatrixElement(0,0);
	return 0;
}
