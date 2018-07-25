#include <cmath>
#include <iostream>
#include <fstream>

#include "Molecule.hh"
#include "Vector.hh"
#include "STOnGOrbit.hh"
#include "Numerics.hh"
#include "GaussianOrbital.hh"

int main(int argc, char* argv[])
{
	/*
	Vector<double> pos(3);
	GaussianOrbital test = GaussianOrbital(3,3,3,1,pos);
	std::cout << test.NuclearFunction(2,0,0,4,3,1,1,-1,3,10) << std::endl;
*/
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
		xaxis[i] = -5 + 10.0*(double)i/100;
		yaxis[i] = -5 + 10.0*(double)i/100;
		zaxis[i] = -5 + 10.0*(double)i/100;
	}
	Molecule mol = Molecule(charges, positions, 2);
	mol.CalculateEnergy();
	mol.GetEnergyLevels().Print();
	mol.GetBasisCoefficients().Print();
	mol.SetXAxis(xaxis);
	mol.SetYAxis(yaxis);
	mol.SetZAxis(zaxis);
//	mol.OutputData(7, "./Paraview/test2");
	
//	Array3D<double> wavefn1 = mol.CalculateWavefunction(0);
//	Array3D<double> wavefn2 = mol.CalculateWavefunction(4);
//	Array3D<double> wavefn3 = wavefn1 * wavefn2;
//	std::cout << Numerics::SimpsonsRule3D(xaxis,yaxis,zaxis,wavefn3) << std::endl;


//	Vector<double> loc(3);
//	loc[0] = 0;
//	loc[1] = 0;
//	loc[2] = 0;
//	STOnGOrbit test1 = STOnGOrbit("./OrbitalData/test", 0, 0, 0, loc);
//	STOnGOrbit test2 = STOnGOrbit("./OrbitalData/test", 0, 0, 2, loc);
//	std::cout << test1.Overlap(test2) << std::endl;
//	std::cout << test2.GetNormaliseConstant() << std::endl;


	return 0;
}
