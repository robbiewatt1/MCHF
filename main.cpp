#include "Molecule.hh"
#include "Vector.hh"
#include <cmath>
#include "STOnGOrbit.hh"
#include "Functions.hh"

int main()
{
	

	Vector<Vector<double>> positions(2);

	Vector<double> pos1(3), pos2(3);
	double test = std::sqrt(2.0);

	pos1[0] = 2.0;
    pos1[1] = 0.0;
    pos1[2] = 0.0;

    pos2[0] = 0.0;
    pos2[1] = 0.0;
    pos2[2] = 0.0;

	positions[0] = pos1;
	positions[1] = pos2;

	Vector<double> charges(2);
	charges[0] = 1;
	charges[1] = 1;
	
	Molecule a = Molecule(charges, positions, 0);

	a.CalculateEnergy();

	a.GetEnergyLevels().Print();
	
/*
	STOnGOrbit orbit1 = STOnGOrbit("./OrbitalData/STO3Test", 0, 0, 2,0,0,1);
	STOnGOrbit orbit2 = STOnGOrbit("./OrbitalData/STO3Test", 0, 0, 1,0,0,0);
	std::cout << orbit1.NuclearOverlap(orbit2,1,0,0,1) << std::endl;
//	std::cout << orbit.GetNormaliseConstant() << std::endl;
*/
	return 0;
}
