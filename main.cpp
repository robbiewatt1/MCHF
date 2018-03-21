#include "Molecule.hh"
#include "Vector.hh"
#include <cmath>
#include "STOnGOrbit.hh"
#include "Functions.hh"
#include "Matrix.hh"
#include "LinearAlgebra.hh"

int main()
{
/*	

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
*/	

	Vector<double> pos1(3), pos2(3);
	pos1[0] = 0;
        pos1[1] = 0;
        pos1[2] = 1;

	pos2[0] = 0;
	pos2[1] = 0;
	pos2[2] = -1;


	STOnGOrbit orbit1 = STOnGOrbit("./OrbitalData/STO3Test", 0, 0, 0, pos1);
	STOnGOrbit orbit2 = STOnGOrbit("./OrbitalData/STO3Test", 0, 0, 0, pos2);
	
	Matrix<double>  energy(2,2);
	Matrix<double>	overlap(2,2);
	Matrix<double> eigV(2,2);
	Vector<double> eigE(2);
	energy[0][0] = orbit1.NuclearOverlap(orbit1,1,pos1) + 0 + orbit1.KineticOverlap(orbit1);
	energy[0][1] = orbit1.NuclearOverlap(orbit2,1,pos1) + 0 + orbit1.KineticOverlap(orbit2);
	energy[1][1] = orbit2.NuclearOverlap(orbit2,1,pos2) + 0 + orbit2.KineticOverlap(orbit2);
	energy[1][0] = energy[0][1];
	
	overlap[0][0] =1;
	overlap[0][1] = orbit1.Overlap(orbit2);
	overlap[1][0] = overlap[0][1];
	overlap[1][1] = 1;
	energy.Print();
	
	LinearAlgebra::GeneralisedEigenSolver(energy, overlap, eigV, eigE);
	
	eigE.Print();
	eigV.Print();	

	return 0;

}
