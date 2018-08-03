#include <cmath>
#include <omp.h>
#include "BoysFunction.hh"
#include "Functions.hh"
#include "Constants.hh"

BoysFunction::BoysFunction()
{
}

BoysFunction::BoysFunction(int maxV, double maxU, int uLength)
{
	// set the v and u axis
	m_uLength = uLength;
	m_vLength = maxV + 1;
	m_uAxis = Vector<double>(uLength);
	m_vAxis = Vector<int>(maxV+1);
	m_boysData = Matrix<double>(maxV + 1, uLength);

	for (int i = 0; i < m_uLength; i++)
	{
		m_uAxis[i] =  (maxU * (double)i / (m_uLength - 1.0));// * (maxU * (double)i / (m_uLength - 1.0)) / maxU;
	}
	for (int i = 0; i < m_vLength; i++)
	{
		m_vAxis[i] = i;
	}
	// Generater the boysfunction data
	BoysGenerator();
}

BoysFunction::~BoysFunction()
{
}

double BoysFunction::Interpolate(int v, double u) const
{
	int loc1 = u /  m_uAxis.End() * (m_uLength - 1.0);
	int iIndex, jIndex;
	double boysValue(0);
	if(loc1 == 0)
	{
		for (int i = 0; i < 3; i++)
		{
			iIndex = loc1 + i;
			double poln(1.0);
			for (int j = 0; j < 3; j++)
			{
				if(j != i)
				{
					jIndex = loc1 + j;
					poln *=  (u - m_uAxis[jIndex]) / (m_uAxis[iIndex] - m_uAxis[jIndex]);
				}
			}
			boysValue += poln * m_boysData[v][iIndex];
		}
	} else if (loc1 == m_uLength - 2)
	{
		for (int i = 0; i < 3; i++)
		{
			iIndex = loc1 + i - 1;
			double poln(1.0);
			for (int j = 0; j < 3; j++)
			{
				if(j != i)
				{
					jIndex = loc1 + j - 1;
					poln *=  (u - m_uAxis[jIndex]) / (m_uAxis[iIndex] - m_uAxis[jIndex]);
				}
			}
			boysValue += poln * m_boysData[v][iIndex];
		}
	} else if (loc1 >= m_uLength - 1)
	{
		std::cerr << "Error: Boys function range too low!" << std::endl;
		std::cerr << " ABORTING" << std::endl;
		std::exit(0);
	} else
	{
		for (int i = 0; i < 4; i++)
		{
			iIndex = loc1 + i - 1;
			if(iIndex < 0)
			{
				std::cerr << "error" << std::endl;
			}
			double poln(1.0);
			for (int j = 0; j < 4; j++)
			{
				if(j != i)
				{
					jIndex = loc1 + j - 1;
					poln *=  (u - m_uAxis[jIndex]) / (m_uAxis[iIndex] - m_uAxis[jIndex]);
				}
			}
			boysValue += poln * m_boysData[v][iIndex];
		}
	}
	return boysValue;
}

void BoysFunction::BoysGenerator()
{
	m_boysData = Matrix<double>(m_vAxis.Length(), m_uAxis.Length());
	#pragma omp parallel for
	for(int i = 0; i < m_uAxis.Length(); i++)
	{
		double result_p;
		if (m_uAxis[i] < 1.0)
		{
			result_p = (1.0 / (2.0 * m_vAxis.End() + 1.0))
					 - (m_uAxis[i] / (2.0 * m_vAxis.End() + 3.0))
					 + (m_uAxis[i] * m_uAxis[i] / (2.0 * (2.0 * m_vAxis.End() + 5.0)))
					 - (m_uAxis[i] * m_uAxis[i] * m_uAxis[i] / (6.0 * (2.0 * m_vAxis.End() + 7.0)))
					 + (m_uAxis[i] * m_uAxis[i] * m_uAxis[i] * m_uAxis[i] / (24.0 * (2.0 * m_vAxis.End() + 11.0)))
					 - (m_uAxis[i] * m_uAxis[i] * m_uAxis[i] * m_uAxis[i] * m_uAxis[i] / (120.0 * (2.0 * m_vAxis.End() + 13.0)));
					 + (m_uAxis[i] * m_uAxis[i] * m_uAxis[i] * m_uAxis[i] * m_uAxis[i] * m_uAxis[i] / (720.0 * (2.0 * m_vAxis.End() + 15.0)));
		} else if(m_uAxis[i] < 20.0)
		{
			double sum = 0;
			for(int j = 0; j < 30; j++)
			{
				sum += Functions::SemiFactorial(2.0 * m_vAxis.End() - 1.0) * std::pow(2.0 * m_uAxis[i], j) 
					  / Functions::SemiFactorial(2.0 * m_vAxis.End() + 2.0 * j + 1.0);
			}
			result_p = std::exp( -1.0 * m_uAxis[i]) * sum;
		} else
		{
			result_p = Functions::SemiFactorial(2 * m_vAxis.End() - 1) * std::sqrt(Constants::pi / std::pow(m_uAxis[i], 2 * m_vAxis.End() + 1)) 
					   / std::pow(2.0, m_vAxis.End() + 1);
		} 
		m_boysData[m_vLength-1][i] = result_p;
		for(int k = m_vLength - 2; k > -1; k--)
		{
			m_boysData[k][i] = (2.0 * m_uAxis[i] * m_boysData[k+1][i] + std::exp(-m_uAxis[i]))
						 / (2.0 * k + 1.0);
		}
	}
}