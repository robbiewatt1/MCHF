#include <cmath>
#include "BoysFunction.hh"
#include "Functions.hh"

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
		m_uAxis[i] =  maxU * (double)i / (m_uLength - 1.0); 
	}
	for (int i = 0; i < m_vLength; i++)
	{
		m_vAxis[i] = i;
	}
	// Generater the boysfunction data
	BoysGenerator();
	m_boysData.Print();
}

BoysFunction::~BoysFunction()
{
}

double BoysFunction::Interpolate(int v, double u) const
{
	// Find the u vlaues that are next to the value u
	int loc1 = u /  m_uAxis.End() * (m_uLength - 1.0);
	int loc0, loc2, loc3;
	if(loc1 == 0)
	{
		loc0 = 0;
		loc2 = loc1 + 1;
		loc3 = loc1 + 2;
		std::cout << "here1" << std::endl;
	} else if (loc1 == m_uLength - 1)
	{		
		loc0 = loc1 - 1;
		loc2 = loc1 + 1;
		loc3 = loc1 + 1;
		std::cout << "here2" << std::endl;
	} else if (loc1 == m_uLength)
	{
		loc0 = loc1 - 1;
		loc2 = loc1 + 1;
		loc3 = loc1 + 2;
		std::cout << "here3" << std::endl;
	} else
	{
		loc0 = loc1 - 1;
		loc2 = loc1 + 1;
		loc3 = loc1 + 2;
	}
	return m_boysData[v][loc1] + 0.5 * u * (m_boysData[v][loc2] - m_boysData[v][loc0] + u 
		   * (2.0 * m_boysData[v][loc0] - 5.0 * m_boysData[v][loc1] + 4.0 * m_boysData[v][loc2] 
		   - m_boysData[v][loc3] + u * (3.0 * (m_boysData[v][loc1] - m_boysData[v][loc2]) 
		   + m_boysData[v][loc3] - m_boysData[v][loc0])));
}

void BoysFunction::BoysGenerator()
{
	m_boysData = Matrix<double>(m_vAxis.Length(), m_uAxis.Length());
	for(int i = 0; i < m_uAxis.Length(); i++)
	{
		double result_p;
		if ( m_uAxis[i] > 1.0)
		{
			double sum = 0;
			for(int j = 0; j < 20; j++)
			{
				sum += Functions::SemiFactorial(2 * m_vAxis.End() - 1) * std::pow(2 * m_uAxis[i], j) 
					  / Functions::SemiFactorial(2 * m_vAxis.End() + 2 * j + 1);
			}
			result_p = std::exp( -1.0 * m_uAxis[i]) * sum;
		} else
		{
			result_p = (1.0 / (2.0 * m_vAxis.End() + 1.0))
					 - (m_uAxis[i] / (2.0 * m_vAxis.End() + 3.0))
					 + (m_uAxis[i] * m_uAxis[i] / (2.0 * (2.0 * m_vAxis.End() + 5.0)))
					 - (m_uAxis[i] * m_uAxis[i] * m_uAxis[i] / (6.0 * (2.0 * m_vAxis.End() + 7.0)))
					 + (m_uAxis[i] * m_uAxis[i] * m_uAxis[i] * m_uAxis[i] / (24.0 * (2.0 * m_vAxis.End() + 11.0)))
					 - (m_uAxis[i] * m_uAxis[i] * m_uAxis[i] * m_uAxis[i] * m_uAxis[i] / (120.0 * (2.0 * m_vAxis.End() + 13.0)))
					 + (m_uAxis[i] * m_uAxis[i] * m_uAxis[i] * m_uAxis[i] * m_uAxis[i] * m_uAxis[i] / (720.0 * (2.0 * m_vAxis.End() + 15.0)));
		}
		m_boysData[m_vLength-1][i] = result_p;
		for(int k = m_vLength - 2; k > -1; k--)
		{
			m_boysData[k][i] = (2 * m_uAxis[i] * m_boysData[k+1][i] + std::exp(-m_uAxis[i]))
						 / (2 * (k) + 1);
		}
	}
}