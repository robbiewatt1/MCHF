import numpy as np
import matplotlib.pyplot as plt

# Load all the data files for atomic code
temp = 1.0
dens = 0.001
fracC1 = 1.0/3.0
fracC0 = 2.0/3.0
alpha = 1.0/137.0

atomicData = np.loadtxt("../OutputData/CarbonEnergy.dat")

radial	= atomicData[850:,0]
potLow  = atomicData[850:,1]
potHigh = atomicData[850:,2]
oscStr	= atomicData[850:,3]

potInf = -18.0
photonEnergy = potHigh - potLow 
dE_dR = np.gradient(photonEnergy) / np.gradient(radial)
Acoeff = 2.0 * photonEnergy**2 * alpha**3 * oscStr


plt.figure(1)
plt.semilogy(radial, np.exp(-(potHigh-potInf)))
plt.ylim(1e-11, 1)
plt.ylabel("boltz factor")
plt.ylabel("r")
plt.figure(2)
plt.semilogy(radial,13.6*photonEnergy)
plt.xlabel("r")
plt.ylabel("Energy")


absorption = (2.0 * dens**2 * fracC0 * fracC1 * np.pi**2 * alpha * oscStr / dE_dR) \
		   * (np.exp(-(potLow-potInf) / temp) - np.exp(-(potHigh-potInf) / temp))



plt.show()
