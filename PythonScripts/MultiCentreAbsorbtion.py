import numpy as np
import matplotlib.pyplot as plt

# Load all the data files for atomic code
atomicData = np.loadtxt("../OutputData/Carbon.dat")

radial	= atomicData[:,0]
potLow  = atomicData[:,1]
potHigh = atomicData[:,2]
oscStr	= atomicData[:,3]

potInf =  potHigh[0]
photonEnergy = potHigh - potLow 

print((potHigh - potInf))

plt.figure(1)
plt.plot(radial, np.exp(-(potHigh - potInf)))
plt.show()