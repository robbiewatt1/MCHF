import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("./test.dat")

plt.figure(1)
plt.plot(data[:,0],data[:,1])
plt.ylim(-1.0,5.0)
plt.show()
