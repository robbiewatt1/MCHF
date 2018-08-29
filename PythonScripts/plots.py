import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("./OutputData/H2+OscStr.dat")

plt.figure(1)
plt.semilogy(data[:, 0], data[:, 1])
plt.show()
print(data)
