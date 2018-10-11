import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("./OutputData/H2l.dat")

plt.figure(1)
plt.plot(data[:, 0], data[:, 1])
plt.ylim(-1.2,1)
plt.show()
print(data)
