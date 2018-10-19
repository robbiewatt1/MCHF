import matplotlib.pyplot as plt
import numpy as np

datal = np.loadtxt("./OutputData/H2l.dat")
datah = np.loadtxt("./OutputData/H2h.dat")

plt.figure(1)
plt.plot(datah[:, 0], datah[:, 1], datal[:, 0], datal[:, 1])
plt.ylim(-1.2,0.3)
plt.xlim(0.5,8)
plt.show()
print(data)
