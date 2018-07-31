import matplotlib.pyplot as plt
import numpy as np


data = np.loadtxt("./BoysData")
print(data[:,9])

plt.figure(1)
plt.plot(data)
plt.show()
