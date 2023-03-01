import numpy as np
import matplotlib.pyplot as plt

rho = np.loadtxt("Waves2D.dat",delimiter=" ")

disk = np.loadtxt("disk.dat",delimiter=" ")

#fig, ax = np.subplots(1)

plt.imshow(rho, cmap='hot')
plt.plot(np.transpose(disk[:,0]),np.transpose(disk[:,1]))
plt.show()
