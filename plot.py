import numpy as np
import matplotlib.pyplot as plt

rho = np.loadtxt("Interphase.dat",delimiter=" ")
rho_noBall = np.loadtxt("Interphase_noBall.dat",delimiter=" ")

disk = np.loadtxt("dots_ellipse.dat",delimiter=" ")

#disk_pi2 = np.loadtxt("dots_ellipse_pi2.dat",delimiter=" ")

disk_norot = np.loadtxt("dots_ellipse_norot.dat",delimiter=" ")
#fig, ax = np.subplots(1)

plt.imshow(rho - rho_noBall, cmap='hot')
plt.plot(np.transpose(disk[:,1]),np.transpose(disk[:,2]),'.')
#plt.plot(np.transpose(disk_pi2[:,1]),np.transpose(disk_pi2[:,2]),'.')
#plt.plot(np.transpose(disk_norot[:,1]),np.transpose(disk_norot[:,2]),'.')
plt.show()
