from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np

fig = plt.figure(figsize=(10,10))

steps = np.arange(0,512,2)
x = np.arange(0,63)

frames = []

for t in steps:
    rho = np.loadtxt("Waves3D_st="+str(t)+".dat",delimiter=" ")
    im = plt.imshow(np.abs(rho),cmap='hot',interpolation=None, vmin = 0, vmax = 0.4) 
    #plt.colorbar(im)    
    frames.append([im])

anim = animation.ArtistAnimation(fig, frames, interval=50, blit = True, repeat_delay=0)
plt.colorbar();
anim.save("IBM_2Dball.gif",writer="imagemagick",fps=10)
plt.show()
