from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np

fig, (ax_wav) = plt.subplots(1)
#fig = plt.figure()
steps = np.arange(0,4000,16)
x = np.arange(0,500) 
#y = 0.2*np.sin(2*np.pi*(x-t/2)/128)
frames = []

for t in steps:
    rho_wav = np.loadtxt("Wav3D_st="+str(t)+".dat",delimiter=" ")
    #rho_wav = np.loadtxt("p_max.dat",delimiter=" ")
    rho_inc = np.loadtxt("Inc3D_st="+str(t)+".dat",delimiter=" ")
    #rho_sct = rho_wav-rho_inc
    im_wav = ax_wav.imshow(rho_wav-rho_inc,cmap='hot',interpolation=None)#, vmin = -1.0, vmax = 1.0) 
    #im_inc = ax_inc.plot(x,rho_inc,color="red")#, vmin = 0, vmax = 1.0) 
    #im_wav = ax_inc.plot(x,0.2*np.sin(2*np.pi*(t/2-x)/128),color="blue")#, vmin = 0, vmax = 1.0) 
    #im_wav = ax_wav.plot(x,rho_wav,color="blue")#, vmin = 0, vmax = 1.0) 
    #ax_inc.set_xlim(-1,x.size+1)
    #ax_inc.set_ylim(-1.0,1.0)
    #ax_inc.grid()
    #im_sct = ax_sct.imshow(np.abs(rho_sct),cmap='hot',interpolation=None)#, vmin = -1.0, vmax = 1.0) 
    #plt.colorbar(im)    
    #frames.append(im_inc)
    frames.append([im_wav])

anim = animation.ArtistAnimation(fig, frames, interval=100, blit = True, repeat=False)
anim.save("scattered_field.gif",writer="pillow",fps=10)

plt.show()
