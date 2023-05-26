from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np

fig, ((ax_wav,ax_inc,ax_tav),(ax_wx,ax_ix,ax_tx)) = plt.subplots(2,3)
#fig = plt.figure()
steps = np.arange(0,800,2)
x = np.arange(0,250) 
#y = 0.2*np.sin(2*np.pi*(x-t/2)/128)
frames = []
rho_tav = np.zeros((64,250))
for t in steps:
    rho_wav = np.loadtxt("Wav3D_st="+str(t)+".dat",delimiter=" ")
    #rho_wav = np.loadtxt("p_max.dat",delimiter=" ")
    rho_inc = np.loadtxt("Inc3D_st="+str(t)+".dat",delimiter=" ")
    rho_sct = rho_wav-rho_inc
    rho_tav = rho_tav + rho_sct 
    im_wav = ax_wav.imshow(rho_wav,cmap='hot')#,interpolation=None,vmin=-1.,vmax=1.)
    ax_wav.set_title("Total field")
    im_inc = ax_inc.imshow(rho_sct,cmap='hot')#,interpolation=None,vmin=-4e-2,vmax=4e-2)
    ax_inc.set_title("Scattered field")
    im_tav = ax_tav.imshow(rho_tav / 1000,cmap='hot')#,interpolation=None,vmin=-6e-3,vmax=6e-3)
    ax_tav.set_title("Time average")
    im_wx = ax_wx.plot(x,rho_wav[:][31],'bs')#,color="red")#, vmin = 0, vmax = 1.0) 
    ax_wx.set_title("Along x-axis")
    im_ix = ax_ix.plot(x,rho_sct[:][31],'bs')#,color="red")#, vmin = 0, vmax = 1.0) 
    ax_ix.set_title("Along x-axis")
    im_tx = ax_tx.plot(x,rho_tav[:][31]/1000,'bs')#,color="green")#, vmin = 0, vmax = 1.0) 
    ax_tx.set_title("Along x-axis")
    #im_wav = ax_wav.plot(x,rho_wav,color="blue")#, vmin = 0, vmax = 1.0) 
    #ax_inc.set_xlim(-1,x.size+1)
    #ax_inc.set_ylim(-1.0,1.0)
    #ax_inc.grid()
    #im_sct = ax_sct.imshow(np.abs(rho_sct),cmap='hot',interpolation=None)#, vmin = -1.0, vmax = 1.0) 
    #plt.colorbar(im)    
    frames.append(im_wx)
    frames.append(im_ix)
    frames.append(im_tx)
    frames.append([im_wav])
    frames.append([im_inc])
    frames.append([im_tav])

anim = animation.ArtistAnimation(fig, frames, interval=10, blit = True, repeat=True)
#anim.save("all_fields_in_quarter.gif",writer="pillow",fps=10)

plt.show()
