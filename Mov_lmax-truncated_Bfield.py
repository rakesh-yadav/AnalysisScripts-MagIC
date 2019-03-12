from magic import *
import numpy as np
import matplotlib.pyplot as P
import pickle
from scipy.interpolate import griddata
from glob import glob

files = sorted(glob("Br_slices/Br_r*"))


counter=0
for file in files[::4]:
 print 'Working on', file
 f1 = open(file, "rb" )
 time = pickle.load(f1)
 radius = pickle.load(f1)
 data = pickle.load(f1)
 f1.close()

 #=========original grid specs=========
 l_max = (data.shape[1]*2-1)/3
 m_max = l_max
 lm_max = l_max*(l_max+1)/2 + l_max + 1
 n_theta_max = data.shape[1]#(3*l_max + 1)/2
 
 #print l_max, lm_max, n_theta_max

 #Create transform object
 sh = SpectralTransforms(l_max, 1, lm_max, n_theta_max)

 #go from real to spectral space
 datalm = sh.spat_spec(data)

 #--------plot---------
 cm = 'RdBu'
 lev = 20
 phi = np.linspace(-np.pi, np.pi, 2*n_theta_max)
 theta = np.linspace(np.pi/2, -np.pi/2, n_theta_max)
 pphi, ttheta = np.mgrid[-np.pi:np.pi:2*n_theta_max*1j,
                    np.pi/2.:-np.pi/2.:n_theta_max*1j]
 delat = 30. ; delon = 60.
 circles = np.arange(delat, 90.+delat, delat).tolist()+\
                  np.arange(-delat, -90.-delat, -delat).tolist()+[0]
 meridians = np.arange(-180+delon, 180, delon)
 xx, yy = hammer2cart(ttheta, pphi)
 vmax = 0.6#2*np.std(data)
 vmin =  -vmax

 cs = np.linspace(vmin, vmax, lev) 

 fig = P.figure(figsize=(8,10))
 ax = fig.add_axes([0.01, 0.5, 0.92, 0.48])
 for lat0 in circles:
    x0, y0 = hammer2cart(lat0*np.pi/180., phi)
    ax.plot(x0, y0, 'k--', linewidth=0.8)
 for lon0 in meridians:
    x0, y0 = hammer2cart(theta, lon0*np.pi/180.)
    ax.plot(x0, y0, 'k--', linewidth=0.8)
 xxout, yyout  = hammer2cart(theta, -np.pi)
 xxin, yyin  = hammer2cart(theta, np.pi)
 ax.plot(xxin, yyin, 'k-',lw=2)
 ax.plot(xxout, yyout, 'k-',lw=2)
 #tangent cylinder
 #North
 x0, y0 = hammer2cart(70*np.pi/180., phi)
 ax.plot(x0, y0, 'k-', linewidth=1.2)
 #South
 x0, y0 = hammer2cart(-70*np.pi/180., phi)
 ax.plot(x0, y0, 'k-', linewidth=1.2)
 #ax.text(3, 1.05, 'Data', fontsize=15)
 lmax_label = r'$\ell_{max}$='+'%.i' % l_max
 ax.text(-2.8, -1.2, lmax_label, fontsize=20)
 #data = symmetrize(data, gr.minc)
 im = ax.contourf(xx, yy, data, cs, cmap=P.get_cmap(cm), extend='both')
 #im.set_clim(surf.min(), surf.max())

 cax = fig.add_axes([0.85, 0.4, 0.015, 0.18])
 mir = fig.colorbar(im, cax=cax)

 ax.axis('off')
 #----------------------


 #=========truncated grid specs=========
 l_max_trunc = 13
 m_max_trunc = l_max_trunc
 lm_max_trunc = l_max_trunc*(l_max_trunc+1)/2 + l_max_trunc + 1
 n_theta_max_trunc = (3*l_max_trunc + 1)/2


 # Get appropriate indices location
 idx_new = np.zeros((l_max_trunc+1, m_max_trunc+1), 'i')
 idx_new[0:l_max_trunc+2, 0] = np.arange(l_max_trunc+1)
 k = l_max_trunc+1
 for m in range(1, l_max_trunc+1):
  for l in range(m, l_max_trunc+1):
   idx_new[l, m] = k
   k +=1

 # Get appropriate spectral coefficients
 datalm_trunc = np.zeros((lm_max_trunc), 'Complex64')
 for l in range(1, l_max_trunc+1):
  for m in range(0, l+1):
   lm = idx_new[l, m]
   datalm_trunc[lm] = datalm[sh.idx[l,m]]

 #define new transform object with the smaller grid
 sh_trunc = SpectralTransforms(l_max_trunc, 1, lm_max_trunc, n_theta_max_trunc)

 #Go from spectral to real
 data_trunc =  sh_trunc.spec_spat(datalm_trunc)

 #---------plot----------
 lev = 40
 phi2 = np.linspace(-np.pi, np.pi, 2*n_theta_max_trunc)
 theta2 = np.linspace(np.pi/2, -np.pi/2, n_theta_max_trunc)
 pphi2, ttheta2 = np.mgrid[-np.pi:np.pi:2*n_theta_max_trunc*1j,
                    np.pi/2.:-np.pi/2.:n_theta_max_trunc*1j]


 xx2, yy2 = hammer2cart(ttheta2, pphi2)

 cs = np.linspace(vmin, vmax, lev)

 ax2 = fig.add_axes([0.01, 0.01, 0.92, 0.48])
 for lat0 in circles:
    x0, y0 = hammer2cart(lat0*np.pi/180., phi)
    ax2.plot(x0, y0, 'k--', linewidth=0.8)
 for lon0 in meridians:
    x0, y0 = hammer2cart(theta, lon0*np.pi/180.)
    ax2.plot(x0, y0, 'k--', linewidth=0.8)
 ax2.plot(xxin, yyin, 'k-',lw=2)
 ax2.plot(xxout, yyout, 'k-',lw=2)
 #ax2.text(3, 1.05, 'Data', fontsize=15)
 #tangent cylinder
 #North
 x0, y0 = hammer2cart(70*np.pi/180., phi)
 ax2.plot(x0, y0, 'k-', linewidth=1.2)
 #South
 x0, y0 = hammer2cart(-70*np.pi/180., phi)
 ax2.plot(x0, y0, 'k-', linewidth=1.2)
 lmax_label = r'$\ell_{max}$='+'%.i' % l_max_trunc
 ax2.text(-2.8, -1.2, lmax_label, fontsize=20)

 #refine the data set for a better plot
 data_trunc = griddata((pphi2.ravel(), ttheta2.ravel()), data_trunc.ravel(), (pphi, ttheta), method='cubic')
 print data_trunc.shape

 im2 = ax2.contourf(xx, yy, data_trunc, cs, cmap=P.get_cmap(cm), extend='both')

 ax2.axis('off')
 #-----------------------
 counter = counter+1
 filename = 'Br_surf_%.4i.png' % counter
 fig.savefig(filename, dpi=100)


#P.show()


