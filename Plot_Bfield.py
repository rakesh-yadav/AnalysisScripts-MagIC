from __future__ import division
import numpy as N
import pylab as P
import pickle
from mpl_toolkits.basemap import Basemap
#from magic import *
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt

#------------Custom colormap------------------
cdict1 = {'red':   ((0.0, 0.0, 0.75),
                   (1/8, 0.5, 0.5),
                   (2/8, 0.0, 0.0),
                   (3/8, 0.0, 0.0),
                   (4/8, 0.0, 0.0),
                   #
                   (5/8, 0.5, 0.5),
                   (6/8, 1.0, 1.0),
                   (7/8, 1.0, 1.0),
                   (8/8, 1.0, 1.0)),

         'green': ((0.0, 0.0, 1.0),
                   (1/8, 1.0, 1.0),
                   (2/8, 0.5, 0.5),
                   (3/8, 0.0, 0.0),
                   (4/8, 0.0, 0.0),
                   #
                   (5/8, 0.0, 0.0),
                   (6/8, 0.5, 0.5),
                   (7/8, 1.0, 1.0),
                   (8/8, 1.0, 1.0)),

         'blue':  ((0.0, 0.0, 1.0),
                   (1/8, 1.0, 1.0),
                   (2/8, 1.0, 1.0),
                   (3/8, 0.5, 0.5),
                   (4/8, 0.0, 0.0),
                   #
                   (5/8, 0.0, 0.0),
                   (6/8, 0.0, 0.0),
                   (7/8, 0.5, 0.5),
                   (8/8, 0.75, 1.0)),
        }
hot_ice = LinearSegmentedColormap('hot_ice', cdict1)
plt.register_cmap(cmap=hot_ice)


fPickle = open('FC_B_Trunc_surf_l10.pickle', 'rb')
theta = pickle.load(fPickle)
phi = pickle.load(fPickle)
Br = pickle.load(fPickle)
Bt = pickle.load(fPickle)
Bp = pickle.load(fPickle)

fPickle.close()


fig = P.figure(figsize=(5,5))

#----------------------------magnetic-------------------------------
data= Br[:,::-1]*2300/5.97
data= N.clip(data,-3000, 3000)
cs = N.linspace(-3000, 3000, 30)
cmap = plt.get_cmap('hot_ice') 

lons = N.linspace(-180, 180, data.shape[0])
lats = N.linspace(-90, 90, data.shape[1])
ax1 = fig.add_axes([0.05, 0.01, 0.9, 0.9])

lons, lats = N.mgrid[-N.pi:N.pi:data.shape[0]*1j,
                    -N.pi/2.:N.pi/2.:data.shape[1]*1j]


map = Basemap(projection='ortho',lat_0=45,lon_0=0,resolution=None)
map.drawmeridians(N.arange(0,360,60), color='w', dashes=[20,10], latmax=90)
map.drawparallels([-60, -30, 30, 60], color='w', dashes=[20,10], latmax=90)
map.drawparallels([0], dashes=[200,1])
x, y = map(lons*180./N.pi, lats*180./N.pi)
im1 = map.contourf(x, y, data, cs, cmap=P.get_cmap(cmap))

im1.set_clim(-3000, 3000)
map.drawmapboundary()

cbar = map.colorbar(im1,location='top',pad="7%", extend='both')
cbar.set_ticks([-3000,-1500,0,1500,3000])
cbar.set_ticklabels(['-3 kG', '-1.5 kG', '0', '1.5 kG', '3 kG'])

fig.savefig('lowres.pdf')

P.show()
