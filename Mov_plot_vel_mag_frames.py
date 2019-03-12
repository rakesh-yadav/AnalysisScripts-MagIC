from __future__ import division
import numpy as N
import pylab as P
from mpl_toolkits.basemap import Basemap
from magic import *
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



mov_frames_mag =  movie.Movie(file="Br_R=C6_mov.tag_zzzzl", iplot=False)
mov_frames_vel =  movie.Movie(file="Vr_R=C1_mov.tag_zzzzl", iplot=False)



for i in range(0,100):
	data= mov_frames_vel.data[i,...]
	data= N.clip(data,-1500, 1500)
	cmap = 'bwr'
	
	rprof = data
	lons = N.linspace(-180, 180, data.shape[0])
	lats = N.linspace(-90, 90, data.shape[1])
	
	fig = P.figure(figsize=(10,10))
	
	#-------------------------------------------------------------------
	ax1 = fig.add_axes([0.01, 0.01, 0.47, 0.47])
	map = Basemap(projection='ortho',lat_0=-50,lon_0=100,resolution=None)
	new_dat = map.transform_scalar(rprof.T, lons, lats,
				   data.shape[0], data.shape[1], masked=True)
	im1 = map.imshow(new_dat, cmap=cmap)
	ax1.text(0.15, 0.1, 'South', color='w', transform=ax1.transAxes, ha='right')
	#para = map.drawparallels([-84.29,0], color='k', dashes=[200,1], linewidth=2)
	#del ax1.lines[2]
	map.drawmeridians(N.arange(0,360,60), color='k', dashes=[20,10], latmax=90)
	map.drawparallels([-60, -30, 30, 60], color='k', dashes=[20,10], latmax=90)
	map.drawmapboundary()
	#cbar = map.colorbar(im,location='top',pad="7%", extend='both')
	#-------------------------------------------------------------------
	ax2 = fig.add_axes([0.01, 0.5, 0.47, 0.47])
	map = Basemap(projection='ortho',lat_0=50,lon_0=100,resolution=None)
	new_dat = map.transform_scalar(rprof.T, lons, lats,
				   data.shape[0], data.shape[1], masked=True)
	im2 = map.imshow(new_dat, cmap=cmap)
	ax2.text(0.69, 1.02, 'Radial velocity', color='w', transform=ax2.transAxes, ha='right', fontsize=18)
	ax2.text(0.15, 0.9, 'North', color='w', transform=ax2.transAxes, ha='right')
	ax2.text(0.15, 0.04, r'$\pm 1500$', color='w', transform=ax2.transAxes, ha='right', fontsize=16)
	ax2.text(0.17, 0., '(Reyn. Num.)', color='w', transform=ax2.transAxes, ha='right', fontsize=10)
	map.drawmeridians(N.arange(0,360,60), color='k', dashes=[20,10], latmax=90)
	map.drawparallels([-60, -30, 30, 60], color='k', dashes=[20,10], latmax=90)
	map.drawmapboundary()
	#-------------------------------------------------------------------
	#----------------------------magnetic-------------------------------
	data= mov_frames_mag.data[i,...]*2300/5.97
	data= N.clip(data,-10000, 10000)
	cmap = plt.get_cmap('hot_ice') 
	
	rprof = data#brsurf
	
	ax3 = fig.add_axes([0.5, 0.01, 0.47, 0.47])
	map = Basemap(projection='ortho',lat_0=-50,lon_0=100,resolution=None)
	new_dat = map.transform_scalar(rprof.T, lons, lats,
				   data.shape[0], data.shape[1], masked=True)
	im3 = map.imshow(new_dat, cmap=cmap)
	
	map.drawmeridians(N.arange(0,360,60), color='w', dashes=[20,10], latmax=90)
	map.drawparallels([-60, -30, 30, 60], color='w', dashes=[20,10], latmax=90)
	map.drawmapboundary()
	#-------------------------------------------------------------------
	ax4 = fig.add_axes([0.5, 0.5, 0.47, 0.47])
	map = Basemap(projection='ortho',lat_0=50,lon_0=100,resolution=None)
	new_dat = map.transform_scalar(rprof.T, lons, lats,
				   data.shape[0], data.shape[1], masked=True)
	im4 = map.imshow(new_dat, cmap=cmap)
	ax4.text(0.77, 1.02, 'Radial magnetic field', color='w', transform=ax4.transAxes, ha='right', fontsize=18)
	ax4.text(1, 0, r'$\pm 10\,kG$', color='w', transform=ax4.transAxes, ha='right', fontsize=16)
	map.drawmeridians(N.arange(0,360,60), color='w', dashes=[20,10], latmax=90)
	map.drawparallels([-60, -30, 30, 60], color='w', dashes=[20,10], latmax=90)
	map.drawmapboundary()
	#-------------------------------------------------------------------
	
	
	filename = 'movie/img%05d.png' % i
	fig.savefig(filename, dpi=100, facecolor='0.1')
	im1.lines = []
	im2.lines = []
	im3.lines = []
	im4.lines = []
	print 'Frame',i+1,'of',mov_frames_vel.data.shape[0],'is Done.'

