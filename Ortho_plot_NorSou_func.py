import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def Ortho_plot(data, south=True, nlongs=8, tilt=20, ind=0, skip=0, cmap='RdBu_r', \
               vmin=-30000, save_path='/home/rakesh/Desktop', dpi=150, \
               label=r'$(\nabla\times{\bf V})_z$', contour=False, \
               data_cont=[], vmin_cont=1000, l_max_trunc_cont=40, 
               ncontour=10, alpha=0.4, concolor='k'):

 data = np.clip(data, vmin, -vmin)

 lons = np.linspace(-180, 180, data.shape[0])
 lats = np.linspace(-90, 90, data.shape[1])

 if south==False:
  fig = plt.figure(figsize=(6,5))
 else:
  fig = plt.figure(figsize=(11,5))


 # North view
 if south==False:
  ax = fig.add_axes([0.001, 0.01, 0.89, 0.97])
 else:
  ax = fig.add_axes([0.005, 0.01, 0.45, 0.97])
 
 map = Basemap(projection='ortho',lat_0=90-tilt,lon_0=90,resolution=None)
 map.drawmeridians(np.arange(0,360,360/nlongs), dashes=[200,1], latmax=90)
 map.drawparallels([-60, -30, 30, 60], dashes=[200,1], latmax=90)
 map.drawparallels([0], dashes=[200,1])
 new_dat = map.transform_scalar(data.T, lons, lats, \
                                data.shape[0], \
                                data.shape[1], masked=True)
 im = map.imshow(new_dat, cmap=plt.get_cmap(cmap))
 map.drawmapboundary()
 
 if contour==True:
  # add location of shtns library for tra=uncating data
  import sys
  sys.path.append('/media/rakesh/Seagate Expansion Drive1/Data2/Polar_vortex_runs/Movies_and_plots/Paper_plots')
  try:
   from ltruncate_shtns_function import shtns_truncate
  except:
   import sys
   sys.exit('******SHTNS load unsuccessful. Try conda deactivate.*******')
  
  # truncate data
  data2 = shtns_truncate(data_cont[:-1,:], l_max_trunc=l_max_trunc_cont).T
  
  llons, llats = np.mgrid[-np.pi:np.pi:data2.shape[0]*1j,
                -np.pi/2.:np.pi/2.:data2.shape[1]*1j]
  x, y = map(llons*180./np.pi, llats*180./np.pi)
  levels = np.linspace(vmin_cont, data2.max(), num=ncontour)
  map.contour(x,y,data2,levels,colors=concolor,linewidths=[1.5], alpha=alpha)
  levels = np.linspace(data2.min(), -vmin_cont, num=ncontour)
  map.contour(x,y,data2,levels,colors=concolor,linewidths=[1.5], linestyles='dashed', alpha=alpha)  
  print('North contours plotted')
  
 # Color bar
 if south==False:
  cax = fig.add_axes([0.87, 0.04, 0.02, 0.8])
 else:
  cax = fig.add_axes([0.90, 0.04, 0.02, 0.8])
 cbar = fig.colorbar(im,cax=cax,extend='both')


 # South view
 if south!=False:
  ax2 = fig.add_axes([0.45, 0.01, 0.45, 0.97])
  map = Basemap(projection='ortho',lat_0=-(90-tilt),lon_0=90,resolution=None)
  map.drawmeridians(np.arange(0,360,360/nlongs), dashes=[200,1], latmax=90)
  map.drawparallels([-60, -30, 30, 60], dashes=[200,1], latmax=90)
  map.drawparallels([0], dashes=[200,1])
  new_dat = map.transform_scalar(data.T, lons, lats, \
                                 data.shape[0], \
                                 data.shape[1], masked=True)
  im2 = map.imshow(new_dat, cmap=plt.get_cmap(cmap))
  map.drawmapboundary()
  
  if contour==True:
   map.contour(x,y,data2,levels,colors=concolor,linewidths=[1.5], alpha=alpha)
   levels = np.linspace(data2.min(), -vmin_cont, num=ncontour)
   map.contour(x,y,data2,levels,colors=concolor,linewidths=[1.5], linestyles='dashed', alpha=alpha)  
   print('South contours plotted') 
 

 if south==False:
  fig.text(0.8, 0.9, label, fontsize=20) 
 else:
  fig.text(0.88, 0.9, label, fontsize=20)
 filename = save_path+'/img%05d.png' % (ind+skip)
 fig.savefig(filename, dpi=dpi)
 im.lines = []
 fig.clf()
 
 print('Frame saved at', filename)


