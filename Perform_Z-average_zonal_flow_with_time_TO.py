import pickle
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy as np
from glob import glob 

ro = 1.538461538
nx = 200
ny = 200

#f1 = open( "zonal_vphi_0000.pickle", "rb" )
#time1 = pickle.load(f1)
#radius1 = pickle.load(f1)/ro
#vphi1 = pickle.load(f1)
#f1.close()

#f2 = open( "zonal_vphi_2100.pickle", "rb" )
#time2 = pickle.load(f2)
#radius2 = pickle.load(f2)/ro
#vphi2 = pickle.load(f2)
#f2.close()



#---------------------------
def grid_stuff(radius,data):
 #circular grid
 theta = np.linspace(np.pi/2, -np.pi/2, data.shape[0])
 rr, ttheta = np.meshgrid(radius, theta)
 xx = rr * np.cos(ttheta)
 yy = rr * np.sin(ttheta)
 
 #regular grid
 rec_x = np.linspace(0, 1, nx)
 rec_y = np.linspace(-1, 1, ny)
 xx2, yy2 = np.meshgrid(rec_x, rec_y)
 
 # array with number of  0 elements
 zero_mask = np.ones(xx2.shape)
 for i in range(xx2.shape[0]):
  for j in range(xx2.shape[1]):
   radius = (xx2[i,j]**2 + yy2[i,j]**2)**0.5
   #print radius
   if radius <0.35 or radius > 1.:
    zero_mask[i,j]=0.
 
 # count number of zero element
 num_zero = np.zeros(xx2.shape[1])
 for i in range(xx2.shape[1]):
  num_zero[i] = np.count_nonzero(zero_mask[:,i]==0)

 return [xx, yy, xx2, yy2, num_zero]
#----------------------------

#grid1 = grid_stuff(radius1,vphi1)
#grid2 = grid_stuff(radius2,vphi2)

#-----------------------------------------------------------
N_temp_data = [0] * nx
S_temp_data = [0] * nx
time_arr = []
files = sorted(glob("zonal_vphi_*.pickle"))[:]
for file in files:
 print 'Working on', file
 f = open( file, "rb" )
 temp_time = pickle.load(f)
 radius = pickle.load(f)/ro
 vphi = pickle.load(f)
 f.close()
 #print temp_time
 [xx,yy,xx2,yy2,num_zero] = grid_stuff(radius,vphi)
 regular_vphi = griddata((xx.ravel(), yy.ravel()), vphi.ravel(), (xx2, yy2), method='nearest')
 #north
 N_z_sum = np.sum(regular_vphi[0:xx2.shape[0]/2,:], axis=0)
 N_vphi_z_avg = N_z_sum/(ny-num_zero)/2.
 #south
 S_z_sum = np.sum(regular_vphi[xx2.shape[0]/2:,:], axis=0)
 S_vphi_z_avg = S_z_sum/(ny-num_zero)/2.
 #print 'averaging done on grid2'

 N_temp_data = np.vstack((N_temp_data,N_vphi_z_avg))
 S_temp_data = np.vstack((S_temp_data,S_vphi_z_avg))
 print temp_time
 time_arr.append(temp_time)

N_data = N_temp_data[1:,...]
S_data = S_temp_data[1:,...]

yy3 = np.linspace(0, 1, nx)
xx3 = time_arr#np.linspace(0, 1, N_data.shape[0])


cut = 0.2
fig = plt.figure(figsize=(10,5))
ax = fig.add_axes([0.15, 0.07, 0.8, 0.9])
ax.set_ylim(-0.35,1)
cs = np.linspace(-100,100, 60)

out_TC_data = (N_data[:,np.int(nx*0.35):].T + S_data[:,np.int(nx*0.35):].T)/2.0
in_TC_N_data = N_data[:,0:np.int(nx*0.35)].T
in_TC_S_data = S_data[:,0:np.int(nx*0.35)].T

im = ax.contourf(xx3, yy3[np.int(nx*0.35):],  out_TC_data, cs, cmap='seismic', extend='both')
im2 = ax.contourf(xx3, yy3[0:np.int(nx*0.35)],  in_TC_N_data, cs, cmap='seismic', extend='both')
im3 = ax.contourf(xx3, -yy3[0:np.int(nx*0.35)],  in_TC_S_data, cs, cmap='seismic', extend='both')

ax.plot(xx3,np.zeros(len(time_arr))+0.35, '--k', lw=3)
ax.plot(xx3,np.zeros(len(time_arr)), '--k', lw=3)

#print vphi_z_avg
#plt.plot(vphi_z_avg)
plt.show()


if False:
 cm = 'RdYlBu_r'
 lev = 100
 to_plot = regular_vphi
 fig = plt.figure(figsize=(4,7))
 ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
 
 vmax =  2*np.std(to_plot)
 vmin =  -vmax
 cs = np.linspace(vmin, vmax, lev)
 
 ax.contourf(xx2, yy2, to_plot, cs, cmap='seismic', extend='both'); ax.axis('off')
 #plt.imshow(to_plot, extent=(xx2.min(), xx2.max(), yy2.max(), yy2.min()))

 plt.show()

 #fig.savefig('e6_NSD_2e9.png', dpi=100)

