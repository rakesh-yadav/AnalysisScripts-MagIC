import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

#print test_data.shape
#juno_data = np.loadtxt('Juno_field.txt').T
#field_mag = np.absolute(juno_data)

# grid for juno field
juno_data = np.loadtxt('Juno_field.txt').T
field_mag = np.absolute(juno_data)
pphi, ttheta = np.mgrid[0:2.*np.pi:1024*1j, 0:np.pi:512*1j]
pphi = pphi.flatten(); ttheta = ttheta.flatten()


def chose_grid(grid,minB,maxB):
 if grid==250:
  data = np.loadtxt('mathematica_1.5k_to_0.25k/0.25k.dat')
 elif grid==500:
  data = np.loadtxt('mathematica_1.5k_to_0.25k/0.5k.dat')
 elif grid==750:
  data = np.loadtxt('mathematica_1.5k_to_0.25k/0.75k.dat')
 elif grid==1000:
  data = np.loadtxt('mathematica_1.5k_to_0.25k/1k.dat')
 elif grid==1250:
  data = np.loadtxt('mathematica_1.5k_to_0.25k/1.25k.dat')
 elif grid==1500:
  data = np.loadtxt('mathematica_1.5k_to_0.25k/1.5k.dat')
 x = data[:,0]
 y = data[:,1]
 z = data[:,2]
 
 #convert to lat, lon
 r = np.sqrt(x**2 + y**2 + z**2)
 ttheta1 = np.arccos(z/r)
 pphi1 = np.arctan2(y,x)+np.pi
 #pphi1, ttheta1 = np.meshgrid(phi1,theta1)
 #print phi1.min(), phi1.max(), theta1.min(), theta1.max()

 interp_data = griddata((pphi, ttheta), field_mag.flatten(), (pphi1, ttheta1), method='nearest')
 print interp_data.shape
 #plt.hist(interp_data)
 #plt.show()

 #new x,y,z with smaller radius
 x = 1.6*np.sin(ttheta1)*np.cos(pphi1)
 y = 1.6*np.sin(ttheta1)*np.sin(pphi1)
 z = 1.6*np.cos(ttheta1)

 trunc_x = x[np.where( ((minB < interp_data) & (interp_data < maxB)) )]
 trunc_y = y[np.where( ((minB < interp_data) & (interp_data < maxB)) )]
 trunc_z = z[np.where( ((minB < interp_data) & (interp_data < maxB)) )]

 return trunc_x,trunc_y,trunc_z

#------------------------
f=open("seeds.csv","w+")
f.write('x coord, y coord, z coord\n')

x1,y1,z1 = chose_grid(1500,50,65)

for i in range(len(x1)):
 to_write = "%f, %f, %f" %(x1[i], y1[i], z1[i])
 f.write(to_write+'\n')

x1,y1,z1 = chose_grid(1250,40,50)
for i in range(len(x1)):
 to_write = "%f, %f, %f" %(x1[i], y1[i], z1[i])
 f.write(to_write+'\n')

x1,y1,z1 = chose_grid(1000,30,40)
for i in range(len(x1)):
 to_write = "%f, %f, %f" %(x1[i], y1[i], z1[i])
 f.write(to_write+'\n')

x1,y1,z1 = chose_grid(750,20,30)
for i in range(len(x1)):
 to_write = "%f, %f, %f" %(x1[i], y1[i], z1[i])
 f.write(to_write+'\n')
 
x1,y1,z1 = chose_grid(500,10,20)
for i in range(len(x1)):
 to_write = "%f, %f, %f" %(x1[i], y1[i], z1[i])
 f.write(to_write+'\n')

x1,y1,z1 = chose_grid(250,0,10)
for i in range(len(x1)):
 to_write = "%f, %f, %f" %(x1[i], y1[i], z1[i])
 f.write(to_write+'\n')

f.close()
