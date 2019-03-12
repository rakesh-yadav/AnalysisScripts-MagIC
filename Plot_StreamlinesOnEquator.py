from magic import *
from scipy.integrate import trapz
from scipy.interpolate import griddata
import pylab as P
import numpy as N
import os
        
cm = 'RdYlBu_r'
lev = 100

gr = MagicGraph()

#circular grid
phi = N.linspace(0., 2.*N.pi, gr.nphi-1)
rr, pphi = N.meshgrid(gr.radius, phi)
xx = rr * N.cos(pphi)
yy = rr * N.sin(pphi)

eq_vx = gr.vr[:,gr.ntheta/2,:]*N.cos(pphi) - gr.vphi[:,gr.ntheta/2,:]*N.sin(pphi)
eq_vy = gr.vr[:,gr.ntheta/2,:]*N.sin(pphi) + gr.vphi[:,gr.ntheta/2,:]*N.cos(pphi)

#regular grid
rec_x = N.linspace(-gr.radius[0], gr.radius[0], 1000)
rec_y = rec_x

xx2, yy2 = N.meshgrid(rec_x, rec_y)

regular_vx = griddata((xx.ravel(), yy.ravel()), eq_vx.ravel(), (xx2, yy2), method='nearest')
regular_vy = griddata((xx.ravel(), yy.ravel()), eq_vy.ravel(), (xx2, yy2), method='nearest')


fig2 = P.figure(figsize=(7,6))
ax2 = fig2.add_axes([0.05, 0.05, 0.9, 0.9])
vmax2 =  2*N.std(regular_vy)
vmin2 =  -vmax2
print vmin2
cs2 = N.linspace(vmin2, vmax2, lev)

ax2.contourf(xx, yy, gr.vphi[:,gr.ntheta/2,:], cs2, cmap='seismic', extend='both')

im = ax2.streamplot(xx2, yy2, regular_vx, regular_vy, density=[8, 8])


#ax1.axis('off')

P.show()

#fig.savefig('e6_NSD_2e9.png', dpi=100)

