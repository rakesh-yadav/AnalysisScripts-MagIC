from magic import *
import numpy as np
import matplotlib.pyplot as P

###real to spec
gr = MagicGraph()
sh = SpectralTransforms(gr.l_max, gr.minc, gr.lm_max, gr.n_theta_max)
data = gr.Br[:,:,30] 
#print data.shape
datalm = sh.spat_spec(data)


#--------plot---------
cm = 'RdBu'
lev = 20
phi = np.linspace(-np.pi, np.pi, gr.nphi)
theta = np.linspace(np.pi/2, -np.pi/2, gr.ntheta)
pphi, ttheta = np.mgrid[-np.pi:np.pi:gr.nphi*1j,
                    np.pi/2.:-np.pi/2.:gr.ntheta*1j]
delat = 30. ; delon = 60.
circles = np.arange(delat, 90.+delat, delat).tolist()+\
                  np.arange(-delat, -90.-delat, -delat).tolist()+[0]
meridians = np.arange(-180+delon, 180, delon)
xx, yy = hammer2cart(ttheta, pphi)
vmax =  2*np.std(data)
vmin =  -vmax

cs = np.linspace(vmin, vmax, lev)

fig = P.figure(figsize=(8,4))
ax = fig.add_axes([0.01, 0.01, 0.87, 0.98])
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
ax.text(3, 1.05, 'Data', fontsize=15)
#ax.text(3.14, 1.05, '100', fontsize=12, color='white')
#ax.text(-2.76, 1.2, 'c', fontsize=20)
lmax_label = 'lmax='+'%.i' % gr.l_max
ax.text(2, -1.2, lmax_label, fontsize=15)
data = symmetrize(data, gr.minc)
im = ax.contourf(xx, yy, data, cs, cmap=P.get_cmap(cm), extend='both')
#im.set_clim(surf.min(), surf.max())

pos = ax.get_position()
l, b, w, h = pos.bounds
cax = fig.add_axes([0.9, 0.46-0.7*h/2., 0.015, 0.7*h])
mir = fig.colorbar(im, cax=cax)
#mir.set_ticks([0.5,1,1.5,2])
#cbytick_obj = P.getp(mir.ax.axes, 'yticklabels')
#P.setp(cbytick_obj, color='white')
#mir.set_ticks([4800, 4400, 4000, 3600, 3200])

ax.axis('off')
#----------------------



#=========new smaller lmax stuff
l_max_trunc = 5
m_max_trunc = l_max_trunc
lm_max_trunc = l_max_trunc*(l_max_trunc+1)/2 + l_max_trunc + 1
n_theta_max_trunc = (3*l_max_trunc + 1)/2


# Get indices location
idx_new = np.zeros((l_max_trunc+1, m_max_trunc+1), 'i')
idx_new[0:l_max_trunc+2, 0] = np.arange(l_max_trunc+1)
k = l_max_trunc+1
for m in range(1, l_max_trunc+1):
	for l in range(m, l_max_trunc+1):
		idx_new[l, m] = k
		k +=1

datalm_trunc = np.zeros((lm_max_trunc), 'Complex64')
for l in range(1, l_max_trunc+1):
	for m in range(0, l+1):
		lm = idx_new[l, m]
		datalm_trunc[lm] = datalm[sh.idx[l,m]]

#define new transform
sh_trunc = SpectralTransforms(l_max_trunc, 1, lm_max_trunc, n_theta_max_trunc)
#spec to real using truncated datalm
data_trunc =  sh_trunc.spec_spat(datalm_trunc)

#---------plot----------
cm = 'RdBu'
lev = 50
phi = np.linspace(-np.pi, np.pi, 2*n_theta_max_trunc)
theta = np.linspace(np.pi/2, -np.pi/2, n_theta_max_trunc)
pphi, ttheta = np.mgrid[-np.pi:np.pi:2*n_theta_max_trunc*1j,
                    np.pi/2.:-np.pi/2.:n_theta_max_trunc*1j]
delat = 30. ; delon = 60.
circles = np.arange(delat, 90.+delat, delat).tolist()+\
                  np.arange(-delat, -90.-delat, -delat).tolist()+[0]
meridians = np.arange(-180+delon, 180, delon)
xx, yy = hammer2cart(ttheta, pphi)
vmax =  2*np.std(data_trunc)
vmin =  -vmax

cs = np.linspace(vmin, vmax, lev)

fig = P.figure(figsize=(8,4))
ax = fig.add_axes([0.01, 0.01, 0.87, 0.98])
ax.plot(xxin, yyin, 'k-',lw=2)
ax.plot(xxout, yyout, 'k-',lw=2)
ax.text(3, 1.05, 'Data', fontsize=15)
lmax_label = 'lmax='+'%.i' % l_max_trunc
ax.text(2, -1.2, lmax_label, fontsize=15)
im = ax.contourf(xx, yy, data_trunc, cs, cmap=P.get_cmap(cm), extend='both')

pos = ax.get_position()
l, b, w, h = pos.bounds
cax = fig.add_axes([0.9, 0.46-0.7*h/2., 0.015, 0.7*h])
mir = fig.colorbar(im, cax=cax)

ax.axis('off')
#-----------------------

P.show()
