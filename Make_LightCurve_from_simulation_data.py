from __future__ import division
import os
from magic import *
import numpy as N
import pylab as P
import multiprocessing
from potential import *

#startTime = datetime.now()


#unit normal vector of plane defined by points a, b, and c
def unit_normal(a, b, c):
    x = N.linalg.det([[1,a[1],a[2]],
             [1,b[1],b[2]],
             [1,c[1],c[2]]])
    y = N.linalg.det([[a[0],1,a[2]],
             [b[0],1,b[2]],
             [c[0],1,c[2]]])
    z = N.linalg.det([[a[0],a[1],1],
             [b[0],b[1],1],
             [c[0],c[1],1]])
    magnitude = N.sqrt(x**2 + y**2 + z**2)
    return (x/magnitude, y/magnitude, z/magnitude)

#area of polygon poly
def poly_area(poly):
    if len(poly) < 3: # not a plane - no area
        return 0

    total = [0, 0, 0]
    for i in range(len(poly)):
        vi1 = poly[i]
        if i is len(poly)-1:
            vi2 = poly[0]
        else:
            vi2 = poly[i+1]
        prod = N.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = N.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2)

# convert from hammer projection to cartesian
def hammer2cart(ttheta, pphi):
    xx = 2*N.sqrt(2) * N.cos(ttheta)*N.sin(pphi/2)\
         /N.sqrt(1+N.cos(ttheta)*N.cos(pphi/2))
    yy = N.sqrt(2) * N.sin(ttheta)\
         /N.sqrt(1+N.cos(ttheta)*N.cos(pphi/2))
    return xx, yy

# convert from spherical to cartesin
def sph2cart(r, t, p):
    x = r*N.cos(p)*N.sin(t)
    y = r*N.sin(p)*N.sin(t)
    z = r*N.cos(t)
    return x, y, z
    
# give angle between two vectors
def angle(p1, p2):
	p1_cart = sph2cart(p1[0], p1[1], p1[2])
	p2_cart = sph2cart(p2[0], p2[1], p2[2])
	
	angle = N.arccos(N.dot(p1_cart, p2_cart)/N.linalg.norm(p1_cart)/N.linalg.norm(p2_cart))
	return angle

#calculate total intensity seen at some far-away point
def intensity(data, view_point, r, areas, dark_coef):
	theta = gr.colatitude
	phi  = N.linspace(0, 2*N.pi, gr.nphi)
	result = 0
	for i in range(gr.ntheta-1):
		for j in range(gr.nphi-2):
			tot_inten_quad = 0
			visi_inten_quad = 0
			limb_dark_inten_quad = 0
			ang_points = angle(view_point, [r, theta[i], phi[j]])
			if ang_points < N.pi/2:
				# data[] is in [phi, theta, r] format and data[:,:,0] is outer boundary
				tot_inten_quad = areas[i]*(data[j,i] + data[j+1,i] + data[j+1,i+1]+data[j,i+1])/4
				visi_inten_quad = tot_inten_quad*N.cos(ang_points)
				#applying limb darkning
				limb_dark_inten_quad = visi_inten_quad*(1-dark_coef*(1-N.cos(ang_points)))
			result += limb_dark_inten_quad
	print 'View point', view_point[1]*180/N.pi, view_point[2]*180/N.pi, 'is done.'
	return [view_point[2], result]

# save longitudinal intensity values for a given inclination
def lat_files(incli, VPs, inten):
	lat_incli = incli*180.0/N.pi
	filename = '%.2f' % lat_incli
	f = open(filename, 'w')
	f = open(filename, 'r+')
	for i in range(len(inten)):
		value = '%.5f %.5f' %(VPs[i], inten[i])
		f.write(value + "\n")
	f.close()

# save png files to make movie
def mov_files(data, incli, VPs, inten, nphi, ntheta,plot_cbar=True):
	from mpl_toolkits.basemap import Basemap
	rprof = data[..., 0]
	rprof = symmetrize(rprof, gr.minc)[:,::-1]
	lons = N.linspace(-180, 180, nphi)
	lats = N.linspace(-90, 90, ntheta)
	if not os.path.exists('movie_light_curve'):
		os.mkdir('movie_light_curve')
	fig = P.figure()
	LC = fig.add_subplot(2, 1, 1)
	LC.set_xlim([-10,370])
	LC.set_ylim([inten.min(),inten.max()])
	LC.set_xlabel('Longitudes')
	LC.set_ylabel('Visible intensity')
	for i in range(len(inten)):
		LC.plot(VPs[i]*180/N.pi,inten[i],'o',mfc='b')
		map = Basemap(projection='ortho', 
				lon_0=my2ortho_lon(VPs[i])*180/N.pi, 
				lat_0=my2ortho_lat(incli)*180/N.pi, resolution=None)
		plot2 = fig.add_subplot(2, 1, 2)
		map.drawparallels(N.linspace(90, -90, 7), dashes=[2,3], linewidth=0.5)
		map.drawmeridians(N.arange(-180, 180, 30), dashes=[2,3], linewidth=0.5)
		map.drawmapboundary()
		new_dat = map.transform_scalar(rprof.T, lons, lats,
					   nphi, ntheta, masked=True)
		im = map.imshow(new_dat, cmap=P.get_cmap('gist_heat'))
		if plot_cbar==True:
			cbar = map.colorbar(im,location='right',pad="5%")
		filename = 'movie_light_curve/img%04d.png' % i
		P.savefig(filename, dpi=100)
		plot2.lines =[]
		if plot_cbar==True:
			fig.delaxes(fig.axes[2])

# function for parallel processing
def calculate(func, args):
    result = func(*args)
    return result

# benchmarking flux function
def bench_flux(theta, phi, phase_angle_phi):
	return N.sin(theta)*N.cos(phi+phase_angle_phi)	

# latitude 0...180 to 90 ... -90
# longitude 0...360 to -180...180
def my2ortho_lat(incli):
	return -(incli - N.pi/2)
def my2ortho_lon(longi):
	return (longi - N.pi)


#-----------------------load data-----------------------
gr = MagicGraph()

r = 1
r /= (1-gr.radratio) # as we give a normalised radius
ind = N.nonzero(N.where(abs(gr.radius-r) \
				== min(abs(gr.radius-r)), 1, 0))
indPlot = ind[0][0]
rad = gr.radius[indPlot] * (1.-gr.radratio)


#----------------------------------------------------------------------
#limb darkning function = I(q)=Io(1-w(1-cosq))
#q  = angle between the normal to the stellar 
#     surface and the line of sight to the observer.
#Io = the intensity of light at the center of the stellar disk (q=0)
#w  = is the wavelength dependent limb darkening coefficient (<1, ~0.3)
dark_coef = 0.5

#----------------------defining areas for a strip ----------------------
theta = gr.colatitude
phi  = N.linspace(0, 2*N.pi, gr.nphi)
areas = N.zeros(gr.ntheta-1)
for i in range(len(areas)):
	p1 = sph2cart(r, theta[i], phi[0])
	p2 = sph2cart(r, theta[i], phi[1])
	p3 = sph2cart(r, theta[i+1], phi[1])
	p4 = sph2cart(r, theta[i+1], phi[0])
	areas[i] = poly_area([p1, p2, p3, p4])
#-----------------------------------------------------------------------

#-------------------------generating main data--------------------------
save_lat_files=True
draw_hammer=False
save_movie=True

temp0, rho0, beta0 = anelprof(gr.radius, gr.strat, gr.polind, gr.g0, gr.g1, gr.g2)
flux = -rho0*temp0*rderavg(gr.entropy, eta=gr.radratio) #heat flux
# MagIC data format -> [phi,theta,r]

#data = N.zeros_like(gr.entropy)
#for i in range(gr.ntheta-1):
	#for j in range(gr.nphi-1):
		##data[j,i,0] = N.sin(theta[i])*N.cos(phi[j])
		#data[j,i,0] = 1
		#if ((100<i<130) & (200<j<260)):
			#data[j,i,0] = 0
		#if ((270<i<320) & (540<j<640)):
			#data[j,i,0] = 0

#--------Truncate flux to higher degree----------------
"""
l_trunc = 5
flux_surf=flux[..., indPlot]
anlc = N.fft.fft(flux_surf, axis=0)/(4.*N.pi*gr.npI)
rm, tm, pm = extrapolate(anlc, 1., gr.minc, l_trunc)
rsurf = N.fft.ifft(rm, axis=0)*gr.npI
trunc_flux_surf = rsurf.real
"""
#-----------------------draw hammer projection--------------------------
if draw_hammer==True:
	rprof = trunc_flux_surf
	rprof = symmetrize(rprof, gr.minc)	
	phi_fig = N.linspace(-N.pi, N.pi, gr.nphi)
	theta_fig = N.linspace(N.pi/2, -N.pi/2, gr.ntheta) #hammer goes from 90 to -90!
	pphi_fig, ttheta_fig = N.mgrid[-N.pi:N.pi:gr.nphi*1j,
						N.pi/2:-N.pi/2:gr.ntheta*1j]
	lon2_fig = pphi_fig * 180/N.pi
	lat2_fig = ttheta_fig * 180/N.pi
	
	delat = 30 ; delon = 60
	circles = N.arange(delat, 90+delat, delat).tolist()+\
			  N.arange(-delat, -90-delat, -delat).tolist() + [0]
	meridians = N.arange(-180+delon, 180, delon)


	fig = P.figure(figsize=(9,4))
	ax = fig.add_axes([0.01, 0.01, 0.87, 0.98])
	
	x, y = hammer2cart(ttheta_fig, pphi_fig)
	#    drawing guide circles and meridians
	for lat0 in circles:
		x0, y0 = hammer2cart(lat0*N.pi/180, phi_fig)
		ax.plot(x0, y0, 'k:', linewidth=0.7)
	for lon0 in meridians:
		x0, y0 = hammer2cart(theta_fig, lon0*N.pi/180)
		ax.plot(x0, y0, 'k:', linewidth=0.7)
	#    drawing outermost boundaries of hammer projection
	xxout, yyout  = hammer2cart(theta_fig, -N.pi)
	xxin, yyin  = hammer2cart(theta_fig, N.pi)
	ax.plot(xxin, yyin, 'k-')
	ax.plot(xxout, yyout, 'k-')
	
	#xdot, ydot  = hammer2cart(ttheta, pphi)
	#ax.plot(xdot, ydot, 'k.')
	
	#    drawing data contaours
	im = ax.contourf(x, y, rprof, 60, cmap=P.get_cmap('gist_heat'))
	pos = ax.get_position()
	l, b, w, h = pos.bounds
	cax = fig.add_axes([0.9, 0.51-0.7*h/2., 0.03, 0.7*h])
	mir = fig.colorbar(im, cax=cax)
	ax.axis('off')
#-----------------------------------------------------------------------


#=======================================================================
#for inclination in N.linspace(0, N.pi, 9):

for inclination in [N.pi/3]:
	view_points = []
	for i in N.linspace(0, 2*N.pi, 60):
		view_points.append([5, inclination, i])
	
	inten_array = N.zeros(len(view_points))
	reshuffled_view_points = N.zeros(len(view_points))
	
	#----------------------parallel processing part---------------------
	PROCESSES = 10
	pool = multiprocessing.Pool(PROCESSES)
	TASKS = [(intensity, (trunc_flux_surf, view_point, r, areas, dark_coef)) for view_point in view_points]
	
	results = [pool.apply_async(calculate, t) for t in TASKS] #assign tasks
	
	for ind, res in enumerate(results):                     #compute tasks
			reshuffled_view_points[ind] = res.get()[0]
			inten_array[ind] = res.get()[1]
	
	#P.plot(reshuffled_view_points*180/N.pi, inten_array)
	#reshuffled_view_points=N.loadtxt('/home/rakesh/90')[:,0]
	#inten_array=N.loadtxt('/home/rakesh/90')[:,1]
	
	if save_lat_files==True:
		lat_files(inclination, reshuffled_view_points, inten_array)
		
	# save png files for movie
	# for movie use -> ffmpeg -i img%05d.png -c:v libx264 out.mp4
	# use this to join mp4: mencoder -oac pcm -ovc copy -idx -o joined.mp4 1.mp4 2.mp4
	# for GIF use -> convert -delay 10 -loop 0 *.png animation.gif
	if save_movie==True:
		mov_files(data, inclination, reshuffled_view_points, 
					inten_array, gr.nphi, gr.ntheta, plot_cbar=False)

	pool.terminate()         # clear the instances of processes
	#-------------------------------------------------------------------



P.show()

