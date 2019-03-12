from magic import *
import os
import matplotlib.pyplot as P
import numpy as N
from scipy.integrate import trapz

os.chdir(os.getcwd())

# which force to calcualte
Lorentz=True
Coriolis=True
Pressure=True
Inertial=True
Buoyancy=True
Viscous=True

# force cancellations
force_can_cal=True

# calculate averages inside/outside TC
inout_TC=True

# print resutls
report=True
#---------------------------------------------
gr = MagicGraph(ivar=1)

#-----------in/out TC mask--------
rr2D = N.zeros((gr.ntheta, gr.nr), dtype='Float32')
th2D = N.zeros_like(rr2D)
for i in range(gr.ntheta):
  th2D[i, :] = gr.colatitude[i]
  rr2D[i, :] = gr.radius
ss2D = rr2D*N.sin(th2D) #Cylindrical radius (theta is from 0 to Pi)
#---------
in_TC_map = N.zeros_like(rr2D)
for i in range(ss2D.shape[0]):
 for j in range(ss2D.shape[1]):
  if ss2D[i,j] < 0.34999999:
   in_TC_map[i,j]=1.
out_TC_map=-(-1+in_TC_map)

temp = N.trapz(in_TC_map, gr.colatitude, axis=0)
area_in =  N.trapz(temp*gr.radius, gr.radius)
temp = N.trapz(out_TC_map, gr.colatitude, axis=0)
area_out = N.trapz(temp*gr.radius, gr.radius)
#--------------------------------


th3D = N.zeros_like(gr.Bphi)
rr3D = N.zeros_like(th3D)
for i in range(gr.nr):
    rr3D[:, :, i] = gr.radius[i]
for i in range(gr.ntheta):
    th3D[:, i, :] = gr.colatitude[i]

#----------get spherically averaged quantities------
def get_avg(data):
  #Checked using total volume of the shell defined by "data". 
  #Produces "5.72537" (with negative sign, 
  #radius array is inverted), analytical value "5.74017" for 0.35 shells.

  Phi_avg = N.trapz(data, dx=2*N.pi/gr.nphi, axis=0)/2./N.pi
  PhiTheta_avg = N.trapz(Phi_avg*N.sin(th2D), th2D, axis=0)/2.
  PhiThetaR_avg = N.trapz(PhiTheta_avg*gr.radius**2, gr.radius)
  #normalize the radial integration
  PhiThetaR_avg = PhiThetaR_avg/N.trapz((1.+N.zeros_like(PhiTheta_avg))*gr.radius**2, gr.radius)
  return PhiThetaR_avg


#----------------get inside TC average--------------
def get_avg_inTC(data):
  integrand = in_TC_map*( N.trapz(data, dx=2*N.pi/gr.nphi, axis=0)/2./N.pi )
  Theta_avg = N.trapz(integrand, gr.colatitude, axis=0)  
  R_avg =  N.trapz(Theta_avg*gr.radius, gr.radius)
  R_avg = R_avg/area_in
  return R_avg

#----------------get outside TC average--------------
def get_avg_outTC(data):
  integrand = out_TC_map*( N.trapz(data, dx=2*N.pi/gr.nphi, axis=0)/2./N.pi )
  Theta_avg = N.trapz(integrand, gr.colatitude, axis=0)
  R_avg =  N.trapz(Theta_avg*gr.radius, gr.radius)
  R_avg =  R_avg/area_out
  return R_avg

#------------------Gradient-------------------------
def grad(data):
  r_compo = rderavg(data, eta=gr.radratio) 
  t_compo = thetaderavg(data)/rr3D
  p_compo = phideravg(data)/rr3D/N.sin(th3D)
  return r_compo, t_compo, p_compo

#-------------------Scalar Laplace Operator----------
def scal_laplace(data):
  result  = ((1/rr3D)**2) * rderavg( (rderavg(data, eta=gr.radratio)*rr3D**2), eta=gr.radratio) + \
  (1/N.sin(th3D)/rr3D**2) * thetaderavg( (N.sin(th3D) * (thetaderavg(data)))) + \
  (1/(N.sin(th3D)**2)/(rr3D**2)) * phideravg(phideravg(data))
  return result

#-------------------divergence operator--------------
def divergence(r_compo, t_compo, p_compo):
  div = rderavg(r_compo*rr3D**2, eta=gr.radratio)/rr3D**2 + \
        thetaderavg(N.sin(th3D)*t_compo)/rr3D/N.sin(th3D) + \
        phideravg(p_compo)/rr3D/N.sin(th3D)
  return div

#-------------------Curl of a vector-----------------
def curl(r_compo, t_compo, p_compo):
  curl_r =  (thetaderavg(p_compo*N.sin(th3D)) - phideravg(t_compo, gr.minc))/rr3D/N.sin(th3D)
  curl_t =  (phideravg(r_compo, gr.minc)/N.sin(th3D) - rderavg(rr3D*p_compo, eta=gr.radratio))/rr3D
  curl_p =  (rderavg(rr3D*t_compo, eta=gr.radratio) - thetaderavg(r_compo))/rr3D
  return curl_r, curl_t, curl_p

#-------------------Advective Derivative------------
def advec_der(r_compo, t_compo, p_compo):

  adv_r = gr.vr * rderavg(r_compo, eta=gr.radratio) + \
          gr.vtheta * thetaderavg(r_compo) / rr3D + \
          gr.vphi * phideravg(r_compo, gr.minc) / N.sin(th3D) / rr3D - \
          (gr.vtheta*t_compo + gr.vphi*p_compo) / rr3D

  adv_t = gr.vr * rderavg(t_compo, eta=gr.radratio) + \
          gr.vtheta * thetaderavg(t_compo) / rr3D + \
          gr.vphi * phideravg(t_compo, gr.minc) / N.sin(th3D) / rr3D + \
          gr.vtheta * r_compo / rr3D - \
          gr.vphi*p_compo * N.arctan(th3D) / rr3D

  adv_p = gr.vr * rderavg(p_compo, eta=gr.radratio) + \
          gr.vtheta * thetaderavg(p_compo) / rr3D + \
          gr.vphi * phideravg(p_compo, gr.minc) / N.sin(th3D) / rr3D+ \
          gr.vphi * r_compo / rr3D + \
          gr.vphi * t_compo * N.arctan(th3D) / rr3D

  return adv_r, adv_t, adv_p


#------------------combo averages---------------
def combo_avg(r_compo, t_compo, p_compo, inout=True):

 if inout==True:
  RMS_r = get_avg(r_compo**2)**0.5
  RMS_t = get_avg(t_compo**2)**0.5
  RMS_p = get_avg(p_compo**2)**0.5
  RMS = ((RMS_r**2) + \
         (RMS_t**2) + \
         (RMS_p**2))**0.5

  RMS_r_inTC = get_avg_inTC(r_compo**2)**0.5
  RMS_t_inTC = get_avg_inTC(t_compo**2)**0.5
  RMS_p_inTC = get_avg_inTC(p_compo**2)**0.5
  RMS_inTC = ((RMS_r_inTC**2) + \
  	      (RMS_t_inTC**2) + \
              (RMS_p_inTC**2))**0.5

  RMS_r_outTC = get_avg_outTC(r_compo**2)**0.5
  RMS_t_outTC = get_avg_outTC(t_compo**2)**0.5	
  RMS_p_outTC = get_avg_outTC(p_compo**2)**0.5
  RMS_outTC = ((RMS_r_outTC**2) + \
               (RMS_t_outTC**2) + \
	       (RMS_p_outTC**2))**0.5

  return RMS, RMS_inTC, RMS_outTC

 else:
  RMS_r = get_avg(r_compo**2)**0.5
  RMS_t = get_avg(t_compo**2)**0.5
  RMS_p = get_avg(p_compo**2)**0.5
  RMS = ((RMS_r**2) + \
         (RMS_t**2) + \
         (RMS_p**2))**0.5

  return RMS
         


##############################################################
#
#
#
#

#----------------------Lorentz force----------------------
# Lorentz force as combination of tension and pressure forces
if Lorentz==True:
  #Current as curl of B
  J_r, J_t, J_p =  curl(gr.Br, gr.Btheta, gr.Bphi)
  #Lorentz force as J cross B
  L_r = J_t*gr.Bphi - gr.Btheta*J_p
  L_t = -J_r*gr.Bphi + gr.Br*J_p
  L_p = J_r*gr.Btheta - gr.Br*J_t

  # convecrt to MagIC units
  L_r = L_r/gr.prmag/gr.ek
  L_t = L_t/gr.prmag/gr.ek
  L_p = L_p/gr.prmag/gr.ek

  RMS_L = combo_avg(L_r, L_t, L_p, inout=inout_TC)

  print 'Lorentz done'

#----------------------Coriolis force--------------------
if Coriolis==True:
  # as -2 zhat cross v

  C_r = 2*gr.vphi*N.sin(th3D)
  C_r = C_r/gr.ek

  C_t = 2*gr.vphi*N.cos(th3D)
  C_t = C_t/gr.ek

  C_p = -2*(gr.vr*N.sin(th3D) + gr.vtheta*N.cos(th3D))
  C_p = C_p/gr.ek

  RMS_C = combo_avg(C_r, C_t, C_p, inout=inout_TC)
  print 'Coriolis done'

#----------------------Pressure force----------------------
if Pressure==True:
	# as -grad(P)
  P_r, P_t, P_p = grad(gr.pre)
  P_r = -P_r; P_t = -P_t; P_p = -P_p
        
  RMS_P = combo_avg(P_r, P_t, P_p, inout=inout_TC)
  print 'Pressure done'

#--------------------inertial force-----------------------
if Inertial==True:
	
  I_r, I_t, I_p = advec_der(gr.vr, gr.vtheta, gr.vphi)
  I_r = -I_r; I_t = -I_t; I_p = -I_p
  RMS_I = combo_avg(I_r, I_t, I_p, inout=inout_TC)
  print 'Inertial done'

#---------------------Bouyancy Force-----------------------
if Buoyancy==True:
  ro = 1./(gr.radratio-1.)
  ri = gr.radratio/(gr.radratio-1.)

  bou = gr.ra/gr.pr*gr.entropy*rr3D/ro

  RMS_bou = get_avg(bou**2)**0.5
  RMS_bou_inTC = get_avg_inTC(bou**2)**0.5
  RMS_bou_outTC = get_avg_outTC(bou**2)**0.5
  print 'Buoyancy done'

#---------------------Viscous Force-----------------------
if Viscous==True:

  GD_r, GD_t, GD_p = grad(divergence(gr.vr, gr.vtheta, gr.vphi))
  Cur_r, Cur_t, Cur_p = curl(gr.vr, gr.vtheta, gr.vphi)
  CurCur_r, CurCur_t, CurCur_p = curl(Cur_r, Cur_t, Cur_p)

  V_r = GD_r - CurCur_r
  V_t = GD_t - CurCur_t
  V_p = GD_p - CurCur_p

  RMS_V = combo_avg(V_r, V_t, V_p, inout=inout_TC)
  print 'Viscous done'


#------------------------------------------------------------
if force_can_cal==True:
  CP_r = C_r + P_r
  CP_t = C_t + P_t
  CP_p = C_p + P_p
  RMS_CP = combo_avg(CP_r, CP_t, CP_p, inout=inout_TC)

  CPB_r = CP_r + bou
  RMS_CPB = combo_avg(CPB_r, CP_t, CP_p, inout=inout_TC)

  CPV_r = CP_r + V_r
  CPV_t = CP_t + V_t
  CPV_p = CP_p + V_p
  RMS_CPV = combo_avg(CPV_r, CPV_t, CPV_p, inout=inout_TC)

  CPVI_r = CPV_r + I_r 
  CPVI_t = CPV_t + I_t 
  CPVI_p = CPV_p + I_p
  RMS_CPVI = combo_avg(CPVI_r, CPVI_t, CPVI_p, inout=inout_TC)

  CVB_r = C_r + V_r + bou
  RMS_CVB = combo_avg(CVB_r, C_t+V_t, C_p+V_p, inout=inout_TC)

  CVI_r = C_r + V_r + I_r
  CVI_t = C_t + V_t + I_t 
  CVI_p = C_p + V_p + I_p
  RMS_CVI = combo_avg(CVI_r, CVI_t, CVI_p, inout=inout_TC)

  CIB_r = C_r + I_r + bou
  RMS_CIB = combo_avg(CIB_r, C_t+I_t, C_p+I_p, inout=inout_TC)


  LC_r = L_r + C_r
  LC_t = L_t + C_t
  LC_p = L_p + C_p 
  RMS_LC = combo_avg(LC_r, LC_t, LC_p, inout=inout_TC)

  LCP_r = LC_r + P_r
  LCP_t = LC_t + P_t
  LCP_p = LC_p + P_p
  RMS_LCP = combo_avg(LCP_r, LCP_t, LCP_p, inout=inout_TC)

  LCPB_r = LCP_r + bou
  RMS_LCPB = combo_avg(LCPB_r, LCP_t, LCP_p, inout=inout_TC)

  LCPBI_r = LCPB_r + I_r
  LCPBI_t = LCP_t + I_t
  LCPBI_p = LCP_p + I_p
  RMS_LCPBI = combo_avg(LCPBI_r, LCPBI_t, LCPBI_p, inout=inout_TC)

  LCPBIV_r = LCPBI_r + V_r
  LCPBIV_t = LCPBI_t + V_t
  LCPBIV_p = LCPBI_p + V_p
  RMS_LCPBIV = combo_avg(LCPBIV_r, LCPBIV_t, LCPBIV_p, inout=inout_TC)

  LCPV_r = LCP_r + V_r
  LCPV_t = LCP_t + V_t
  LCPV_p = LCP_p + V_p
  RMS_LCPV = combo_avg(LCPV_r, LCPV_t, LCPV_p, inout=inout_TC)


#-------------------------------------------------------------
if report==True:
  if ((inout_TC==True) & (force_can_cal==True)):
    print RMS_L[0], RMS_L[1], RMS_L[2], \
          RMS_C[0], RMS_C[1], RMS_C[2], \
          RMS_P[0], RMS_P[1], RMS_P[2], \
          RMS_bou, RMS_bou_inTC, RMS_bou_outTC, \
          RMS_I[0], RMS_I[1], RMS_I[2], \
          RMS_V[0], RMS_V[1], RMS_V[2], \
          RMS_CP[0], RMS_CP[1], RMS_CP[2], \
          RMS_CPB[0], RMS_CPB[1], RMS_CPB[2], \
          RMS_CPV[0], RMS_CPV[1], RMS_CPV[2], \
          RMS_CPVI[0], RMS_CPVI[1], RMS_CPVI[2], \
          RMS_CVB[0], RMS_CVB[1], RMS_CVB[2], \
          RMS_CVI[0], RMS_CVI[1], RMS_CVI[2], \
          RMS_CIB[0], RMS_CIB[1], RMS_CIB[2], \
          RMS_LC[0], RMS_LC[1], RMS_LC[2], \
          RMS_LCP[0], RMS_LCP[1], RMS_LCP[2], \
          RMS_LCPB[0], RMS_LCPB[1], RMS_LCPB[2], \
          RMS_LCPBI[0], RMS_LCPBI[1], RMS_LCPBI[2], \
          RMS_LCPBIV[0], RMS_LCPBIV[1], RMS_LCPBIV[2], \
          RMS_LCPV[0], RMS_LCPV[1], RMS_LCPV[2]
  elif ((inout_TC==True) & (force_can_cal==False)):
    print RMS_L[0], RMS_L[1], RMS_L[2], \
          RMS_C[0], RMS_C[1], RMS_C[2], \
          RMS_P[0], RMS_P[1], RMS_P[2], \
          RMS_bou, RMS_bou_inTC, RMS_bou_outTC, \
          RMS_I[0], RMS_I[1], RMS_I[2], \
          RMS_V[0], RMS_V[1], RMS_V[2]
  elif ((inout_TC==False) & (force_can_cal==False)):
    print RMS_L, RMS_C, RMS_P, RMS_bou, RMS_I, RMS_V
   
#--------------For plotting------------
if False:
  th = N.linspace(0, N.pi, gr.ntheta)
  rr, tth = N.meshgrid(gr.radius, th)
  xx = rr * N.sin(tth)
  yy = rr * N.cos(tth)

  #-----figure 1----
  data = (L_p).mean(axis=0)
  vmax =  data.max()/3
  vmin =  -vmax

  cs = N.linspace(vmin, vmax, 30)

  fig = P.figure(figsize=(4,5))
  ax = fig.add_axes([0.02, 0.02, 0.6, 0.88])
  im = ax.contourf(xx, yy, data, cs, cmap=P.get_cmap('RdYlBu_r'), extend='both')

  pos = ax.get_position()
  l, b, w, h = pos.bounds
  cax = fig.add_axes([0.75, 0.46-0.7*h/2., 0.03, 0.7*h])
  mir = fig.colorbar(im, cax=cax)
  ax.axis('off')

if False:
  #---figure 2-----
  data = (C_p).mean(axis=0)
  vmax =  data.max()/3
  vmin =  -vmax
  cs = N.linspace(vmin, vmax, 30)
  fig2 = P.figure(figsize=(4,5))
  ax2 = fig2.add_axes([0.02, 0.02, 0.6, 0.88])
  im2 = ax2.contourf(xx, yy, data, cs, cmap=P.get_cmap('RdYlBu_r'), extend='both')

  pos = ax2.get_position()
  l, b, w, h = pos.bounds
  cax = fig2.add_axes([0.75, 0.46-0.7*h/2., 0.03, 0.7*h])
  mir = fig2.colorbar(im2, cax=cax)
  ax2.axis('off')
	
if False:
#---figure 3-----
  data = (L_r).mean(axis=0)
  vmax =  data.max()/3
  vmin =  -vmax

  cs = N.linspace(vmin, vmax, 30)

  fig3 = P.figure(figsize=(4,5))
  ax3 = fig3.add_axes([0.02, 0.02, 0.6, 0.88])
  im3 = ax3.contourf(xx, yy, data, cs, cmap=P.get_cmap('RdYlBu_r'), extend='both')

  pos = ax3.get_position()
  l, b, w, h = pos.bounds
  cax = fig3.add_axes([0.75, 0.46-0.7*h/2., 0.03, 0.7*h])
  mir = fig3.colorbar(im3, cax=cax)
  ax3.axis('off')

P.show()

