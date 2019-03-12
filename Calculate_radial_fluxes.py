from __future__ import division
from magic import *
import pylab as P
import numpy as N
from stress_ten import ST
from glob import glob

#background = N.loadtxt('scond.dat')
#back_entropy = N.zeros_like(gr.entropy)
#for i in range(gr.nr):
# back_entropy[:, :, i] = background[i,1]

L_cond = []; L_conv = []
L_KE   = []; L_visc = []
L_poyn = []; L_res  = []

for g in glob("G*.tag*"):
 gr=MagicGraph(tag=g[-6:])

 S = ST(gr)
 
 Fr = (1+gr.radratio)/(2*(-1+N.exp(gr.strat/gr.polind)))
 Di = gr.pr/gr.ra/Fr
 mu_o = 4*N.pi*1e-7

 #phi_lin = N.linspace(-N.pi, N.pi, gr.nphi)
 #theta_lin = gr.colatitude
 th2D = N.zeros((gr.ntheta, gr.nr), 'f')
 #rr2D = N.zeros_like(th2D)
 for i in range(gr.ntheta):
  th2D[i, :] = gr.colatitude[i]
  #rr2D[i, :] = gr.radius

 #conductive flux
 temp0, rho0, beta0 = anelprof(gr.radius, gr.strat, gr.polind, gr.g0, gr.g1, gr.g2)
 L_cond_Phi = N.trapz((-rho0*temp0*rderavg(gr.entropy, eta=gr.radratio)/gr.pr), dx=2*N.pi/gr.nphi, axis=0)
 L_cond_PhiTheta = (gr.radius**2)*N.trapz(L_cond_Phi*N.sin(th2D), dx=N.pi/gr.ntheta, axis=0)
 L_cond.append(L_cond_PhiTheta)
 
 #convective flux
 L_conv_Phi = N.trapz((rho0*temp0*gr.vr*gr.entropy), dx=2*N.pi/gr.nphi, axis=0)
 L_conv_PhiTheta = (gr.radius**2)*N.trapz(L_conv_Phi*N.sin(th2D), dx=N.pi/gr.ntheta, axis=0)
 L_conv.append(L_conv_PhiTheta)
 
 #Kinetic energy flux
 L_KE_Phi = N.trapz((0.5*(-rho0*(gr.vr**2 + gr.vtheta**2 + gr.vphi**2)*gr.vr)), dx=2*N.pi/gr.nphi, axis=0)
 L_KE_PhiTheta = Di*(gr.radius**2)*N.trapz(L_KE_Phi*N.sin(th2D), dx=N.pi/gr.ntheta, axis=0)
 L_KE.append(L_KE_PhiTheta)

 #Viscous energy flux
 L_visc_Phi = N.trapz((-rho0*(gr.vr*S.rr + gr.vphi*S.rp + gr.vtheta*S.rt)), dx=2*N.pi/gr.nphi, axis=0)
 L_visc_PhiTheta = Di*(gr.radius**2)*N.trapz(L_visc_Phi*N.sin(th2D), dx=N.pi/gr.ntheta, axis=0)
 L_visc.append(L_visc_PhiTheta)

 #Poynting flux
 if gr.prmag!=0.0:
  L_poyn_PhiAvg =(-Di/gr.ek/gr.prmag)* \
                 ((gr.vr*gr.Br + gr.vphi*gr.Bphi + gr.vtheta*gr.Btheta)*gr.Br -\
                 (gr.Br**2 + gr.Bphi**2 + gr.Btheta**2)*gr.vr).mean(axis=0)
  L_poyn_PhiThetaAvg = 4*N.pi*(gr.radius**2)*L_poyn_PhiAvg.mean(axis=0)
  L_poyn.append(L_poyn_PhiThetaAvg)
 
 #Resistive flux
 if gr.prmag!=0.0:
  L_res_PhiAvg = (Di/gr.ek/gr.prmag/gr.prmag)* \
                 (-0.5*rderavg((gr.Br**2 + gr.Bphi**2 + gr.Btheta**2), eta=gr.radratio) + \
                 gr.Br*rderavg(gr.Br, eta=gr.radratio) + \
                 gr.Btheta*thetaderavg(gr.Br, order=4)/S.radius + \
                 gr.Bphi*phideravg(gr.Br)/S.radius/S.theta - \
                 (gr.Btheta**2 + gr.Bphi**2)/S.radius).mean(axis=0)
  L_res_PhiThetaAvg = 4*N.pi*(gr.radius**2)*L_res_PhiAvg.mean(axis=0)
  L_res.append(L_res_PhiThetaAvg)

 print g, 'is done.'

L_cond = (P.asarray(L_cond)).mean(axis=0)
L_conv = (P.asarray(L_conv)).mean(axis=0)
L_KE = (P.asarray(L_KE)).mean(axis=0)
L_visc = (P.asarray(L_visc)).mean(axis=0)
if len(L_poyn)!=0:
 L_poyn = (P.asarray(L_poyn)).mean(axis=0)
 L_res = (P.asarray(L_res)).mean(axis=0)

P.plot(gr.radius, L_cond)
P.plot(gr.radius, L_conv)
P.plot(gr.radius, L_KE)
P.plot(gr.radius, L_visc)
if len(L_poyn)!=0:
 P.plot(gr.radius, L_poyn)
 P.plot(gr.radius, L_res)

#---total----
if len(L_poyn)!=0:
 P.plot(gr.radius, L_cond+L_conv+L_KE+L_visc+L_poyn+L_res, '--k')
 P.legend((r'$L_{cond}$', r'$L_{conv}$', r'$L_{KE}$', r'$L_{visc}$', r'$L_{poyn}$', r'$L_{res}$', r'$L_{tot}$'), loc='best')
else:
 P.plot(gr.radius, L_cond+L_conv+L_KE+L_visc, '--k')
 P.legend((r'$L_{cond}$', r'$L_{conv}$', r'$L_{KE}$', r'$L_{visc}$', r'$L_{tot}$'), loc='best')
P.show()
