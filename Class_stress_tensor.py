from magic import *
import numpy as N


class ST:
	
 def __init__(self, gr):

  self.radius = N.zeros((gr.npI, gr.ntheta, gr.nr), dtype='f')
  self.theta = N.zeros((gr.npI, gr.ntheta, gr.nr), dtype='f')
  for i in range(gr.ntheta):
   self.theta[:, i, :] = gr.colatitude[i]
  for i in range(gr.nr): 
   self.radius[:, :, i] = gr.radius[i]

  div_v = rderavg(gr.vr*self.radius**2, eta=gr.radratio)/self.radius**2 + \
          thetaderavg(N.sin(self.theta)*gr.vtheta, order=4)/self.radius/N.sin(self.theta) + \
          phideravg(gr.vphi)/self.radius/N.sin(self.theta)
  #--diagonal components---------
  self.rr = 2.0*rderavg(gr.vr, eta=gr.radratio) - \
            2.0*div_v/3.0
  #self.tt = (thetaderavg(gr.vtheta, order=4) + \
  #          gr.vr)*2.0/self.radius - \
  #          2.0*div_v/3.0
  #self.pp = (phideravg(gr.vphi)/N.sin(self.theta) + \
  #          gr.vr + \
  #          gr.vtheta/N.tan(self.theta))*2.0/self.radius - \
  #          2.0*div_v/3.0
  
  #---off diaigonal components---
  self.rt = thetaderavg(gr.vr, order=4)/self.radius + \
            self.radius*rderavg(gr.vtheta/self.radius, eta=gr.radratio)
  self.rp = phideravg(gr.vr)/self.radius/N.sin(self.theta) + \
            self.radius*rderavg(gr.vphi/self.radius, eta=gr.radratio)
  #self.tp = phideravg(gr.vtheta)/self.radius/N.sin(self.theta) + \
  #          N.sin(self.theta)*thetaderavg(gr.vphi/N.sin(self.theta), order=4)/self.radius

