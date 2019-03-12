from pylab import sqrt, loadtxt, zeros, ones, linspace, pi
#from magic import anelprof

class LD:
	
	def __init__(self, which_data='MB'):

		if ( which_data=='FSH' or which_data=='FSD' or \
			which_data=='FSM' or  which_data=='FSDyn' or \
			which_data=='NSH' or which_data=='NSD'
			):
			#-----------------------
			data_set1={'FSH':'free_hyd', 'FSD':'free_dip', 
						'FSM':'free_multi', 'FSDyn':'free_dyn', 
						'NSH':'no-slip_hyd', 'NSD':'no-slip_dyn'}
			data=loadtxt(data_set1[which_data])
			self.Ra=data[:, 0]
			self.Rac=ones(len(self.Ra))
			self.shell_vol=ones(len(self.Ra))*14.59
			self.E=data[:,1]
			self.Pr=data[:, 2]
			self.Pm=data[:, 3]
			self.Nrho = data[:, 4]
			self.ekin_pol = data[:, 5]
			self.ekin_tor = data[:, 6]
			self.ekin_pola = data[:, 7]
			self.ekin_tora = data[:, 8]
			self.emag_pol = data[:, 13]
			self.emag_tor = data[:, 14]
			self.emag_pola = data[:, 15]
			self.emag_tora = data[:, 16]
			self.Rm = data[:, 17]
			self.Rol = data[:, 19]
			self.Dip = data[:,21]
			self.DipCMB = data[:, 22]
			self.Els = data[:, 23]
			self.ElsCMB = data[:,24]
			self.Nu = data[:, 25]
			self.dlV = data[:,26]
			self.dlVc = data[:,41]
			self.dlB=data[:,28]
			self.dmB=data[:,29]
			self.diptotal=data[:,30]
			self.Dipl11=data[:,31]
			self.diptotl11 = data[:,32]
			self.dip3=data[:,33]
			self.buoyancy = data[:, 34]
			self.ohm_diss = data[:, 35]
			self.fohm = data[:, 36]
			#-------------------------------
			self.Numod=(self.Nu-1.)*(self.E/self.Pr)		
			self.RaQ=self.Ra*self.Numod*(self.E**2.)/self.Pr	
			#-------------------------------
			self.Ekin=self.ekin_pol+self.ekin_tor
			self.Emag=self.emag_pol+self.emag_tor
			self.Ekin_convec=self.Ekin-self.ekin_tora
			self.EkinENZ=self.Ekin/self.Ekin_convec
			#--------------------------------
			self.Ro=self.E*sqrt((2.*self.Ekin)/self.shell_vol)
			self.zonal_Ro=self.E*sqrt((2.*self.ekin_tora)/self.shell_vol)
			self.convec_Ro=self.E*sqrt((2.*self.Ekin_convec)/self.shell_vol)
			#--------------------------------
			if ( which_data=='FSD' or which_data=='FSM' or which_data=='NSD'):
				self.Els=(2.*self.Emag*(self.E*self.Pm))/(self.shell_vol)
				self.Els_pol=(2.*self.emag_pol*(self.E*self.Pm))/(self.shell_vol)
				self.Els_tor=(2.*self.emag_tor*(self.E*self.Pm))/(self.shell_vol)
				self.Lo=(self.Els*self.E/self.Pm)**0.5
				self.LoCMB=(self.ElsCMB*self.E/self.Pm)**0.5
				self.Lo_pol=(self.Els_pol*self.E/self.Pm)**0.5
				self.Lo_tor=(self.Els_tor*self.E/self.Pm)**0.5
				self.Lofohm=self.Lo/(self.fohm**0.5)
				self.Lofohm_pol=self.Lo_pol/(self.fohm**0.5)
				self.Lofohm_tor=self.Lo_tor/(self.fohm**0.5)
				self.tmag=(self.Emag)/(self.ohm_diss)
		elif ((which_data=='Lucia')):
			#------------------------
			data=loadtxt(which_data)
			self.Nrho = data[:,0]
			self.Ra = data[:,1]
			self.E = data[:,2]
			self.Pm = data[:,3]
			self.shell_vol=ones(len(self.Ra))*8.11
			self.Pr = ones(len(self.Ra))
			self.Rol= data[:,4]
			self.u_Rol = data[:,5]
			self.Nu_IC = data[:,6]
			self.Nu_OC = data[:,7]
			self.Nu = self.Nu_IC
			self.ekin_pol = data[:, 8]
			self.ekin_tor = data[:, 9]
			self.ekin_pola = data[:, 10]
			self.ekin_tora = data[:, 11]
			#-------------------------------
			self.Numod=(self.Nu-1.)*(self.E/self.Pr)		
			self.RaQ=self.Ra*self.Numod*(self.E**2.)/self.Pr	
			#-------------------------------
			self.Ekin=self.ekin_pol+self.ekin_tor
			self.Ekin_convec=self.Ekin-self.ekin_tora
			self.EkinENZ=self.Ekin/self.Ekin_convec
			#--------------------------------
			self.Ro=self.E*sqrt((2.*self.Ekin)/self.shell_vol)
			self.zonal_Ro=self.E*sqrt((2.*self.ekin_tora)/self.shell_vol)
			self.convec_Ro=self.E*sqrt((2.*self.Ekin_convec)/self.shell_vol)
		elif ((which_data=='Thomas_FSH') or (which_data=='Thomas_NSH') or \
				(which_data=='Thomas_FSDyn')):
			#------------------------
			data=loadtxt(which_data)
			self.Ra = data[:,0]
			self.E = data[:,1]
			self.Nrho = 0*data[:,2]
			self.Pr = data[:,3]
			self.shell_vol=ones(len(self.Ra))*51.3
			self.ekin_pol = data[:, 8]
			self.ekin_tor = data[:, 9]
			self.ekin_pola = data[:, 10]
			self.ekin_tora = data[:, 11]
			self.Rol = data[:, 15]
			self.Nu = data[:, 16]
			#-------------------------------
			self.Numod=(self.Nu-1.)*(self.E/self.Pr)		
			self.RaQ=self.Ra*self.Numod*(self.E**2.)/self.Pr	
			#-------------------------------
			self.Ekin=self.ekin_pol+self.ekin_tor
			self.Ekin_convec=self.Ekin-self.ekin_tora
			self.EkinENZ=self.Ekin/self.Ekin_convec
			#--------------------------------
			self.Ro=self.E*sqrt((2.*self.Ekin)/self.shell_vol)
			self.zonal_Ro=self.E*sqrt((2.*self.ekin_tora)/self.shell_vol)
			self.convec_Ro=self.E*sqrt((2.*self.Ekin_convec)/self.shell_vol)
