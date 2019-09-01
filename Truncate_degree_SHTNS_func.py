import shtns
import numpy as np


def shtns_truncate(data, l_max_trunc=20):

    # assumed input data shape => nphi,ntheta

    # make sure data is float64, requirement for shtns
    data = np.array(data, dtype=np.float64)

    #=========original grid specs=========
    l_max = int((data.shape[1]*2)/3)
    m_max = l_max
    lm_max = int(l_max*(l_max+1)/2 + l_max + 1)
    n_theta_max = data.shape[1]#(3*l_max + 1)/2
    print('Input ===> l_max={}, lm_max={}, n_theta_max={}'.format(l_max, lm_max, n_theta_max))
    #=====================================

    #Create transform object
    sh = shtns.sht(l_max, m_max, norm=shtns.sht_orthonormal | shtns.SHT_NO_CS_PHASE)
    #print('SHTNS lm_max = {}'.format(sh.nlm))

    polar_opt_threshold = 1e-10

    # setup grid for shtns
    nlat = int(l_max*(3./2./2.)*2.)
    nphi = 2*nlat
    nlat, nphi = sh.set_grid(nlat, nphi, polar_opt=polar_opt_threshold)
    #print('data shape of input grid needed: nlat={}, nphi={}'.format(nlat,nphi))


    # real space to specral
    datalm = sh.analys(data.T)
    #print('SHTNS spec data shape (lm_max): {}'.format(datalm.shape))


    #=========truncated grid specs=========
    m_max_trunc = l_max_trunc
    lm_max_trunc = int( l_max_trunc*(l_max_trunc+1)/2 + l_max_trunc + 1 )
    n_theta_max_trunc = int( (3*l_max_trunc + 1)/2 )
    print('[Truncated values] l_max={}, lm_max={}, n_theta_max={}'.format(l_max_trunc, lm_max_trunc, n_theta_max_trunc))
    #=======================================

    #define new transform object with the smaller grid
    sh_trunc = shtns.sht(l_max_trunc, m_max_trunc, mres=1, norm=shtns.sht_orthonormal | shtns.SHT_NO_CS_PHASE)
    #print('SHTNS truncated lm_max = {}'.format(sh_trunc.nlm))

    # set up truncated grid for shtns
    nlat = max(int(l_max_trunc*(3./2./2.)*2.),384)
    nphi = 2*nlat
    nlat, nphi = sh_trunc.set_grid(nlat, nphi, polar_opt=polar_opt_threshold)
    print('======SHTNS will output truncated data on: nlat={}, nphi={}'.format(nlat,nphi))

    # define empty array to be filled later
    datalm_trunc = np.zeros(lm_max_trunc, dtype=np.complex128)

    # get the l,m values from the main datalm
    # and assign to smaller lmax data
    for l in range(1, l_max_trunc+1):
     for m in range(0, l+1):
      #print(sh_trunc.idx(l,m), sh.idx(l,m))
      datalm_trunc[sh_trunc.idx(l,m)] = datalm[sh.idx(l,m)]

    # return from spectral to real space
    data_trunc = sh_trunc.synth(datalm_trunc)

    return data_trunc

"""
#### Example usage with Magic G file

from magic import MagicGraph
import matplotlib.pyplot as plt

G =  MagicGraph()

data = G.Br[:,:,5]
print('Magic input data: nphi={},ntheta={}'.format(data.shape[0], data.shape[1]))


data_trunc = shtns_truncate(data,l_max_trunc=50)

plt.imshow(data_trunc)
plt.show()

"""
