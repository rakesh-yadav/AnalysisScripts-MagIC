import numpy as np
import matplotlib.pyplot as plt

lmax=853
difnu=90
ldif=100
ldifexp=10
Ek=7e-7

degrees=np.arange(lmax)[1:]
mod_ek=Ek*(1+np.zeros_like(degrees))

for ind, degree in enumerate(degrees):
    #print(ind, degree)
    if degree>ldif:
        fac =1+ difnu*( ((degree+1-ldif)/(lmax+1-ldif))**ldifexp)
        mod_ek[ind] = mod_ek[ind]*fac


fig = plt.figure(figsize=(9,5.5))
ax = fig.add_axes([0.18, 0.15, 0.67, 0.8])
ax.tick_params(axis='both', which='major', labelsize=15)
ax.tick_params(which='both', width=1.5)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4)

ax.set_xlabel('$\ell$', fontsize=30)
ax.set_ylabel('Effective Ekman number', fontsize=20)
ax.loglog(degrees,mod_ek)
plt.show()
