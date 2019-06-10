#!/usr/bin/env python3
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
from scipy.optimize import curve_fit

import sys

def phi(MHI, Mstar, alpha, phi_0):
    #Mass Function
    return np.log(10.)*phi_0* (MHI/Mstar)**(alpha+1.) * np.exp(-MHI/Mstar)

def wfunc(W,Wstar,alpha,beta,nstar):
    func = np.log(10) * nstar * (W/Wstar)**alpha * np.exp(-(W/Wstar)**beta)
    return func

def func_V0(mag):
    a = 2.25799719
    b = 0.40588194
    c = 4.94486961
    d = -1.88611724
    x=mag
    return a*np.exp(-x*b-c)+d*x

def func_rt(mag):
    a = 0.04264749
    b = 1.12975029
    x=mag
    return a*x +b

def polyex(inmag,alpha,Ropt,r):
    V0 = func_V0(inmag)
    rt = Ropt*func_rt(inmag)
    func = V0*(1.0-np.exp(-r/rt))*(1.0+alpha*r/rt)
    return func

Glist = ascii.read('Glist.txt')
print('Loaded Galaxies')

beams = Glist['beams']
MHI = Glist['MHI']
dist_MPC = Glist['MHI']
DHI = Glist['DHI']
Rs = (DHI/2.)*0.18
RHI = DHI/2.
vflat = Glist['vflat']
alpha = Glist['alpha']
mag = Glist['Mag']
ropt = Glist['Ropt']
dx = Glist['dx']
RHI5 = Glist['RHI5']
Mstar = Glist['Mstar']
slope = Glist['slope']

vtest_RHI = polyex(mag,alpha,ropt,RHI)
vtest_OPT = polyex(mag,alpha,ropt,ropt)
vtest_RHI5 = polyex(mag,alpha,ropt,RHI5)

N = len(Glist)
z = 0.03
c = 3E5 #km/s
H = 70 #km/s

def limcalc(array):
    a = np.min(array)
    b = np.max(array)
    num = 50
    bins = np.linspace(a,b,num+1)
    bw = (b-a)/num

    V_tot = (4./3.)*np.pi*((c/H) * ((1.+z)**2.-1.)/((1.+z)**2.+1.))**3. #in Mpc

    W = 1./(bw*V_tot)*np.ones_like(MHI)
    print('Bin Width [dex]',bw)
    return(W,bins,a,b)

Mrange = np.linspace(7.0,10.8462,100)
Analytic = phi(10.**Mrange,10.**9.96,-1.33,0.0048)

Vrange = np.linspace(10**0.9761,10**2.8,1000)
VAnalytic = wfunc(Vrange,10**2.06,0.14,1.46,10**(-1.17))

#####################################################
# PLOTTING STUFF
#####################################################

plt.figure(1)

gs = gridspec.GridSpec(4, 1)
ax1 = plt.subplot(gs[0:3,:])
ax1.set_ylabel('N')

lims = limcalc(MHI)

yplt,xplt = np.histogram(MHI,weights=lims[0],bins=lims[1]) 
zplt = yplt * 0
for i,x in enumerate(zplt):
    zplt[i] = (xplt[i]+xplt[i+1])/2.
ax1.hist(MHI,weights=lims[0],histtype='step',bins=lims[1])

ax1.scatter(zplt,yplt,color='C3',marker="+")
ax1.plot(Mrange,Analytic,color='C1')
ax1.set_yscale('log')
ax1.set_xlim(lims[2],lims[3])

ax2 = plt.subplot(gs[3:,:])
ax2.hist(MHI,histtype='step',color='black',bins=lims[1])
ax2.set_yscale('log')
ax2.set_ylabel('N')
ax2.set_xlabel('HI Mass')
ax2.set_xlim(lims[2],lims[3])

plt.tight_layout()
plt.xlabel('HI mass')
plt.savefig('HIMF_volume.png',bbox_inches='tight')
#####################################################################

plt.figure(2)

gs = gridspec.GridSpec(4, 1)
ax1 = plt.subplot(gs[0:3,:])
ax1.set_ylabel('N')
ax1.set_xlabel('Vcirc')

lims = limcalc(np.log10(vflat))

yplt,xplt = np.histogram(np.log10(vflat),weights=lims[0],bins=lims[1])           
zplt = yplt * 0
for i,x in enumerate(zplt):
    zplt[i] = (xplt[i]+xplt[i+1])/2.

ax1.scatter(zplt,yplt,color='C3',marker="+")

ax1.hist(np.log10(vflat),weights=lims[0],histtype='step',bins=lims[1],color='C4',label='0.5 Msun/Pc$^2$')
ax1.plot(np.log10(Vrange),VAnalytic,color='C1',label='Analytic')
ax1.set_yscale('log')
#ax1.set_xlim(lims[2],lims[3])
#ax1.set_xlim(1,3.0)
#ax1.legend()

ax2 = plt.subplot(gs[3:,:])
#ax2.hist(np.log10(vtest_RHI5),histtype='step',bins=lims[1])
ax2.hist(np.log10(vflat),histtype='step',bins=lims[1])
ax2.set_yscale('log')
ax2.set_ylabel('N')
ax2.set_xlabel('Vcirc')
plt.tight_layout()
plt.savefig('Velocity_volume.png',bbox_inches='tight')
#####################################################################

plt.figure(3)

lim = np.float(sys.argv[1])

#dist_MPC_lim = np.random.shuffle(dist_MPC[beams>lim])
#beams_lim = beams[beams>lim]


dist_MPC_lim = np.hstack(dist_MPC[beams>lim])
dist_MPC_lim = np.random.shuffle(dist_MPC_lim)

MHI_lim = np.hstack(MHI[beams>lim])
vflat_lim = np.hstack(vflat[beams > lim])
DHI_lim = np.hstack(DHI[beams > lim])

inclination = np.random.random_sample(size=np.shape(MHI_lim))*90.
rng = np.random.random_sample(size=np.shape(MHI_lim))*1.


a = np.min(MHI_lim)
av = np.min(np.log10(vflat))
b = np.max(MHI_lim)
bv = np.max(np.log10(vflat))

num = 50
numv = 50

bins = np.linspace(a,b,num+1)
binsv = np.linspace(av,bv,num+1)

bw = (b-a)/num
bwv = (bv-av)/numv
if True:
    inc_frac = inclination * 2. / 100.
    inc_frac[inc_frac > 1.] = 1.
    
    MHI_lim = MHI_lim[inc_frac >= rng]
    vflat_lim = vflat_lim[inc_frac >= rng]
    DHI_lim = DHI_lim[inc_frac >= rng]
    inc_frac = inc_frac[inc_frac >= rng]
else:
    inc_frac = np.ones_like(MHI_lim)

V_max = (4./3.)*np.pi*((DHI_lim/1000)*206265/(lim*30))**3. #in Mpc
V_tot = (4./3.)*np.pi*((c/H) * ((1.+z)**2.-1.)/((1.+z)**2.+1.))**3. #in Mpc

W = 1./(bw*V_max*inc_frac)
for i,mass in enumerate(MHI_lim):
    if V_max[i] > V_tot:
        W[i] = 1./(bw*inc_frac[i]*V_tot)

plt.figure(3)
gs = gridspec.GridSpec(4, 1)
ax1 = plt.subplot(gs[0:3,:])
ax1.set_ylabel('N')

lims = limcalc(MHI_lim)

yplt,xplt = yplt,xplt= np.histogram(MHI_lim,weights=W,bins=lims[1]) 
zplt = yplt * 0
for i,x in enumerate(zplt):
    zplt[i] = (xplt[i]+xplt[i+1])/2.

#ax1.axvline(x=9.3566)
ax1.scatter(zplt,yplt,color='C3',marker="+")
#ax1.hist(MHI_lim,weights=W,histtype='step',bins=lims[1])
ax1.plot(Mrange,Analytic,color='C1')
ax1.set_yscale('log')
ax1.set_xlim(7,11)


ax2 = plt.subplot(gs[3:,:])
ax2.hist(MHI_lim,histtype='step',color='black',bins=lims[1])
ax2.set_yscale('log')
ax2.set_ylabel('N')
ax2.set_xlabel('HI Mass')
ax2.set_xlim(7,11)
plt.tight_layout()
plt.savefig('HIMF_observed.png',bbox_inches='tight')

#W = 1./(bw*V_max*inc_frac)
#for i,mass in enumerate(MHI_lim):
#    if V_max[i] > V_tot:
#       W[i] = 1./(bw*inc_frac[i]*V_tot)

W = 1./(bwv*V_max*inc_frac)
for i,mass in enumerate(MHI_lim):
    if V_max[i] > V_tot:
        W[i] = 1./(bwv*inc_frac[i]*V_tot)

plt.figure(4)
gs = gridspec.GridSpec(4, 1)
ax3 = plt.subplot(gs[0:3,:])
ax3.set_ylabel('N')

lims = limcalc(np.log10(vflat))
#ax1.axvline(x=9.3566)
yplt,xplt = yplt,xplt= np.histogram(np.log10(vflat_lim),weights=W,bins=lims[1]) 
zplt = yplt * 0
for i,x in enumerate(zplt):
    #print(i,x,zplt[i],xplt[i])
    zplt[i] = (xplt[i]+xplt[i+1])/2.


#ax3.hist(np.log10(vflat_lim),weights=W,histtype='step',bins=lims[1])
#popt,pcov = curve_fit(vfunc,zplt,yplt)
ax3.scatter(zplt,yplt,color='C3',marker="+")
ax3.plot(np.log10(Vrange),VAnalytic,color='C1')
ax3.set_yscale('log')
#ax3.set_xlim(7,11)

ax4 = plt.subplot(gs[3:,:])
ax4.hist(np.log10(vflat_lim),histtype='step',color='black',bins=lims[1])
ax4.set_yscale('log')
ax4.set_ylabel('N')
ax4.set_xlabel('Vflat')
#ax4.set_xlim(7,11)
plt.tight_layout()
plt.savefig('VFUNC_observed.png',bbox_inches='tight')

plt.show()

#####################################################
# PLOTTING STUFF
#####################################################

