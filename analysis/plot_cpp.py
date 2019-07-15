#!/usr/bin/env python3
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
from scipy.optimize import curve_fit
from matplotlib.ticker import ScalarFormatter 

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

print('Beginning')
fname_list = np.loadtxt('list.dat',dtype=np.str)
Ntot = len(fname_list)
Array_list = []
beams = []
MHI = []
dist_MPC = []
DHI = []
Rs = []
RHI = []
vflat = []
alpha = []
mag =  []
ropt =  []
dx = []
RHI5 = []
Mstar = []
slope = []
for i in range(0,Ntot):
    Array_list.append(ascii.read(fname_list[i]))
    print('Loading Galaxy '+str(i+1)+' of '+str(Ntot))
    MHI_temp, DHI_temp,vflat_temp,beams_temp = np.loadtxt(fname_list[i],skiprows=1,usecols=(0,1,3,11),unpack=True)
    MHI.append( MHI_temp)
    DHI.append( DHI_temp)
    vflat.append( vflat_temp)
    beams.append( beams_temp)

N = len(Array_list[0])
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
VAnalytic = wfunc(Vrange,10**2.04,0.14,1.46,10**(-1.17))

#####################################################
# PLOTTING STUFF
#####################################################
lim = 4#np.float(sys.argv[1])

MHI_lim = []
vflat_lim = []
DHI_lim = []

inclination = []
rng = []

for i in range(0,Ntot):
    print(str(i+1)+" of "+str(len(fname_list)))
    MHI_lim.append(np.hstack(MHI[i][beams[i]>lim]))
    vflat_lim.append(np.hstack(vflat[i][beams[i] > lim]))
    DHI_lim.append(np.hstack(DHI[i][beams[i] > lim]))

    inclination.append(np.random.random_sample(size=np.shape(MHI_lim[i]))*90.)
    rng.append(np.random.random_sample(size=np.shape(MHI_lim[i]))*1.)

#####################################################################
import matplotlib as mpl
label_size=21.5
lw=1.5
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = lw
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = lw
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = lw
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.minor.width'] = lw
mpl.rcParams['axes.linewidth'] = lw
mpl.rcParams['font.monospace'] = 'Courier'
mpl.rcParams['legend.scatterpoints'] = '3'
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

a = np.min(MHI)
b = np.max(MHI)
num = 25
bins = np.linspace(a,b,num+1)
bw = (b-a)/num

av = np.min(np.log10(vflat))
bv = np.max(np.log10(vflat))
numv = 25
binsv = np.linspace(av,bv,numv+1)
bwv = (bv-av)/numv

V_tot = (4./3.)*np.pi*((c/H) * ((1.+z)**2.-1.)/((1.+z)**2.+1.))**3. #in Mpc

#####**********
zplt = []#np.array([])
yplt = []#np.array([])

fig, ax1 = plt.subplots(figsize=(16,10))
ax1.set_ylabel('Number Density\nn(V$_{rot}$)[Mpc$^{-3}$ dex$^{-1}$]',fontsize=label_size)
ax1.set_xlabel('Rotational Velocity, V$_{rot}$[log km s$^{-1}$]',fontsize=label_size)

for j in range(0,len(fname_list)):
    W = np.float64(np.ones_like(MHI[j])/(bw*V_tot))
    yplt_j,xplt_j = np.histogram(MHI[j],weights=W,bins=bins)
    zplt_j = yplt_j * 0
    print(str(j+1)+" of "+str(len(fname_list)))
    for i,x in enumerate(zplt_j):
        zplt_j[i] = (xplt_j[i]+xplt_j[i+1])/2.
    ax1.scatter(zplt_j,yplt_j,marker="+",zorder=3,color='C2',s=100)
    zplt.append(zplt_j)
    yplt.append(yplt_j)
ax1.plot(Mrange,Analytic,zorder=2,color='black')
ax1.set_yscale('log')
plt.tight_layout()
plt.savefig('MFUNC_universe.png',bbox_inches='tight')

#####**********
zplt = []#np.array([])
yplt = []#np.array([])

fig, ax2 = plt.subplots(figsize=(16,10))
ax2.set_ylabel('Number Density\nn(V$_{rot}$)[Mpc$^{-3}$ dex$^{-1}$]',fontsize=label_size)
ax2.set_xlabel('Rotational Velocity, V$_{rot}$[log km s$^{-1}$]',fontsize=label_size)

for j in range(0,len(fname_list)):
    V_max = (4./3.)*np.pi*((DHI_lim[j]/1000)*206265/(lim*30))**3. #in Mpc
    #for i,mass in enumerate(MHI_lim[j]):
    W = np.float64(np.ones_like(MHI_lim[j])/(bw*V_max))
    for i,mass in enumerate(W):
        if V_max[i] > V_tot:
            W[i] = 1./(bw*V_tot)


    yplt_j,xplt_j = np.histogram(MHI_lim[j],weights=W,bins=bins)
    zplt_j = yplt_j * 0
    
    print(str(j+1)+" of "+str(len(fname_list)))
    for i,x in enumerate(zplt_j):
        zplt_j[i] = (xplt_j[i]+xplt_j[i+1])/2.
    ax2.scatter(zplt_j,yplt_j,marker="+",zorder=3,color='C2',s=100)
    zplt.append(zplt_j)
    yplt.append(yplt_j)
ax2.plot(Mrange,Analytic,zorder=2,color='black')
ax2.set_yscale('log')
plt.tight_layout()
plt.savefig('MFUNC_observed.png',bbox_inches='tight')
###******************    
 
zplt = []#np.array([]) 
yplt = []#np.array([]) 
#fig, ax3 = plt.subplots(figsize=(16,10)) 
#ax3.set_ylabel('Number Density\nn(V$_{rot}$)[Mpc$^{-3}$ dex$^{-1}$]',fontsize=label_size) 
#ax3.set_xlabel('Rotational Velocity, V$_{rot}$[log km s$^{-1}$]',fontsize=label_size) 
 
for j in range(0,len(fname_list)): 
    W = np.float64(np.ones_like(np.log10(vflat[j]))/(bwv*V_tot)) 
    yplt_j,xplt_j = np.histogram(np.log10(vflat[j]),weights=W,bins=binsv) 
    zplt_j = yplt_j * 0 
    print(str(j+1)+" of "+str(len(fname_list))) 
    for i,x in enumerate(zplt_j): 
        zplt_j[i] = 10**((xplt_j[i]+xplt_j[i+1])/2.) 
    ax3.scatter(zplt_j,yplt_j,marker="+",zorder=3,color='C2',s=100) 
    zplt.append(zplt_j) 
    yplt.append(yplt_j) 
ax3.plot(Vrange,VAnalytic,zorder=2,color='black') 
ax3.set_yscale('log') 
ax3.set_xscale('log') 
plt.tight_layout() 
plt.savefig('VFUNC_universe.png',bbox_inches='tight') 
ax3.set_xticks([10,20,30,50,70,100,200,400]) 
ax3.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter()) 
ax3.set_xlim(10,500) 
ax3.set_ylim(1E-5,1E0)  
plt.savefig('VFUNC_universe.png',bbox_inches='tight')

###### Observed velocity function:

zplt = []#np.array([]) 
yplt = []#np.array([]) 
fig, ax4 = plt.subplots(figsize=(16,10)) 
ax4.set_ylabel('Number Density\nn(V$_{rot}$)[Mpc$^{-3}$ dex$^{-1}$]',fontsize=label_size) 
ax4.set_xlabel('Rotational Velocity, V$_{rot}$[log km s$^{-1}$]',fontsize=label_size) 

for j in range(0,len(fname_list)):
    V_max = (4./3.)*np.pi*((DHI_lim[j]/1000)*206265/(lim*30))**3. #in Mpc
    W = np.float64(np.ones_like(np.log10(vflat_lim[j]))/(bwv*V_max)) 
    for i,mass in enumerate(W):
        if V_max[i] > V_tot:
            W[i] = 1./(bw*V_tot)
    yplt_j,xplt_j = np.histogram(np.log10(vflat_lim[j]),weights=W,bins=binsv) 
    zplt_j = yplt_j * 0 
    print(str(j+1)+" of "+str(len(fname_list))) 
    for i,x in enumerate(zplt_j): 
        zplt_j[i] = 10**((xplt_j[i]+xplt_j[i+1])/2.) 
    zplt.append(zplt_j) 
    yplt.append(yplt_j) 

###**#*#*#*#**#* Money plot here:

minarray = np.array([])
medarray = np.array([])
iqarray = np.array([])
Narray = np.array([])
for j, y in enumerate(yplt[0]):
     tempmin = 10.
     jth_array = np.array([])
     for i in range(0,Ntot):
         if yplt[i][j] != 0:
             tempmin = np.min((yplt[i][j],tempmin))
             jth_array = np.append(jth_array,yplt[i][j])
     iqarray = np.append(iqarray,np.subtract(*np.percentile(jth_array, [75, 25])))
     Narray  = np.append(Narray,len(jth_array))
     medarray = np.append(medarray,np.median(jth_array))
     minarray = np.append(minarray,tempmin)

maxarray = np.array([])
for j, y in enumerate(yplt[0]):
     tempmax = 0
     for i in range(0,Ntot):
         if yplt[i][j] != 0:
             tempmax = np.max((yplt[i][j],tempmax))
     maxarray = np.append(maxarray,tempmax)

zplt_test = zplt[0]

plt.figure(1,figsize=(20,10))
gs = gridspec.GridSpec(4, 1)
ax4 = plt.subplot(gs[0:3,:])

ax4.plot(zplt_test[Narray == max(Narray)],minarray[Narray == max(Narray)],'--',color='C0')
ax4.plot(zplt_test,maxarray, '--',color='C0')
ax4.fill_between(zplt_test,medarray+iqarray,medarray,alpha=0.25,color='C2')
ax4.fill_between(zplt_test,medarray-iqarray,medarray,alpha=0.25,color='C2')
ax4.plot(zplt_test,medarray,color='C2')
ax4.plot(Vrange,VAnalytic,color='black')
ax4.set_ylabel('Number Density\nn(V$_{rot}$)[Mpc$^{-3}$ dex$^{-1}$]',fontsize=label_size) 
ax4.set_xlabel('Rotational Velocity, V$_{rot}$[log km s$^{-1}$]',fontsize=label_size) 
ax4.set_yscale('log')
ax4.set_xscale('log')

"""
ax4.plot(Vrange,wfunc(Vrange,1.73831153e+02, -7.09819856e-01,  2.01825507e+00,  2.95476876e-02),'--',alpha=0.25)
ax4.plot(Vrange,wfunc(Vrange,8.61583702e+01, 7.00729651e-01, 1.32571312e+00, 7.27689129e-02),'--',alpha=0.25)c
"""
ax4.set_ylim(1E-5,1E0)

ax4.set_xticks([10,20,30,50,70,100,200,400])
ax4.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax5 = plt.subplot(gs[3:,:])
ax5.set_xlabel('V$_{rot}$ [log10(km h$^{-1})$]',fontsize=label_size)
ax5.set_ylabel('N bins $\\ne$ 0',fontsize=label_size)
plt.tight_layout()

log_zplt = np.log10(zplt[0])
ax5.bar(log_zplt,Narray,color='white',edgecolor='black',width=log_zplt[1]-log_zplt[0])
ax5.set_xticks([])
ax5.set_yscale('log')
#ax5.set_xscale('log')

plt.savefig('Money.png')
plt.show()

#####################################################
# PLOTTING STUFF
#####################################################

