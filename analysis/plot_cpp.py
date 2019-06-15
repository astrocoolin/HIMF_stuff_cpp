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

fname_list = np.loadtxt('list_3.dat',dtype=np.str)
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

    beams.append( np.hstack(Array_list[i]['beams']))
    MHI.append( np.hstack(Array_list[i]['MHI']))
    #dist_MPC.append( np.hstack(Array_list[i]['MHI']))
    DHI.append( np.hstack(Array_list[i]['DHI']))
    #Rs.append( (DHI[i]/2.)*0.18)
    #RHI.append( DHI[i]/2.)
    vflat.append( np.hstack(Array_list[i]['vflat']))
    #alpha.append( np.hstack(Array_list[i]['alpha']))
    #mag.append( np.hstack(Array_list[i]['Mag']))
    #ropt.append( np.hstack(Array_list[i]['Ropt']))
    #dx.append( np.hstack(Array_list[i]['dx']))
    #RHI5.append( np.hstack(Array_list[i]['RHI5']))
    #Mstar.append( np.hstack(Array_list[i]['Mstar']))
    #slope.append( np.hstack(Array_list[i]['slope']))

    #vtest_RHI = polyex(mag,alpha,ropt,RHI)
    #vtest_OPT = polyex(mag,alpha,ropt,ropt)
    #vtest_RHI5 = polyex(mag,alpha,ropt,RHI5)

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
'''
plt.figure(1)

ax1 = plt.subplot()
ax1.set_ylabel('N')

lims = limcalc(MHI)

for i in range(0,len(fname_list)): 
    yplt,xplt = np.histogram(MHI[i],weights=lims[0][i],bins=lims[1]) 
    zplt = yplt * 0 
    for j,x in enumerate(zplt): 
        zplt[j] = (xplt[j]+xplt[j+1])/2. 
    ax1.scatter(zplt,yplt,marker="+")

ax1.plot(Mrange,Analytic,color='C1')
ax1.set_yscale('log')
ax1.set_xlim(lims[2],lims[3])
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

for i in range(0,len(fname_list)):
    yplt,xplt = np.histogram(np.log10(vflat[i]),weights=lims[0][i],bins=lims[1])           
    zplt = yplt * 0
    print(str(i+1)+" of "+str(len(fname_list)))
    for j,x in enumerate(zplt):
        zplt[j] = (xplt[j]+xplt[j+1])/2.
    ax1.scatter(zplt,yplt,marker="+")

ax1.plot(np.log10(Vrange),VAnalytic,color='C1',label='Analytic')
ax1.set_yscale('log')

#####################################################################
plt.figure(3)
'''
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

'''
a = np.min(MHI)
b = np.max(MHI)

num = 50

bins = np.linspace(a,b,num+1)

bw = (b-a)/num

ax1 = plt.subplot()
ax1.set_ylabel('N')

ax1.scatter(zplt,yplt,color='C3',marker="+")
ax1.plot(Mrange,Analytic,color='C1')
ax1.set_yscale('log')
ax1.set_xlim(7,11)


plt.tight_layout()
plt.savefig('HIMF_observed.png',bbox_inches='tight')

av = np.min(np.log10(vflat))
bv = np.max(np.log10(vflat))
numv = 100
binsv = np.linspace(av,bv,num+1)
bwv = (bv-av)/numv
'''
#####################################################################

plt.figure(4)
ax3 = plt.subplot()
ax3.set_ylabel('N')

av = np.min(np.log10(vflat))
bv = np.max(np.log10(vflat))
numv = 25
binsv = np.linspace(av,bv,numv+1)
bwv = (bv-av)/numv

V_tot = (4./3.)*np.pi*((c/H) * ((1.+z)**2.-1.)/((1.+z)**2.+1.))**3. #in Mpc
W = []
zplt = []#np.array([])
yplt = []#np.array([])

for j in range(0,len(fname_list)):
    V_max = (4./3.)*np.pi*((DHI_lim[j]/1000)*206265/(lim*30))**3. #in Mpc
    W.append(1./(bwv*V_max))
    for i,mass in enumerate(MHI_lim[j]):
        if V_max[i] > V_tot:
            W[j][i] = 1./(bwv*V_tot)

    yplt_j,xplt_j = np.histogram(np.log10(vflat_lim[j]),weights=W[j],bins=binsv)
    zplt_j = yplt_j * 0
    
    print(str(j+1)+" of "+str(len(fname_list)))
    for i,x in enumerate(zplt_j):
        zplt_j[i] = (xplt_j[i]+xplt_j[i+1])/2.
    ax3.scatter(zplt_j,yplt_j,marker="+")
    zplt.append(zplt_j)
    yplt.append(yplt_j)
    
ax3.plot(np.log10(Vrange),VAnalytic,color='C1')
ax3.set_yscale('log')
plt.tight_layout()
plt.savefig('VFUNC_observed.png',bbox_inches='tight')

plt.show()


from matplotlib.ticker import ScalarFormatter


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
     #print(minarray[j],zplt[0][j])

maxarray = np.array([])
for j, y in enumerate(yplt[0]):
     tempmax = 0
     for i in range(0,Ntot):
         if yplt[i][j] != 0:
             tempmax = np.max((yplt[i][j],tempmax))
     maxarray = np.append(maxarray,tempmax)

zplt_test = 10**zplt[0]

plt.figure(1,figsize=(20,10))
gs = gridspec.GridSpec(4, 1)
ax1 = plt.subplot(gs[0:3,:])

ax1.plot(zplt_test[-18:],minarray[-18:],'--')
#plt.plot(zplt[0],minarray, '--')
ax1.plot(zplt_test,maxarray, '--')
#plt.fill_between(zplt[0],maxarray,medarray,alpha=0.25,color='C0')
#plt.fill_between(zplt[0],minarray,medarray,where=(zplt[0]>1.44387287),alpha=0.25,color='C0')
ax1.fill_between(zplt_test,medarray+iqarray,medarray,alpha=0.25,color='C9')
ax1.fill_between(zplt_test,medarray-iqarray,medarray,alpha=0.25,color='C9')
#plt.fill_between(zplt[0],minarray,medarray,alpha=0.5)
ax1.plot(zplt_test,medarray,color='C0')
ax1.plot(Vrange,VAnalytic,color='C3')
ax1.set_xlabel('V$_{rot}$ [log10(km h$^{-1})$]')
ax1.set_ylabel('n(V$_{rot}$) (Mpc$^{-3}$ dex$^{-1}$ ')
ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.plot(Vrange,wfunc(Vrange,1.73831153e+02, -7.09819856e-01,  2.01825507e+00,  2.95476876e-02),'--',alpha=0.25)
ax1.plot(Vrange,wfunc(Vrange,8.61583702e+01, 7.00729651e-01, 1.32571312e+00, 7.27689129e-02),'--',alpha=0.25)

ax1.set_ylim(1E-5,1E0)

ax1.set_xticks([10,20,30,50,70,100,200,400])
ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax2 = plt.subplot(gs[3:,:])
ax2.set_xlabel('V$_{rot}$ [log10(km h$^{-1})$]')
ax2.set_ylabel('N bins $\\ne$ 0')
plt.tight_layout()
#ax2.set_yscale('log')
ax2.bar(zplt[0],Narray,color='white',edgecolor='black',width=zplt[0][1]-zplt[0][0])
#ax2.set_xscale('log')
#ax2.set_xticks([10,20,30,50,70,100,200,400])
ax2.set_xticks([])
#ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

plt.savefig('Money.png')
plt.show()

plt.figure(2)
plt.plot(zplt_test,maxarray)
plt.plot(zplt_test,medarray)
plt.plot(zplt_test,minarray)
plt.scatter(zplt_test[-30:],minarray[-30:])
plt.yscale('log') ; plt.xscale('log')
#plt.plot(Vrange,wfunc(Vrange,1.88720516e+02, -8.10314251e-01,  3.43987351e+00,2.17682759e-02))
plt.plot(Vrange,wfunc(Vrange,1.73831153e+02, -7.09819856e-01,  2.01825507e+00,  2.95476876e-02))
#plt.plot(Vrange,wfunc(Vrange,1.1562869 , 2.09399563, 0.94451162, 2.09412864))
plt.plot(Vrange,wfunc(Vrange,8.61583702e+01, 7.00729651e-01, 1.32571312e+00, 7.27689129e-02))
plt.ylim(np.min(medarray),np.max(maxarray))
plt.show()
#####################################################
# PLOTTING STUFF
#####################################################

