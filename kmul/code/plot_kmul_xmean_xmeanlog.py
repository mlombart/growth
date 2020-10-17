#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sympy as sym
import numpy as np
from scipy.special import legendre
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset


#change the current directory to the directory
#where the running script file (.py) exists
os.chdir(os.path.dirname(os.path.abspath(__file__)))
path1 = os.getcwd()
# print(path1)
# print(os.path.realpath('..'))
path2 = os.path.realpath('..')


#options for plots
plt.rcParams["font.size"]= 16
plt.rcParams['lines.linewidth'] = 3
plt.rcParams["legend.columnspacing"] = 0.5

marker_style = dict( marker='o',markersize=12, markerfacecolor='white', linestyle='',markeredgewidth=2)
savefig_options=dict(bbox_inches='tight')


#legendre polynomials
def LegendreP(i,x):
   coeffs=legendre(i)
   res = 0
   for j in range(0,i+1):
      res = res+ coeffs[j]*x**j

   return res

#ghtilde
def ghtilde(massgrid,massbins,gij,gh_mean,theta,k,j,x):
   xij = 2/(massgrid[j+1]-massgrid[j])*(x-massbins[j])
   res1 = 0 
   res2 = 0
   if k==0:
      res2 = theta[j]*(gij[j]*LegendreP(k,xij)-gh_mean[j])+gh_mean[j]
   else:
      for i in range(k+1):
         res1 = res1+ gij[j,i]*LegendreP(i,xij)   

      res2 = theta[j]*(res1-gh_mean[j])+gh_mean[j]
   return res2


def I(massgrid,j):
   res= np.logspace(np.log10(massgrid[j]),np.log10(massgrid[j+1]),num=100)
   return res


data_solkmulDL_tend = np.genfromtxt(path2+'/data/list_solkmul_tau=tend.csv',delimiter=",")



nbins=20

#data k0
################################
k0=0
path_data_k0 = path2+'/data/simu/dp/nbins='+str(nbins)+'/kmax='+str(k0)+'/'
massgridk0 = np.loadtxt(path_data_k0+'grainmassgrid.txt')
massbinsk0 = np.loadtxt(path_data_k0+'grainmassbins.txt')
massbinsk0meanlog = [np.sqrt(massgridk0[i]*massgridk0[i+1]) for i in range(nbins)]
gjt_xmean_k0 = np.loadtxt(path_data_k0+'gjt_xmean.txt')
gjt_xmeanlog_k0 = np.loadtxt(path_data_k0+'gjt_xmeanlog.txt')
timek0 = np.loadtxt(path_data_k0+'time.txt')
################################

#data k1
################################
k1=1
path_data_k1 = path2+'/data/simu/dp/nbins='+str(nbins)+'/kmax='+str(k1)+'/'
massgridk1 = np.loadtxt(path_data_k1+'grainmassgrid.txt')
massbinsk1 = np.loadtxt(path_data_k1+'grainmassbins.txt')
massbinsk1meanlog = [np.sqrt(massgridk1[i]*massgridk1[i+1]) for i in range(nbins)]
gjt_xmean_k1 = np.loadtxt(path_data_k1+'gjt_xmean.txt')
gjt_xmeanlog_k1 = np.loadtxt(path_data_k1+'gjt_xmeanlog.txt')
timek1 = np.loadtxt(path_data_k1+'time.txt')
################################

#data k2
################################
k2=2
path_data_k2 = path2+'/data/simu/dp/nbins='+str(nbins)+'/kmax='+str(k2)+'/'
massgridk2 = np.loadtxt(path_data_k2+'grainmassgrid.txt')
massbinsk2 = np.loadtxt(path_data_k2+'grainmassbins.txt')
massbinsk2meanlog = [np.sqrt(massgridk2[i]*massgridk2[i+1]) for i in range(nbins)]
gjt_xmean_k2 = np.loadtxt(path_data_k2+'gjt_xmean.txt')
gjt_xmeanlog_k2 = np.loadtxt(path_data_k2+'gjt_xmeanlog.txt')
timek2 = np.loadtxt(path_data_k2+'time.txt')
################################

#data k3
################################
k3=3
path_data_k3 = path2+'/data/simu/qp/nbins='+str(nbins)+'/kmax='+str(k3)+'/'
massgridk3 = np.loadtxt(path_data_k3+'grainmassgrid.txt')
massbinsk3 = np.loadtxt(path_data_k3+'grainmassbins.txt')
massbinsk3meanlog = [np.sqrt(massgridk3[i]*massgridk3[i+1]) for i in range(nbins)]
gjt_xmean_k3 = np.loadtxt(path_data_k3+'gjt_xmean.txt')
gjt_xmeanlog_k3 = np.loadtxt(path_data_k3+'gjt_xmeanlog.txt')
timek3 = np.loadtxt(path_data_k3+'time.txt')
################################

xmin = np.float(massgridk0[0])
xmax = np.float(massgridk0[-1])

yminloglog = 10**(-10)
ymaxloglog = 10**3

x=np.logspace(np.log10(xmin),np.log10(xmax),num=1000)

#plot gj_xmeanlog
fig,ax = plt.subplots()
ax.set_ylim(yminloglog,ymaxloglog)
ax.set_xlim(xmin,xmax)
ax.set_xscale('log')
ax.set_yscale('log')
ax.plot(data_solkmulDL_tend[:,0],data_solkmulDL_tend[:,1],'--',c='C0',label='Analytic')
ax.plot(massbinsk0meanlog,gjt_xmeanlog_k0[-1,:],markeredgecolor='black',label=r'$k=0$',**marker_style)
ax.plot(massbinsk1meanlog,gjt_xmeanlog_k1[-1,:],markeredgecolor='C3',label=r'$k=1$',**marker_style)
ax.plot(massbinsk2meanlog,gjt_xmeanlog_k2[-1,:],markeredgecolor='C1',label=r'$k=2$',**marker_style)
ax.plot(massbinsk3meanlog,gjt_xmeanlog_k3[-1,:],markeredgecolor='C2',label=r'$k=3$',**marker_style)
ax.set_xlabel(r'mass $x$')
ax.set_ylabel(r'mass density $g(x,\tau)$')
ax.set_title(r'$\tau=%d$' %(timek3[-1]))
ax.legend(loc='lower left',ncol=2)

# plt.savefig(path2+'/plots/kmul_tend_loglog_xmeanlog.pdf',**savefig_options)
# plt.savefig(path2+'/plots/kmul_tend_loglog_xmeanlog.png',**savefig_options)


plt.show()
