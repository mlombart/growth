#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sympy as sym
import numpy as np
from scipy.special import legendre,iv
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


#solution kconst g(x,0)=x exp(-x)
def solKconstDL(x,tau):
   res = 4.*x/((2+tau)**2)*np.exp(-(1-tau/(2.+tau))*x)
   return res


nbins_DGGQ=21
nbins=20

#load data
#data k2
################################
k2=2
Q2=3

#data DGGQ
path_data_DGGQ_k2 = path2+'/data/perf/DGGQ/dp/nbins='+str(nbins_DGGQ)+'/kmax='+str(k2)+'/Q='+str(Q2)+'/'
massgridk2_DGGQ = np.loadtxt(path_data_DGGQ_k2+'grainmassgrid.txt')
massbinsk2_DGGQ = np.loadtxt(path_data_DGGQ_k2+'grainmassbins.txt')
massbinsk2meanlog_DGGQ = [np.sqrt(massgridk2_DGGQ[i]*massgridk2_DGGQ[i+1]) for i in range(nbins_DGGQ)]
gjt_xmean_k2_DGGQ = np.loadtxt(path_data_DGGQ_k2+'gjt_xmean.txt')
gjt_xmeanlog_k2_DGGQ = np.loadtxt(path_data_DGGQ_k2+'gjt_xmeanlog.txt')
timek2_DGGQ = np.loadtxt(path_data_DGGQ_k2+'time.txt')
time_perf_DGGQ = np.loadtxt(path_data_DGGQ_k2+'perf.txt',skiprows=4,usecols=4)

#data DG
path_data_DG_k2 = path2+'/data/perf/DG/dp/nbins='+str(nbins)+'/kmax='+str(k2)+'/'
massgridk2 = np.loadtxt(path_data_DG_k2+'grainmassgrid.txt')
massbinsk2 = np.loadtxt(path_data_DG_k2+'grainmassbins.txt')
massbinsk2meanlog = [np.sqrt(massgridk2[i]*massgridk2[i+1]) for i in range(nbins)]
gjt_xmean_k2 = np.loadtxt(path_data_DG_k2+'gjt_xmean.txt')
gjt_xmeanlog_k2 = np.loadtxt(path_data_DG_k2+'gjt_xmeanlog.txt')
timek2 = np.loadtxt(path_data_DG_k2+'time.txt')
time_perf_DG = np.loadtxt(path_data_DG_k2+'perf.txt',skiprows=4,usecols=5)
################################


xmin = np.float(massgridk2[0])
xmax = np.float(massgridk2[-1])

yminloglog = 10**(-20)
ymaxloglog = 10**(-1)

x=np.logspace(np.log10(xmin),np.log10(xmax),num=1000)


plt.figure(1)
plt.ylim(yminloglog,ymaxloglog)
plt.xlim(xmin,xmax)
plt.loglog(x,solKconstDL(x,timek2[-1]),'--',c='C0',label='Analytic')
plt.loglog(massbinsk2meanlog_DGGQ,gjt_xmeanlog_k2_DGGQ[-1,:],markeredgecolor='C1',label=r'$k=2,\,Q=%d,\, Liu\; et.\; al\;  (2019)$' %(Q2),**marker_style)
plt.loglog(massbinsk2meanlog,gjt_xmeanlog_k2[-1,:],markeredgecolor='C2',label=r'$k=2$',**marker_style)
plt.xlabel(r'mass $x$')
plt.ylabel(r'mass density $g(x,\tau)$')
plt.title(r'$\tau=%.2f$' %(timek2[-1]))
plt.legend(loc='lower left',ncol=1)
plt.tight_layout()

print('Improvement in computational time by a factor of %d.' %(time_perf_DGGQ/time_perf_DG))

# plt.savefig(path2+'/plots/kconst_tend_loglog_xmeanlog_DGvsDGGQ.pdf',**savefig_options)
# plt.savefig(path2+'/plots/kconst_tend_loglog_xmeanlog_DGvsDGGQ.png',**savefig_options)


plt.show()
