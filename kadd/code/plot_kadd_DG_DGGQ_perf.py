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


data_solkaddDL_tend = np.genfromtxt(path2+'/data/list_solkadd_tau=3.csv',delimiter=",")


nbins_DGGQ=21
nbins=20


#load data
#data k1
################################
k1=1
Q1=16

#data DGGQ
path_data_DGGQ_k1 = path2+'/data/perf/DGGQ/dp/nbins='+str(nbins_DGGQ)+'/kmax='+str(k1)+'/Q='+str(Q1)+'/'
massgridk1_DGGQ = np.loadtxt(path_data_DGGQ_k1+'grainmassgrid.txt')
massbinsk1_DGGQ = np.loadtxt(path_data_DGGQ_k1+'grainmassbins.txt')
massbinsk1meanlog_DGGQ = [np.sqrt(massgridk1_DGGQ[i]*massgridk1_DGGQ[i+1]) for i in range(nbins_DGGQ)]
gjt_xmean_k1_DGGQ = np.loadtxt(path_data_DGGQ_k1+'gjt_xmean.txt')
gjt_xmeanlog_k1_DGGQ = np.loadtxt(path_data_DGGQ_k1+'gjt_xmeanlog.txt')
timek1_DGGQ = np.loadtxt(path_data_DGGQ_k1+'time.txt')
time_perf_DGGQ = np.loadtxt(path_data_DGGQ_k1+'perf.txt',skiprows=4,usecols=4)

#data DG
path_data_DG_k1 = path2+'/data/perf/DG/dp/nbins='+str(nbins)+'/kmax='+str(k1)+'/'
massgridk1 = np.loadtxt(path_data_DG_k1+'grainmassgrid.txt')
massbinsk1 = np.loadtxt(path_data_DG_k1+'grainmassbins.txt')
massbinsk1meanlog = [np.sqrt(massgridk1[i]*massgridk1[i+1]) for i in range(nbins)]
gjt_xmean_k1 = np.loadtxt(path_data_DG_k1+'gjt_xmean.txt')
gjt_xmeanlog_k1 = np.loadtxt(path_data_DG_k1+'gjt_xmeanlog.txt')
timek1 = np.loadtxt(path_data_DG_k1+'time.txt')
time_perf_DG = np.loadtxt(path_data_DG_k1+'perf.txt',skiprows=4,usecols=5)
################################




xmin = np.float(massgridk1[0])
xmax = np.float(massgridk1[-1])

yminloglog = 10**(-20)
ymaxloglog = 10**(-1)

x=np.logspace(np.log10(xmin),np.log10(xmax),num=1000)


plt.figure(1)
plt.ylim(yminloglog,ymaxloglog)
plt.xlim(xmin,xmax)
plt.loglog(data_solkaddDL_tend[:,0],data_solkaddDL_tend[:,1],'--',c='C0',label='Analytic')
plt.loglog(massbinsk1meanlog_DGGQ,gjt_xmeanlog_k1_DGGQ[-1,:],markeredgecolor='C1',label=r'$k=2,\,Q=%d,\, Liu\; et.\; al\;  (2019)$' %(Q1),**marker_style)
plt.loglog(massbinsk1meanlog,gjt_xmeanlog_k1[-1,:],markeredgecolor='C2',label=r'$k=2$',**marker_style)
plt.xlabel(r'mass $x$')
plt.ylabel(r'mass density $g(x,\tau)$')
plt.title(r'$\tau=%.2f$' %(timek1[-1]))
plt.legend(loc='lower left',ncol=1)
plt.tight_layout()

print('Improvement in computational time by a factor of %2.f.' %(time_perf_DGGQ/time_perf_DG))
# plt.savefig(path2+'/plots/kadd_tend_loglog_xmeanlog_DGvsDGGQ.pdf',**savefig_options)
# plt.savefig(path2+'/plots/kadd_tend_loglog_xmeanlog_DGvsDGGQ.png',**savefig_options)


plt.show()
