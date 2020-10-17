#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sympy as sym
import numpy as np
from scipy.special import legendre
from matplotlib import pyplot as plt

#change the current directory to the directory
#where the running script file (.py) exists
os.chdir(os.path.dirname(os.path.abspath(__file__)))
path1 = os.getcwd()
# print(path1)
# print(os.path.realpath('..'))
path2 = os.path.realpath('..')


#options for plots
plt.rcParams["font.size"]= 20
plt.rcParams['lines.linewidth'] = 3
savefig_options=dict(bbox_inches='tight')


nbins=20


#load data
#data k0
################################
k0=0
path_data_k0 = path2+'/data/simu/dp/nbins='+str(nbins)+'/kmax='+str(k0)+'/'
timek0 = np.loadtxt(path_data_k0+'time.txt')
errabs_k0 = np.loadtxt(path_data_k0+'/err_abs_moments.txt')
################################

#data k1
################################
k1=1
path_data_k1 = path2+'/data/simu/dp/nbins='+str(nbins)+'/kmax='+str(k1)+'/'
timek1 = np.loadtxt(path_data_k1+'time.txt')
errabs_k1 = np.loadtxt(path_data_k1+'/err_abs_moments.txt')
################################

#data k2
################################
k2=2
path_data_k2 = path2+'/data/simu/qp/nbins='+str(nbins)+'/kmax='+str(k2)+'/'
timek2 = np.loadtxt(path_data_k2+'time.txt')
errabs_k2 = np.loadtxt(path_data_k2+'/err_abs_moments.txt')
################################

#data k3
################################
k3=3
path_data_k3 = path2+'/data/simu/qp/nbins='+str(nbins)+'/kmax='+str(k3)+'/'
timek3 = np.loadtxt(path_data_k3+'time.txt')
errabs_k3 = np.loadtxt(path_data_k3+'/err_abs_moments.txt')
################################



plt.figure(1)
plt.ylim(10**(-7),10**(-4))
plt.xlim(0.01,timek0[-1])
plt.loglog(timek0[1:],errabs_k0[:,1],':',c='black',label=r'$k=0$')
plt.loglog(timek1[1:],errabs_k1[:,1],'-.',c='C3',label=r'$k=1$')
plt.loglog(timek2[1:],errabs_k2[:,1],'--',c='C1',label=r'$k=2$')  
plt.loglog(timek3[1:],errabs_k3[:,1],'-',c='C2',label=r'$k=3$')
plt.xlabel(r'time $\tau$')
plt.ylabel(r'numerical error $e_{M_1,N}$')
plt.legend(loc='upper left',ncol=2)
plt.tight_layout()


# plt.savefig(path2+'/plots/kadd_err_M1.pdf',**savefig_options)
# plt.savefig(path2+'/plots/kadd_err_M1.png',**savefig_options)

plt.show()