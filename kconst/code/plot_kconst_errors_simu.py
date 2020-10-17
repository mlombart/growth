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

#data k0
################################
k0=0
path_data_k0 = path2+'/data/simu/dp/nbins='+str(nbins)+'/kmax='+str(k0)+'/'
timek0 = np.loadtxt(path_data_k0+'time.txt')
errL1_cont_k0 = np.loadtxt(path_data_k0+'errL1_cont.txt')
errL1_dis_k0 = np.loadtxt(path_data_k0+'errL1_dis.txt')
################################

#data k1
################################
k1=1
path_data_k1 = path2+'/data/simu/dp/nbins='+str(nbins)+'/kmax='+str(k1)+'/'
timek1 = np.loadtxt(path_data_k1+'time.txt')
errL1_cont_k1 = np.loadtxt(path_data_k1+'errL1_cont.txt')
errL1_dis_k1 = np.loadtxt(path_data_k1+'errL1_dis.txt')
################################

#data k2
################################
k2=2
path_data_k2 = path2+'/data/simu/dp/nbins='+str(nbins)+'/kmax='+str(k2)+'/'
timek2 = np.loadtxt(path_data_k2+'time.txt')
errL1_cont_k2 = np.loadtxt(path_data_k2+'errL1_cont.txt')
errL1_dis_k2 = np.loadtxt(path_data_k2+'errL1_dis.txt')
################################

#data k3
################################
k3=3
path_data_k3 = path2+'/data/simu/qp/nbins='+str(nbins)+'/kmax='+str(k3)+'/'
timek3 = np.loadtxt(path_data_k3+'time.txt')
errL1_cont_k3 = np.loadtxt(path_data_k3+'errL1_cont.txt')
errL1_dis_k3 = np.loadtxt(path_data_k3+'errL1_dis.txt')
################################


plt.figure(1)
plt.ylim(10**(-5),10)
plt.xlim(timek0[1],timek0[-1])
plt.loglog(timek0[1:],errL1_cont_k0,':',c='black',label=r'$k=0$')
plt.loglog(timek1[1:],errL1_cont_k1,'-.',c='C3',label=r'$k=1$')
plt.loglog(timek2[1:],errL1_cont_k2,'--',c='C1',label=r'$k=2$')  
plt.loglog(timek3[1:],errL1_cont_k3,'-',c='C2',label=r'$k=3$')
plt.xlabel(r'time $\tau$')
plt.ylabel(r'numerical error $e_{\mathrm{c},N}$')
plt.legend(loc='lower left',ncol=2)
plt.tight_layout()

# plt.savefig(path2+'/plots/kconst_errL1cont.pdf',**savefig_options)
# plt.savefig(path2+'/plots/kconst_errL1cont.png',**savefig_options)

# plt.close(1)

plt.figure(2)
plt.ylim(10**(-5),10)
plt.xlim(timek0[1],timek0[-1])
plt.loglog(timek0[1:],errL1_dis_k0,':',c='black',label=r'$k=0$')
plt.loglog(timek1[1:],errL1_dis_k1,'-.',c='C3',label=r'$k=1$')
plt.loglog(timek2[1:],errL1_dis_k2,'--',c='C1',label=r'$k=2$')  
plt.loglog(timek3[1:],errL1_dis_k3,'-',c='C2',label=r'$k=3$')
plt.xlabel(r'time $\tau$')
plt.ylabel(r'numerical error $e_{\mathrm{d},N}$')
plt.legend(loc='lower left',ncol=2)
plt.tight_layout()

# plt.savefig(path2+'/plots/paper_kconst_errL1dis.pdf',**savefig_options)
# plt.savefig(path2+'/plots/kconst_errL1dis.png',**savefig_options)

# plt.close(21)

plt.show()
