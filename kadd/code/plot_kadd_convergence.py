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



listnbins = [5,10,20,40,80]
k0=0
k1=1
k2=2
k3=3

#load data
path_data_k0 = path2+'/data/EOC/dp/'
path_data_k1 = path2+'/data/EOC/dp/'
path_data_k2 = path2+'/data/EOC/qp/'
path_data_k3 = path2+'/data/EOC/qp/'

massgridk0 = np.loadtxt(path_data_k0+'nbins='+str(20)+'/kmax='+str(k0)+'/grainmassgrid.txt')

ndecades = np.log10(massgridk0[-1])-np.log10(massgridk0[0])

Nbins_dec = np.array(listnbins)/ndecades


list_errL1cont_k0=[]
list_errL1cont_k1=[]
list_errL1cont_k2=[]
list_errL1cont_k3=[]

list_errL1dis_k0=[]
list_errL1dis_k1=[]
list_errL1dis_k2=[]
list_errL1dis_k3=[]

for i in listnbins:
   list_errL1cont_k0.append(np.loadtxt(path_data_k0+'nbins='+str(i)+'/kmax='+str(k0)+'/errL1_cont.txt'))
   list_errL1cont_k1.append(np.loadtxt(path_data_k1+'nbins='+str(i)+'/kmax='+str(k1)+'/errL1_cont.txt'))
   list_errL1cont_k2.append(np.loadtxt(path_data_k2+'nbins='+str(i)+'/kmax='+str(k2)+'/errL1_cont.txt'))
   list_errL1cont_k3.append(np.loadtxt(path_data_k3+'nbins='+str(i)+'/kmax='+str(k3)+'/errL1_cont.txt'))

   list_errL1dis_k0.append(np.loadtxt(path_data_k0+'nbins='+str(i)+'/kmax='+str(k0)+'/errL1_dis.txt'))
   list_errL1dis_k1.append(np.loadtxt(path_data_k1+'nbins='+str(i)+'/kmax='+str(k1)+'/errL1_dis.txt'))
   list_errL1dis_k2.append(np.loadtxt(path_data_k2+'nbins='+str(i)+'/kmax='+str(k2)+'/errL1_dis.txt'))
   list_errL1dis_k3.append(np.loadtxt(path_data_k3+'nbins='+str(i)+'/kmax='+str(k3)+'/errL1_dis.txt'))


plt.figure(1)
plt.loglog(Nbins_dec,list_errL1cont_k0,'o',c='black',label=r'$k=0$')  
plt.loglog(Nbins_dec,list_errL1cont_k1,'o',c='C3',label=r'$k=1$')
plt.loglog(Nbins_dec,list_errL1cont_k2,'o',c='C1',label=r'$k=2$')
plt.loglog(Nbins_dec,list_errL1cont_k3,'o',c='C2',label=r'$k=3$')  

plt.loglog(Nbins_dec[1:],1*Nbins_dec[1:]**(-1),':',c='black')
plt.text(6,0.2,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-1}$',color='black')

plt.loglog(Nbins_dec[1:],0.7*Nbins_dec[1:]**(-2),':',c='C3') 
plt.text(6,0.02,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-2}$',color='C3')

plt.loglog(Nbins_dec[1:],0.5*Nbins_dec[1:]**(-3),':',c='C1')
plt.text(6,0.0025,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-3}$',color='C1')

plt.loglog(Nbins_dec[1:],0.1*Nbins_dec[1:]**(-4),':',c='C2') 
plt.text(6,1*10**(-4),r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-4}$',color='C2')

# plt.title(precision+' time=%.2f' %(listtabtimek0[0][-1]))
plt.xlabel(r'$N_{\mathrm{bins}/\mathrm{decade}}$')
plt.ylabel(r'$e_{\mathrm{c},N}$')
plt.xlim(xmax=17)
plt.legend(loc='lower left')

# plt.savefig(path2+'/plots/kadd_errL1cont_convergence.pdf',**savefig_options)
plt.savefig(path2+'/plots/kadd_errL1cont_convergence.png',**savefig_options)

# plt.close(1)


plt.figure(2)
plt.loglog(Nbins_dec,list_errL1dis_k0,'o',c='black',label=r'$k=0$')  
plt.loglog(Nbins_dec,list_errL1dis_k1,'o',c='C3',label=r'$k=1$')
plt.loglog(Nbins_dec,list_errL1dis_k2,'o',c='C1',label=r'$k=2$')
plt.loglog(Nbins_dec,list_errL1dis_k3,'o',c='C2',label=r'$k=3$')  

plt.loglog(Nbins_dec[1:],0.6*Nbins_dec[1:]**(-2),':',c='black')
plt.text(6,0.03,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-2}$',color='black')

plt.loglog(Nbins_dec[1:],0.4*Nbins_dec[1:]**(-2),':',c='C3') 
plt.text(6,0.0015,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-2}$',color='C3')

plt.loglog(Nbins_dec[1:],0.5*Nbins_dec[1:]**(-4),':',c='C1')
plt.text(6,4*10**(-4),r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-4}$',color='C1')

plt.loglog(Nbins_dec[1:],0.1*Nbins_dec[1:]**(-4),':',c='C2') 
plt.text(3,2*10**(-5),r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-4}$',color='C2')

plt.xlabel(r'$N_{\mathrm{bins}/decade}$')
plt.ylabel(r'$e_{\mathrm{d},N}$')
plt.xlim(xmax=17)
plt.legend(loc='lower left')

# plt.savefig(path2+'/plots/kadd_errL1dis_convergence.pdf',**savefig_options)
plt.savefig(path2+'/plots/kadd_errL1dis_convergence.png',**savefig_options)

# plt.close(2)

plt.show()

