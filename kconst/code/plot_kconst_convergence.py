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
path_data = path2+'/data/EOC/dp/'
path_data_k3 = path2+'/data/EOC/qp/'

massgridk0 = np.loadtxt(path_data+'nbins='+str(5)+'/kmax='+str(k0)+'/grainmassgrid.txt')
massgridk1 = np.loadtxt(path_data+'nbins='+str(5)+'/kmax='+str(k1)+'/grainmassgrid.txt')
massgridk2 = np.loadtxt(path_data+'nbins='+str(5)+'/kmax='+str(k2)+'/grainmassgrid.txt')
massgridk3 = np.loadtxt(path_data_k3+'nbins='+str(5)+'/kmax='+str(k3)+'/grainmassgrid.txt')

#define array number of bins per decade
ndecades = np.log10(massgridk0[-1])-np.log10(massgridk0[0])

Nbins_dec = np.array(listnbins)/ndecades

listtabtimek0=[]
listtabtimek1=[]
listtabtimek2=[]
listtabtimek3=[]

listtaberrL1k0=[]
listtaberrL1k1=[]
listtaberrL1k2=[]
listtaberrL1k3=[]

listtaberrL1disk0=[]
listtaberrL1disk1=[]
listtaberrL1disk2=[]
listtaberrL1disk3=[]

listtaberrL1_xmeanlog_k0=[]
listtaberrL1_xmeanlog_k1=[]
listtaberrL1_xmeanlog_k2=[]
listtaberrL1_xmeanlog_k3=[]

for i in listnbins:
   listtaberrL1k0.append(np.loadtxt(path_data+'nbins='+str(i)+'/kmax='+str(k0)+'/errL1.txt'))
   listtaberrL1k1.append(np.loadtxt(path_data+'nbins='+str(i)+'/kmax='+str(k1)+'/errL1.txt'))
   listtaberrL1k2.append(np.loadtxt(path_data+'nbins='+str(i)+'/kmax='+str(k2)+'/errL1.txt'))
   listtaberrL1k3.append(np.loadtxt(path_data_k3+'nbins='+str(i)+'/kmax='+str(k3)+'/errL1.txt'))

   listtaberrL1disk0.append(np.loadtxt(path_data+'nbins='+str(i)+'/kmax='+str(k0)+'/errL1_massbins.txt'))
   listtaberrL1disk1.append(np.loadtxt(path_data+'nbins='+str(i)+'/kmax='+str(k1)+'/errL1_massbins.txt'))
   listtaberrL1disk2.append(np.loadtxt(path_data+'nbins='+str(i)+'/kmax='+str(k2)+'/errL1_massbins.txt'))
   listtaberrL1disk3.append(np.loadtxt(path_data_k3+'nbins='+str(i)+'/kmax='+str(k3)+'/errL1_massbins.txt'))

   listtaberrL1_xmeanlog_k0.append(np.loadtxt(path_data+'nbins='+str(i)+'/kmax='+str(k0)+'/errL1_xmeanlog.txt'))
   listtaberrL1_xmeanlog_k1.append(np.loadtxt(path_data+'nbins='+str(i)+'/kmax='+str(k1)+'/errL1_xmeanlog.txt'))
   listtaberrL1_xmeanlog_k2.append(np.loadtxt(path_data+'nbins='+str(i)+'/kmax='+str(k2)+'/errL1_xmeanlog.txt'))
   listtaberrL1_xmeanlog_k3.append(np.loadtxt(path_data_k3+'nbins='+str(i)+'/kmax='+str(k3)+'/errL1_xmeanlog.txt'))

   listtabtimek0.append(np.loadtxt(path_data+'nbins='+str(i)+'/kmax='+str(k0)+'/time.txt'))
   listtabtimek1.append(np.loadtxt(path_data+'nbins='+str(i)+'/kmax='+str(k1)+'/time.txt'))
   listtabtimek2.append(np.loadtxt(path_data+'nbins='+str(i)+'/kmax='+str(k2)+'/time.txt'))
   listtabtimek3.append(np.loadtxt(path_data_k3+'nbins='+str(i)+'/kmax='+str(k3)+'/time.txt'))

# check same time for each order k
# print('same time k0 ? ->',listtabtimek0[0][-1]==listtabtimek0[1][-1]==listtabtimek0[2][-1]==listtabtimek0[3][-1])
# print('same time k1 ? ->',listtabtimek1[0][-1]==listtabtimek1[1][-1]==listtabtimek1[2][-1]==listtabtimek1[3][-1])
# print('same time k2 ? ->',listtabtimek2[0][-1]==listtabtimek2[1][-1]==listtabtimek2[2][-1]==listtabtimek2[3][-1])
# print('same time k3 ? ->',listtabtimek3[0][-1]==listtabtimek3[1][-1]==listtabtimek3[2][-1]==listtabtimek3[3][-1])
# print('same time all k ? ->',listtabtimek0[0][-1]==listtabtimek1[0][-1]==listtabtimek2[0][-1]==listtabtimek3[0][-1])

#select errors at tend=0.01
#error L1 continuous
errL1k0 = [listtaberrL1k0[j][-1] for j in range(len(listnbins))]
errL1k1 = [listtaberrL1k1[j][-1] for j in range(len(listnbins))]
errL1k2 = [listtaberrL1k2[j][-1] for j in range(len(listnbins))]
errL1k3 = [listtaberrL1k3[j][-1] for j in range(len(listnbins))]

#error L1 discrte at xmean
errL1disk0 = [listtaberrL1disk0[j][-1] for j in range(len(listnbins))]
errL1disk1 = [listtaberrL1disk1[j][-1] for j in range(len(listnbins))]
errL1disk2 = [listtaberrL1disk2[j][-1] for j in range(len(listnbins))]
errL1disk3 = [listtaberrL1disk3[j][-1] for j in range(len(listnbins))]

#error L1 discrete at xmeanlog
errL1_xmeanlog_k0 = [listtaberrL1_xmeanlog_k0[j][-1] for j in range(len(listnbins))]
errL1_xmeanlog_k1 = [listtaberrL1_xmeanlog_k1[j][-1] for j in range(len(listnbins))]
errL1_xmeanlog_k2 = [listtaberrL1_xmeanlog_k2[j][-1] for j in range(len(listnbins))]
errL1_xmeanlog_k3 = [listtaberrL1_xmeanlog_k3[j][-1] for j in range(len(listnbins))]




plt.figure(1)
plt.loglog(Nbins_dec,errL1k0,'o',c='black',label=r'$k=0$')  
plt.loglog(Nbins_dec,errL1k1,'o',c='C3',label=r'$k=1$')
plt.loglog(Nbins_dec,errL1k2,'o',c='C1',label=r'$k=2$')
plt.loglog(Nbins_dec,errL1k3,'o',c='C2',label=r'$k=3$')  

plt.loglog(Nbins_dec[1:],1*Nbins_dec[1:]**(-1),':',c='black')
plt.text(6,0.2,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-1}$',color='black')

plt.loglog(Nbins_dec[1:],0.7*Nbins_dec[1:]**(-2),':',c='C3') 
plt.text(6,0.02,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-2}$',color='C3')

plt.loglog(Nbins_dec[1:],0.5*Nbins_dec[1:]**(-3),':',c='C1')
plt.text(6,0.0025,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-3}$',color='C1')

plt.loglog(Nbins_dec[1:],0.1*Nbins_dec[1:]**(-4),':',c='C2') 
plt.text(6,1*10**(-4),r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-4}$',color='C2')

plt.xlabel(r'$N_{\mathrm{bins}/\mathrm{decade}}$')
plt.ylabel(r'$e_{\mathrm{c},N}$')
plt.xlim(xmax=17)
plt.legend(loc='lower left')
plt.tight_layout()

# plt.savefig(path2+'/plots/kconst_errL1_convergence.pdf',**savefig_options)
# plt.savefig(path2+'/plots/kconst_errL1_convergence.png',**savefig_options)

# plt.close(1)


# plt.figure(2)
# plt.loglog(Nbins_dec,errL1disk0,'o',c='black',label=r'$k=0$')  
# plt.loglog(Nbins_dec,errL1disk1,'o',c='C3',label=r'$k=1$')
# plt.loglog(Nbins_dec,errL1disk2,'o',c='C1',label=r'$k=2$')
# plt.loglog(Nbins_dec,errL1disk3,'o',c='C2',label=r'$k=3$')   
# plt.xlabel(r'$N_{\mathrm{bins}}/decade$')
# plt.ylabel(r'$e_{\mathrm{c},N}$')
# plt.xlim(xmax=17)
# plt.legend(loc='lower left')
# plt.close(2)

plt.figure(3)
plt.loglog(Nbins_dec,errL1_xmeanlog_k0,'o',c='black',label=r'$k=0$')  
plt.loglog(Nbins_dec,errL1_xmeanlog_k1,'o',c='C3',label=r'$k=1$')
plt.loglog(Nbins_dec,errL1_xmeanlog_k2,'o',c='C1',label=r'$k=2$')
plt.loglog(Nbins_dec,errL1_xmeanlog_k3,'o',c='C2',label=r'$k=3$')   

plt.loglog(Nbins_dec[1:],1*Nbins_dec[1:]**(-2),':',c='black')
plt.text(6,0.05,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-2}$',color='black')

plt.loglog(Nbins_dec[1:],0.3*Nbins_dec[1:]**(-2),':',c='C3') 
plt.text(6,0.0015,r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-2}$',color='C3')

plt.loglog(Nbins_dec[1:],0.5*Nbins_dec[1:]**(-4),':',c='C1')
plt.text(6,4*10**(-4),r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-4}$',color='C1')

plt.loglog(Nbins_dec[1:],0.1*Nbins_dec[1:]**(-4),':',c='C2') 
plt.text(3,2*10**(-5),r'$\propto N_{\mathrm{bins}/\mathrm{dec}}^{-4}$',color='C2')

plt.xlabel(r'$N_{\mathrm{bins}/decade}$')
plt.ylabel(r'$e_{\mathrm{d},N}$')
plt.xlim(xmax=17)
plt.legend(loc='lower left')
plt.tight_layout()

# plt.savefig(path2+'/plots/kconst_errL1_xmeanlog_convergence.pdf',**savefig_options)
# plt.savefig(path2+'/plots/kconst_errL1_xmeanlog_convergence.png',**savefig_options)

# plt.close(3)

plt.show()
