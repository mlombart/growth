#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sympy as sym
import numpy as np
from scipy.special import legendre
from matplotlib import pyplot as plt

#options for plots
plt.rcParams["font.size"]= 16
plt.rcParams['lines.linewidth'] = 3

savefig_options=dict(bbox_inches='tight')


#legendre polynomials
def LegendreP(i,x):
   coeffs=legendre(i)
   res = 0
   for j in range(0,i+1):
      res = res+ coeffs[j]*x**j

   return res

# y=sym.Symbol('y')
# print(y**0)
# print('legendre(i)=',np.shape(legendre(3)))
# print(legendre(2)[0])
# print(range(0,3))
# print(LegendreP(2,y))
# print("leg=",LegendreP(1,np.array([0.5,0.5])))

#ghtilde
def ghtilde(massgrid,massbins,gij,gh_mean,theta,k,j,x):
   xij = 2/(massgrid[j+1]-massgrid[j])*(x-massbins[j])
   res1 = 0 
   if k==0:
      res1 = theta[j]*(gij[j]*LegendreP(k,xij)-gh_mean[j])+gh_mean[j]
   else:
      for i in range(k+1):
         res1 = res1+ theta[j]*(gij[j,i]*LegendreP(i,xij)-gh_mean[j])+gh_mean[j]  
   return res1


def I(massgrid,j):
   res= np.logspace(np.log10(massgrid[j]),np.log10(massgrid[j+1]),num=100)
   return res




#solution kconst g(x,0)=x exp(-x)
def solKconstDL(x,tau):
   res = 4.*x/((2+tau)**2)*np.exp(-(1-tau/(2.+tau))*x)
   return res


marker_style = dict( marker='o',markersize=8, markerfacecolor='white', linestyle='',markeredgewidth=1.2)



nbins=20
wp=8

#load data fortran
if wp==8:
   path_data='/Users/maxime/Travail/These_ENS_lyon/travail/growth_code/tests_code/check_algo_fortran/kconst_fortran/simu/dp/'
   precision = 'dp'
elif wp==16:
   path_data='/Users/maxime/Travail/These_ENS_lyon/travail/growth_code/tests_code/check_algo_fortran/kconst_fortran/simu/qp/'
   precision = 'qp'

path_data_k3 = '/Users/maxime/Travail/These_ENS_lyon/travail/growth_code/tests_code/check_algo_fortran/kconst_fortran/simu/qp/'

path_plots ='/Users/maxime/Travail/These_ENS_lyon/travail/growth_code/tests_code/check_algo_fortran/plots_kconst/DG/simu/'

path_plots_these = '/Users/maxime/Dropbox/laibe_lombart/redaction_these/figures/'

path_plots_paper = '/Users/maxime/Dropbox/laibe_lombart/Papers/Paper_1/figures/'


#data k0
################################
k0=0
massgridk0 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k0)+'/grainmassgrid.txt')
massbinsk0 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k0)+'/grainmassbins.txt')

massbinsk0meanlog = [np.sqrt(massgridk0[i]*massgridk0[i+1]) for i in range(nbins)]

gij_k0_data = np.genfromtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k0)+'/gijLeg.txt')
ghmean_k0 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k0)+'/ghmean.txt')
theta_k0 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k0)+'/theta.txt')

timek0 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k0)+'/time.txt')

gij_k0 = np.reshape(gij_k0_data,(len(timek0),nbins,k0+1))
################################


#data k1
################################
k1=1
massgridk1 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k1)+'/grainmassgrid.txt')
massbinsk1 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k1)+'/grainmassbins.txt')

massbinsk1meanlog = [np.sqrt(massgridk1[i]*massgridk1[i+1]) for i in range(nbins)]

gij_k1_data = np.genfromtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k1)+'/gijLeg.txt')
ghmean_k1 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k1)+'/ghmean.txt')
theta_k1 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k1)+'/theta.txt')

timek1 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k1)+'/time.txt')

gij_k1 = np.reshape(gij_k1_data,(len(timek1),nbins,k1+1))
################################

#data k2
################################
k2=2
massgridk2 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k2)+'/grainmassgrid.txt')
massbinsk2 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k2)+'/grainmassbins.txt')

massbinsk2meanlog = [np.sqrt(massgridk2[i]*massgridk2[i+1]) for i in range(nbins)]

gij_k2_data = np.genfromtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k2)+'/gijLeg.txt')
ghmean_k2 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k2)+'/ghmean.txt')
theta_k2 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k2)+'/theta.txt')

timek2 = np.loadtxt(path_data+'nbins='+str(nbins)+'/kmax='+str(k2)+'/time.txt')

gij_k2 = np.reshape(gij_k2_data,(len(timek2),nbins,k2+1))

################################

#data k3 t0
################################
k3=3
massgridk3 = np.loadtxt(path_data_k3+'nbins='+str(nbins)+'/kmax='+str(k3)+'/grainmassgrid.txt')
massbinsk3 = np.loadtxt(path_data_k3+'nbins='+str(nbins)+'/kmax='+str(k3)+'/grainmassbins.txt')

massbinsk3meanlog = [np.sqrt(massgridk3[i]*massgridk3[i+1]) for i in range(nbins)]

gij_k3_data = np.genfromtxt(path_data_k3+'nbins='+str(nbins)+'/kmax='+str(k3)+'/gijLeg.txt')
ghmean_k3 = np.loadtxt(path_data_k3+'nbins='+str(nbins)+'/kmax='+str(k3)+'/ghmean.txt')
theta_k3 = np.loadtxt(path_data_k3+'nbins='+str(nbins)+'/kmax='+str(k3)+'/theta.txt')

timek3 = np.loadtxt(path_data_k3+'nbins='+str(nbins)+'/kmax='+str(k3)+'/time.txt')

gij_k3 = np.reshape(gij_k3_data,(len(timek3),nbins,k3+1))

################################

xmin = np.float(massgridk2[0])
xmax = np.float(massgridk2[-1])

print(xmin)
print(xmax)

print(np.log10(xmin))

ymint0linlog=-0.01
ymaxt0linlog=0.4

ymintendlinlog=-10**(-5)
ymaxtendlinlog=10** (-5)

yminloglog = 10**(-16)
ymaxloglog = 1

x=np.logspace(np.log10(xmin),np.log10(xmax),num=1000)



fig, axes = plt.subplots(4,2,figsize=(10,12),sharex='col', gridspec_kw={'hspace': 0, 'wspace': 0.05})


for (m,n), subplot in np.ndenumerate(axes):
   axes[m,n].set_xscale('log')
   axes[m,n].ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
   axes[m,n].set_xlim(xmin,xmax)
   axes[m,n].autoscale(enable=True, axis="y", tight=False)
   axes[m,n].axvline(massgridk0[0],c='grey',alpha=0.3)
   axes[m,n].axhline(0,c='grey',alpha=0.3,linestyle='--')
   axes[m,n].legend()
   axes[m,1].yaxis.tick_right()
   for j in range(nbins):
      axes[m,n].axvline(massgridk2[j+1],ymin=-0.01,c='grey',alpha=0.3)

axes[0,0].plot(x,solKconstDL(x,timek0[0]),'--',c='C0',label='Analytic')
axes[0,1].plot(x,solKconstDL(x,timek0[-1]),'--',c='C0',label='Analytic')
axes[1,0].plot(x,solKconstDL(x,timek1[0]),'--',c='C0',label='Analytic')
axes[1,1].plot(x,solKconstDL(x,timek1[-1]),'--',c='C0',label='Analytic')
axes[2,0].plot(x,solKconstDL(x,timek2[0]),'--',c='C0',label='Analytic')
axes[2,1].plot(x,solKconstDL(x,timek2[-1]),'--',c='C0',label='Analytic')
axes[3,0].plot(x,solKconstDL(x,timek3[0]),'--',c='C0',label='Analytic')
axes[3,1].plot(x,solKconstDL(x,timek3[-1]),'--',c='C0',label='Analytic')



for j in range(nbins):
   axes[0,0].plot(I(massgridk0,j),ghtilde(massgridk0,massbinsk0,gij_k0[0,:,:],ghmean_k0[0,:],theta_k0[0,:],k0,j,I(massgridk0,j)),c='black',label=(r'$g_j(x,\tau)$' if j==0 else '_'))
   axes[0,1].plot(I(massgridk0,j),ghtilde(massgridk0,massbinsk0,gij_k0[-1,:,:],ghmean_k0[-1,:],theta_k0[-1,:],k0,j,I(massgridk0,j)),c='black',label=(r'$g_j(x,\tau)$' if j==0 else '_'))
   axes[1,0].plot(I(massgridk1,j),ghtilde(massgridk0,massbinsk0,gij_k1[0,:,:],ghmean_k1[0,:],theta_k1[0,:],k1,j,I(massgridk0,j)),c='C3',label=(r'$g_j(x,\tau)$' if j==0 else '_'))
   axes[1,1].plot(I(massgridk1,j),ghtilde(massgridk0,massbinsk0,gij_k1[-1,:,:],ghmean_k1[-1,:],theta_k1[-1,:],k1,j,I(massgridk1,j)),c='C3',label=(r'$g_j(x,\tau)$' if j==0 else '_'))
   axes[2,0].plot(I(massgridk2,j),ghtilde(massgridk2,massbinsk2,gij_k2[0,:,:],ghmean_k2[0,:],theta_k2[0,:],k2,j,I(massgridk2,j)),c='C1',label=(r'$g_j(x,\tau)$' if j==0 else '_'))
   axes[2,1].plot(I(massgridk2,j),ghtilde(massgridk2,massbinsk2,gij_k2[-1,:,:],ghmean_k2[-1,:],theta_k2[-1,:],k2,j,I(massgridk2,j)),c='C1',label=(r'$g_j(x,\tau)$' if j==0 else '_'))
   axes[3,0].plot(I(massgridk3,j),ghtilde(massgridk3,massbinsk3,gij_k3[0,:,:],ghmean_k3[0,:],theta_k3[0,:],k3,j,I(massgridk3,j)),c='C2',label=(r'$g_j(x,\tau)$' if j==0 else '_'))
   axes[3,1].plot(I(massgridk3,j),ghtilde(massgridk3,massbinsk3,gij_k3[-1,:,:],ghmean_k3[-1,:],theta_k3[-1,:],k3,j,I(massgridk3,j)),c='C2',label=(r'$g_j(x,\tau)$' if j==0 else '_'))
   

axes[0,0].plot([], [], ' ', label=r'$k=0$')
axes[1,0].plot([], [], ' ', label=r'$k=1$')
axes[2,0].plot([], [], ' ', label=r'$k=2$')
axes[3,0].plot([], [], ' ', label=r'$k=3$')

axes[0,0].legend()
axes[1,0].legend()
axes[2,0].legend()
axes[3,0].legend()

axes[0,0].set_title(r'$\tau=%d$' %(timek3[0]))
axes[0,1].set_title(r'$\tau=%d$' %(timek3[-1]))

axes[1,0].yaxis.get_offset_text().set_visible(False)
axes[1,1].yaxis.get_offset_text().set_visible(False)
axes[2,0].yaxis.get_offset_text().set_visible(False)
axes[2,1].yaxis.get_offset_text().set_visible(False)
axes[3,0].yaxis.get_offset_text().set_visible(False)
axes[3,1].yaxis.get_offset_text().set_visible(False)

for j in range(4):
   axes[j,0].set_ylabel(r'mass density $g$')
axes[3,0].set_xlabel(r'mass $x$')
axes[3,1].set_xlabel(r'mass $x$')

plt.savefig(path_plots+'paper_kconst_linlog.pdf',**savefig_options)
plt.savefig(path_plots_these+'paper_kconst_linlog.pdf',**savefig_options)
plt.savefig(path_plots_paper+'paper_kconst_linlog.pdf',**savefig_options)


plt.show()

