# Grain growth for astrophysics with Discontinuous Galerkin schemes
by Maxime Lombart and Guillaume Laibe

Manuscript is submitted for publication in the Monthly Notices of the Royal Astronomical Society.

This repository contains the data and code used to reproduce all figures shown in the paper.

## Abstract
Depending on their sizes, dust grains store more or less charges, catalyse more or less chemi- cal reactions, intercept more or less photons and stick more or less efficiently to form embryos of planets. Hence the need for an accurate treatment of dust coagulation and fragmentation in numerical modelling. However, existing algorithms for solving the coagulation equation are over-diffusive in the conditions of 3D simulations. We address this challenge by developping a high-order solver based on the Discontinuous Galerkin method. This algorithm conserves mass to machine precision and allows to compute accurately the growth of dust grains over several orders of magnitude in size with a very limited number of dust bins.

## Results
Python 3 is used to produce figures for constant, additive and multiplicative coagulation kernels.

### Benchmarks
#### Constant kernel

##### Positivity, mass conservation and accuracy

<p align="middle">
   <img src="./kconst/plots/kconst_linlog.png" width="80%">
</p>

<div class="row">
   <img src="./kconst/plots/kconst_tend_loglog_xmeanlog.png" width="49%">
   <img src="./kconst/plots/kconst_err_M1.png" width="49%">
</div>

###### Stability in time

<p align="middle">
   <img src="./kconst/plots/kconst_errL1cont.png" width="49%" />
   <img src="./kconst/plots/kconst_errL1dis.png" width="49%" />
</p>

###### Convergence analysis
<p align="middle">
   <img src="./kconst/plots/kconst_errL1_convergence.png" width="49%"/>
   <img src="./kconst/plots/kconst_errL1_xmeanlog_convergence.png" width="49%"/>
</p>


###### Computational efficiency
<div class="row">
   <img src="./kconst/plots/kconst_tend_loglog_xmeanlog_DGvsDGGQ.png" width="400">
   <p>Improvement by a factor ~ 4.</p>
</div>

   
