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
![kconst linlog](./kconst/plots/kconst_linlog.png)
![kconst_tend_loglog_xmeanlog](./kconst/plots/kconst_tend_loglog_xmeanlog.png)
![kconst_errL1_cont](./kconst/plots/kconst_errL1cont.png)
<img src="./kconst/plots/kconst_errL1cont.png" width="100">
![kconst_errL1_dis](./kconst/plots/kconst_errL1dis.png)
![kconst_err_M1](./kconst/plots/kconst_err_M1.png)
![kconst_errL1_convergence](./kconst/plots/kconst_errL1_convergence.png)
![kconst_errL1_xmeanlog_convergence](./kconst/plots/kconst_errL1_xmeanlog_convergence.png)
![kconst_tend_loglog_xmeanlog_DGvsDGGQ](./kconst/plots/kconst_tend_loglog_xmeanlog_DGvsDGGQ.png)

