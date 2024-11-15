# BrainScales
Hierarchies of magnitude and time scales of spontaneous activity fluctuations of whole-brain models

The following functions perform (pseudo-)analytical or numerical estimations of the statistic of different whole-brain models.

Hopf model:

- StochSim_HopfNet.m : stochastic numerical simulations of the network of N hopf nodes.
- StochSim_HopfNet_stimulated.m : stochastic numerical simulations of the network of N hopf nodes under external inputs
- HopfModel_LNA.m :  calculates the covariance, the lagged-covariance, the power spectral density, and the cross-spectrum using the linear approximation

Wilson-Cowan model:

- WCcoupled_StochSim.m : stochastic numerical simulations of the network of N E-I nodes with delayed interactions.
- WCcoupled_LNA_function.m :  calculates the covariance, the lagged-covariance, the power spectral density, and the cross-spectrum using the linear approximation

Linear model:

- LangevinSystem.m : Ornstein–Uhlenbeck process

Connectome:

- DK_connectome.mat : connectivity matrix (DK parcellation)

Regional heterogeneity map:

- Tw1Tw2Map.mat : connectivity matrix for each parcellation.

