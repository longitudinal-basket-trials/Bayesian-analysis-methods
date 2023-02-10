# Bayesian-analysis-methods

JAGS code together with R functions to implement methods detailed in 'Bayesian borrowing for basket trials with longitudinal outcomes'.

All code is based upon code from: https://github.com/BasketTrials/Bayesian-analysis-models

Files contained in this repository can be used to reproduce simulation results reported in the paper entitled: Bayesian borrowing for basket trials with longitudinal outcomes (Whitehead et al, 2023)

The files "HD-L.txt", "BHM-L.txt", "EXNEX-L.txt" and "Stratified-L.txt" are for implementing the proposed Bayesian models through JAGS and may be used for analysing basket trial data with continuous longitudinal outcomes.

In particular, "BHM-L.txt" is the (extended 3-level) hierarchical model which assumes the subgroup-specific parameters to be fully exchangeable.

"EXNEX-L.txt" is the borrowing method proposed by Neuenschwander et al. (2016), extended for longitudinal endpoints, where subgroup-specific parameters are partially exchangeable.

"HD-L.txt" is the borrowing method proposed by Zheng and Wason (2022), extended for longitudinal endpoints, which requires no assumption of exchangeability of parameters between subgroups.

"Stratified-L.txt" is the approach of no borrowing or independent Bayesian analysis (a 2-level BHM for longitudinal data) in each subgroup.

Numerical simulation results reported in Section 4 of the paper can be reproduced by calling these models in R with the R2JAGS package. Example R code to implement each of the models via R2JAGS has been provided in the file "Implementation.R".

Contact: l.whitehead2@ncl.ac.uk
