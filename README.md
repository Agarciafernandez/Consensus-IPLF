# Consensus-IPLF

This repository contains the implementation of the consensus iterated posterior linearisation filter (IPLF) for distributed state estimation using a sensor network. The algorithm is described in [1]. The (centralised) IPLF is described in [2] for additive Gaussian models and in [3] for dynamic and measurement models described by their conditional moments.

The code enables the choice of several consensus strategies [4]: consensus on information (CI), consensus on measurement (CM) and hybrid consensus on measurements and consensus on information (HCMCI).

Consensus_IPLF.m runs the consensus IPLF.
Consensus_IPLF.m runs the centralised IPLF.

File circ_vmrnd.m comes from the circular statistics toolbox for Matlab [5]. The consensus for linear models and graph description parts are based on a code provided by Prof. Giorgio Battistelli and written by Dr. Xiangxiang Dong.




[1] A. F. García-Fernández, G. Batisttelli, "Consensus iterated posterior linearisation filter for distributed state estimation", IEEE Signal Processing Letters, 2025.

[2] Á. F. García-Fernández, L. Svensson, M. R. Morelande and S. Särkkä, "Posterior Linearization Filter: Principles and Implementation Using Sigma Points," in IEEE Transactions on Signal Processing, vol. 63, no. 20, pp. 5561-5573, Oct.15, 2015, doi: 10.1109/TSP.2015.2454485.

[3] F. Tronarp, Á. F. García-Fernández and S. Särkkä, "Iterative Filtering and Smoothing in Nonlinear and Non-Gaussian Systems Using Conditional Moments," in IEEE Signal Processing Letters, vol. 25, no. 3, pp. 408-412, March 2018, doi: 10.1109/LSP.2018.2794767.

[4] G. Battistelli, L. Chisci, G. Mugnai, A. Farina and A. Graziano, "Consensus-Based Linear and Nonlinear Filtering," in IEEE Transactions on Automatic Control, vol. 60, no. 5, pp. 1410-1415, May 2015.

[5] P. Berens, “CircStat: A MATLAB toolbox for circular statistics,” Journal of Statistical Software, vol. 31, pp. 1–21, Sep. 2009


