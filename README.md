# Quantum-Preparation-Games

This repository is a companion to the article:  

Quantum Preparation Games 
by M. Weilenmann, E.A. Aguilar, and N. Navascués,
IQOQI Vienna - Austrian Academy of Sciences,
3rd November 2020   


The main files are: 


(1) OneShotExample.m   

This contains a loop to create a plot of a e1-e2 tradeoff, with different experimental strategies (Global, LOCC, and Fixed Measurements).  Within the loop the reader can see how each strategy is implemented in practice.      

Requires auxiliary file: paulis.m     

External requirements: QETLAB, CVX, MOSEK  


(2) BoundedConfigurationExample.m     

This file also contains a loop to create a e1-e2 tradeoff plot. There are N measurement rounds, with state_dims specifying the size of the configuration space at round k. The optimization is not guaranteed to be optimal, but in practice the errors behave nicely. Sampling over different random initializations may therefore yield different results. The main optimization function is CS.optimize_omega_state.          

Requires auxiliary files: QI.m, CS.m     

External requirements: YALMIP, MOSEK   


External Requirements::     
QETLAB:     
http://www.qetlab.com/Installation  
In particular: [http://www.qetlab.com/PartialTranspose]   

MOSEK:     
https://www.mosek.com/resources/downloads     
(Academic license needed.)     
https://www.mosek.com/resources/academic-license      

CVX:     
http://cvxr.com/cvx/download/      

YALMIP:     
https://yalmip.github.io/download/
