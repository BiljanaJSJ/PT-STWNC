Parallel Tempering via Simulated Tempering Without Normalizing Constants

Biljana Jonoska Stojkova and David A. Campbell


In this paper we develop a new general Bayesian methodology that simultaneously estimates 
parameters of interest and the marginal likelihood of the model. 
The proposed methodology builds on Simulated Tempering, which is a powerful algorithm 
that enables sampling from multi-modal distributions.  However, Simulated Tempering 
comes with the practical limitation of needing to specify a prior for the temperature 
along a chosen discretization schedule that will allow calculation of normalizing constants 
at each temperature.  Our proposed model defines the prior for the temperature so as to 
remove the need for calculating normalizing constants at each temperature and thereby enables 
a continuous temperature schedule, while preserving the sampling efficiency of the Simulated 
Tempering algorithm.  The resulting algorithm simultaneously estimates parameters while 
estimating marginal likelihoods through thermodynamic integration.  We illustrate the 
applicability of the new algorithm to different examples involving mixture models of Gaussian 
distributions and ordinary differential equation models.  