# Bayesian-changepoint
Bayesian modeling of a changepoint problem with R implementation

Bayesian Course Project: Bayesian modeling for a changepoint problem
Yaohui Zeng (yaohui-zeng@uiowa.edu), May 13, 2015

Taken from the OpenBUGS example:

  Stagnant: a changepoint problem and an illustration of how NOT
  to do MCMC! http://mathstat.helsinki.fi/openbugs/Examples/Stagnant.html

* M1_Discrete.R: 
 * Model 1: discrete prior for the index of the changepoint, k.
* M2_Continuous.R:
 * Model 2: continuous prior for the changepoint value, x_k.
* Stan_analysis.R
 * RStan implementation

# The implementations reproduce results from OpenBUGS manual.
