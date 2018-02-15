# PoissonMixR

Runs Poisson mixture models using variational inference (VI) or expectation maximization (EM) algorithms. 
Users provide count data that they suspect has been generated from a Poisson distribution with a certain parameter “lambda”.
The model uses a maximum likelihood EM approach or a Bayesian variational inference algorithm with Gamma priors on the “lambda”
parameter. The arguments are the data (X), expected number of mixture components (K) and the number of iterations to run the 
algorithm, default is 1000.

