#' Expectation Maximization Poisson Mixtures
#'
#' This function allows you to run an ML - EM algorithm for learning Poisson mixtures
#' @param x this is the data
#' @param K this is the number 
#' @param iter number of iterations defaults to 1000
#' @keywords em
#' @export
#' @examples
#' #generate some Poisson mixture data
#' data <- c(rpois(100,2), rpois(100,10))
#' 
#' #run the ML-EM algorithm to learn the mixture model
#' t1 <- EM_PMM_algo(data,2,iter = 100)
#' 
#' #plot the EM objective function, should be monotonically increasing and monitors convergence
#' plot(t1$sum_of_rows_ln_sum,
#'      type = 'l',
#'      main = "convergence of EM objective function",
#'      bty = 'n')
#' 
#' #check the cluster assignments (using original data this is not a prediction, which can also be done)
#' cluster_indicator <- c()
#' for(i in 1:length(data)){
#'   cluster_indicator <- c(cluster_indicator,which.max(t1$norm_phi[i,]))
#' }
#' plot(data,cluster_indicator,
#'      pch = 16,
#'      bty = 'n')

EM_PMM_algo <- function(x,K,iter = 1000){
  # browser()
  n <- length(x)
  ############
  # Initialize
  ############
  pi_vec <- rep(NA,K) 
  lambda_vec <- rep(NA,K)
  nj_vec <- rep(NA,K)
  phi_iofj <- matrix(NA,nrow = n, ncol = K) #raw phi matrix
  norm_phi <- matrix(NA,nrow = n, ncol = K) #normalized phi matrix
  sum_of_rows <- rep(NA, times = n) #likelihood for one iter
  sum_of_rows_ln <- rep(NA, times = n) #log likelihood for one iter
  sum_of_rows_ln_sum <- rep(NA, times = iter) #store log likelihoods
  
  #initialize pi
  pi_vec <- runif(K)
  pi_vec <- pi_vec/sum(pi_vec)
  
  #initialize lambda
  lambda_vec <- runif(K,0,20)
  
  for(iter in 1:iter){
    #########
    # E Step
    #########
    for(i in 1:n){
      for(j in 1:K){
        phi_iofj[i,j] <- dpois(x[i], lambda_vec[j])*pi_vec[j]
      }
    }
    for(i in 1:n){
      norm_phi[i,] <- phi_iofj[i,]/sum(phi_iofj[i,])
      sum_of_rows[i] <- sum(phi_iofj[i,])
    } 
    
    ################
    # Log Likelihood
    ################
    sum_of_rows_ln <- log(sum_of_rows, base = exp(1))
    sum_of_rows_ln_sum[iter] <- sum(sum_of_rows_ln)
    
    #########
    # M Step
    #########
    nj_vec <- colSums(norm_phi) #update nj
    numerator_lambda <- norm_phi*x
    lambda_vec <- colSums(numerator_lambda)/nj_vec #update lambda
    pi_vec <- nj_vec/n # update pi
  }
  return(list(pi_vec=pi_vec, #probability distribution on Poisson models
              lambda_vec=lambda_vec, #K Poisson parameters
              norm_phi=norm_phi, #cluster assignment distribution
              sum_of_rows_ln_sum = sum_of_rows_ln_sum #objective function
              ))
}
