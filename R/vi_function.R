#' Variational Inference Poisson Mixture
#'
#' This function allows you to run a VI algorithm for learning Poisson mixtures
#' @param x this is the data
#' @param K this is the number 
#' @param iter number of iterations defaults to 1000
#' @param alpha  hyperparameter of probability distribution of Poisson models defaults to 0.1
#' @param a Gamma prior parameter defaults to 0.1
#' @param b Gamma prior parameter defaults to 0.25
#' @keywords vi
#' @export
#' @examples
#' #generate some Poisson mixture data
#' data <- c(rpois(100,2), rpois(100,10))
#' 
#' #run the VI algorithm to learn the mixture model
#' t2 <- VI_PMM_algo(data,2,iter = 1000)
#' 
#' #plot the VI objective function, should be monotonically increasing and monitors convergence
#' plot(t2$L_vi,
#'      type = 'l',
#'      main = "convergence of VI objective function",
#'      bty = 'n')
#'      
#' #check the cluster assignments (using original data this is not a prediction, which can also be done)
#' cluster_indicator <- c()
#' for(i in 1:length(data)){
#'   cluster_indicator <- c(cluster_indicator,which.max(t2$norm_phi[i,]))
#' }
#' plot(data,cluster_indicator,
#'      pch = 16,
#'      bty = 'n')

VI_PMM_algo <- function(x,K,iter = 1000, alpha = 0.1, a = 4.5, b = 0.25){
  # browser()
  n <- length(x)
  ############
  # Initialize
  ############
  # alpha <- 0.1
  # a <- 4.5
  # b <- 0.25
  alpha_k <- rep(NA,K)
  a_k <- rep(NA,K)
  b_k <- rep(NA,K)
  nk_vec <- rep(NA,K)
  phi_iofk <- matrix(NA,nrow = n, ncol = K) #phi matrix
  norm_phi <- matrix(NA,nrow = n, ncol = K) #normalized phi matrix
  phi_xi_mat <- matrix(NA,nrow = n, ncol = K) #phi_xi
  L_vi <- c()
  L_indiv <- 0
  
  #initialize alpha_k
  # alpha_k <- runif(K,1,100) + alpha
  alpha_k <- rgamma(K,1,100) + alpha
  
  #initialize a_k, b_k
  a_k <- rgamma(K,1,100) + a
  b_k <- rgamma(K,1,100) + b
  
  for(iter in 1:iter){
    print(iter)
    #############
    # a. q(c_i) #
    #############
    # browser()
    for(i in 1:n){
      for(j in 1:K){
        phi_iofk[i,j] <- exp((-a_k[j]/b_k[j]) 
                             + (x[i]*(digamma(a_k[j]) - log(b_k[j])))
                             + (digamma(alpha_k[j]) - digamma(sum(alpha_k))))
      }
    }
    for(i in 1:n){
      for(j in 1:K){
        norm_phi[i,j] <- phi_iofk[i,j]/sum(phi_iofk[i,])
      }
    }
    #############
    #e.VI obj fn# 
    #############
    # browser()
    L_indiv <- -sum(((a - a_k)*(digamma(a_k) - log(b_k)))
                    + ((b_k - b)*(a_k/b_k))
                    + ((alpha - alpha_k)*(digamma(alpha_k) - digamma(sum(alpha_k)))))
    L_vi <-  c(L_vi,L_indiv)
    # browser()
    #############
    # b. nk_vec #
    #############
    for(j in 1:K){
      nk_vec[j] <- sum(norm_phi[,j]) #update nk
    }
    
    #############
    # c. q(pi) #
    #############
    # for(j in 1:K){
    #   alpha_k[j] <- alpha + nk_vec[j] #update alpha_k
    # }
    alpha_k <- alpha + nk_vec
    
    #############
    #d.q(lambda)# 
    #############
    for(i in 1:n){
      for(j in 1:K){
        phi_xi_mat[i,j] <- norm_phi[i,j]*x[i]
      }
    }
    for(j in 1:K){
      a_k[j] <- a + sum(phi_xi_mat[,j])
      b_k[j] <- b + nk_vec[j]
    }
  }
  return(list(L_vi = L_vi, #objective function
              a_k = a_k, #parameter a for gamma
              b_k = b_k, #parameter b for gamma
              alpha_k = alpha_k, #parameter for probability of Poisson models
              norm_phi = norm_phi #cluster assignment distribution
              ))
}
