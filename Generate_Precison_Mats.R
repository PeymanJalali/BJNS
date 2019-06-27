#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# function: Generate_Precison_Mats
# Porpuse: Generate K precision matrices with random sparse adjacency matrices for all precision matrices. The matrices 
# follow the patter that is given in figure 5 of the paper, where, Graphs in the same row share the same sparsity pattern
# at the bottom right block, whereas graphs in the same column share 
# the same pattern at remaining locations. See page figure 5 on page 22 of the paper.
#
# inputs
# 1. K         # number of the precison matrices to be produced
#
# 2. p         # the size of the precision matrices. All must be equal
#
# 3. rho       # if one intends to add misspecification in which case, some noise will be added to the matrix pattern, explained above
#
# 4. ua and ub # we generate the values of our edges uniformly from uniform(-ub, -ua)UNION uniform(ua, ub)
#
# 5. sparsity  # the sparsity level of the matrices produced
#
#
# Output:
# 1. a list of K positive definite precision matrices
#
#  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gen_Precison_Mats <- function(K, p, rho, ua, ub, sparsity){
  subnetSize = c(p/2, p/2)      # subnet size
  
  subnet_adj_lower_righ_blocks <- vector("list", 2)
  rand.ints <- sample(1:(p/2*(p/2-1)/2), 2*round(sparsity*p/2*(p/2-1)/2), replace = F)
  rand.ints <- split(rand.ints, 1:2)
  temp.unif <- rep(1, p/2*(p/2-1)/2) 
  
  for(i in 1:2){
    temp.ints <- rep(0, p/2*(p/2-1)/2)
    temp.ints[rand.ints[[i]]] <- 1
    temp_subnet = Make_sym_matrix(temp.ints*temp.unif, rep(0,p/2))
    weights <- matrix(0, nrow(temp_subnet), ncol(temp_subnet))
    upperTriangle(weights, diag = F) <- runif((nrow(temp_subnet)*(nrow(temp_subnet) - 1))/2, ua, ub)*(2*rbinom((nrow(temp_subnet)*(nrow(temp_subnet) - 1))/2, 1, 0.5) - 1)
    weights <- weights + t(weights)
    subnet_adj_lower_righ_blocks[[i]] = temp_subnet*weights
  }
  
  n_sub_nets <- round((K + 0.5)/2)
  subnet_adj_lower_not_righ_blocks = vector("list", n_sub_nets)
  
  rand.ints <- sample(1:(p*(p-1)/2), n_sub_nets*round(sparsity*p*(p-1)/2), replace = F)
  rand.ints <- split(rand.ints, 1:n_sub_nets)
  temp.unif <- rep(1, p*(p-1)/2)
  
  
  for(i in 1:n_sub_nets){
    temp.ints <- rep(0, p*(p-1)/2)
    temp.ints[rand.ints[[i]]] <- 1
    temp_subnet = Make_sym_matrix(temp.ints*temp.unif, rep(0,p))
    weights <- matrix(0, nrow(temp_subnet), ncol(temp_subnet))
    upperTriangle(weights, diag = F) <- runif((nrow(temp_subnet)*(nrow(temp_subnet) - 1))/2, ua, ub)*(2*rbinom((nrow(temp_subnet)*(nrow(temp_subnet) - 1))/2, 1, 0.5) - 1)
    weights <- weights + t(weights)
    subnet_adj_lower_not_righ_blocks[[i]] = temp_subnet*weights
  }
  
  
  n_sub_nets <- round((K + 0.5)/2)
  group <- cbind(rep(1:n_sub_nets, rep(2, round((K + 0.5)/2))), rep(c(1,2), round((K + 0.5)/2)))
  
  if(K%%2 != 0){
    group <- group[-nrow(group),]
  }
  
  Omega <- vector("list", K)
  for (k in 1:K){
    Omega[[k]] = subnet_adj_lower_not_righ_blocks[[group[k, 1]]]
    Omega[[k]][(subnetSize[1] + 1):p, (subnetSize[1] + 1):p] = subnet_adj_lower_righ_blocks[[group[k, 2]]]
    if(rho!=0){
      temp_off_diags <- lowerTriangle(Omega[[k]])
      temp_off_diags[sample(which(temp_off_diags==0), round(rho*length(which(temp_off_diags==0))))] <- runif(round(rho*length(which(temp_off_diags==0))), ua, ub)*(2*rbinom(round(rho*length(which(temp_off_diags==0))), 1, 0.5) - 1)
      Omega[[k]] <- Make_sym_matrix(temp_off_diags, rep(0,p))
    }
    ee <- min((eigen(Omega[[k]], only.values = T)$values))
    diag(Omega[[k]]) <- ifelse(ee < 0, -ee + 0.1, 0.1)
  }
  return(Omega)
}
