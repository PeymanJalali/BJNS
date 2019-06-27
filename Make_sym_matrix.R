# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# function: Make_sym_matrix
# pupose: make a symmetric matrix given a vector of off-digonal elements and a vector of diagonal elements, with appropriate sizes 
#
# inputs:
#
# 1.      a vector of off-digonal elements
#
# 2.      a vector of digonal elements
#
# ouput:
# 
# 1. a symmetrix matrix
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Make_sym_matrix <- function(off_diags, diags){
  p <- length(diags)
  q <- length(off_diags)
  if(q != p*(p-1)/2){
    print("inapproprite sizes for the vectors off_diags and diags")
  } else {
    temp <- matrix(0, p, p)
    temp[lower.tri(temp, diag=FALSE)] <- off_diags
    out <- temp + t(temp)
    out <- out + diag(diags)
    return(out)
  }
}
