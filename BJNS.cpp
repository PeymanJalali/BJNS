#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>


using namespace Rcpp;
using namespace arma;
using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//    function: BJNS
// Purpose: to estimate K graphical models, using the prior function S^2 and an exponential
//          prior on the diagonal elements of the precision matrices
//
// input:
// 1. Y              a list of K p-variate data profiles of possibly different sample sizes
//
// 2. matlist,       list of the matrices used for decomposing the precision matrices (for further info, read section 2 of the paper)
// 
// 3. nmc,           total number of samples generated from the MCMC after sufficient burnins (2000 samples is the default)
//
// 4. burnin,        number burnin sample (2000 samples is the default)
//
// 
// Output:
// 1. an array containing the estimated adjacency matrices of the K precision matrices
// ** one can also export Psis, after the burning stage, in order to perform further Bayesian inference
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// [[Rcpp::export]]
List BJNS( List Y, List matlist, int nmc = 2000, int burnin = 2000){
  
  
  mat tempY = Y(1);
  const int p = tempY.n_cols, K = Y.length();
  const int M = matlist.length();
  
  cube Ss = zeros(p, p, K);
  
  vec n = zeros(K);
  for(int m =0; m<K; m++){
    mat tempY = Y(m);
    n(m) = tempY.n_rows;
    Ss.slice(m) = tempY.t()*tempY;
  }
  // pre-selected values of the hyper parameters  r and s; 
  const double r = 1e-4;
  const double s = 1e-8;
  
  // if one or several of the groups have very small sample sizes (n_k < 0.1*p), we recommend using the following values for r and s
  if(n.min() < 0.1*p){
    const double r = 10;
    const double s = 1;
  } 
  
  
  
  // intilize the matrices present in the full decomposion of the precision matrices
  cube Psis = zeros(p,p, M);
  for(int m=0; m < M; m++){
    mat temp = randn(p,p);
    Psis.slice(m) = (temp + temp.t())/sqrt(2);
    Psis.slice(m).diag() = ones(p);
    if(m < K)
      Psis.slice(m).diag() = ones(p);
    if(m > K-1)
      Psis.slice(m).diag() = zeros(p);
  }
  
  // Votes is an array containging the votes we assign to each vector theta_{jk} as described in the paper
  cube Votes = zeros(p, p, M+1); 
  
  // P contains the model selection probabilities
  vec P = zeros(M+1);
  P(M) = 1;
  
  
  vec ints = zeros(M+1);
  for(int t=0; t<(M+1);t++)
    ints(t) = t;
  
  vec temp_diags = zeros(K);
  vec temp_offdiags = zeros(K);
  
  
  // the following nested for loop is a fully vectorized algorithm that updates the vectors theta_{jk}s and the diagonal elements \theta_{jj}^k 
  // using their full conditional posterior destributions. The details of the Gibbs sampler is written in Section 4 of the paper.
  for(int it = 0; it < (nmc + burnin + 1); it++){
    for(int j=0; j<(p-1); j++){
      for(int k=(j+1); k<p; k++){
        
        mat tempj = zeros(p, K);
        mat tempk = zeros(p,K);
        
        // preliminary quantites as detailed in the Gibbs sampler which are required to calculate the posterio mean and variance of \theta_{jk}s
        for(int mt=0; mt<M; mt++){
          vec templist = matlist(mt);
          for(int m=0; m<K;m++){
            if(any(templist == m)){
              tempj.col(m) += Psis.slice(mt).col(j);
              tempk.col(m) += Psis.slice(mt).col(k);
            }
          }
        }
        
        for(int m=0; m<K; m++){
          temp_diags(m) = Ss(j,j,m) + Ss(k,k,m);
          temp_offdiags(m) = -sum(tempj.col(m).t()*Ss.slice(m).col(k) + tempk.col(m).t()*Ss.slice(m).col(j)) + tempj(k,m)*temp_diags(m);
        }
        
        vec SigmainvMu = zeros(M);
        vec sigma2inv = zeros(M);
        
        for(int mt=0; mt<M; mt++){
          vec templist = matlist(mt);
          for(int m=0; m<K; m++){
            if(any(templist == m)){
              sigma2inv(mt) += temp_diags(m);
              SigmainvMu(mt) += temp_offdiags(m);
            }
          }
          sigma2inv(mt) += R::rgamma(r+0.5, 1/(0.5*Psis(j,k,mt)*Psis(j,k,mt) + s));
          P(mt) = sqrt(2*datum::pi/(sigma2inv(mt)))*exp(0.5*SigmainvMu(mt)*SigmainvMu(mt)/sigma2inv(mt));
        }
        
        P(M) = 1;
        
        P = P/sum(P);
        
        if(P.has_nan()==true){
          double maxdisp = (SigmainvMu%SigmainvMu/sigma2inv).max();
          P(M) = exp(0.5*(-maxdisp));
          for(int mt=0; mt<M; mt++) 
            P(mt) = sqrt(2*datum::pi/(sigma2inv(mt)))*exp(0.5*(SigmainvMu(mt)*SigmainvMu(mt)/sigma2inv(mt) -maxdisp));
          P=P/sum(P);
        }
        
        // now we select the space that theta_{jk} lives in
        int rand = (Rcpp::RcppArmadillo::sample(ints,1,false, P))(0);
        Psis.tube(j,k) = zeros(M);
        Psis.tube(k,j) = zeros(M);
        
        // if rand = M, the theta_{jk} is a zero vector otherwise if rand!= M, by the identifiability constriant, theta_{jk} has only one 
        //non-zero coordiante the index of that coordiante is rand. we update that coordinate using it's full conditonal posterior distribution
        if(rand != (M)){
          Psis(j,k,rand) = Psis(k,j,rand) = R::rnorm(SigmainvMu(rand)/sigma2inv(rand),1/(sigma2inv(rand)));
        }
        
        // we record the vote, given to the vector theta_{jk}
        if(it > burnin)
          Votes(j,k,rand) += 1;
        
        // finally we update the diagonal elements theta_{jj}^k
        
        if(k==p-1){
          for(int m=0; m <K; m++){
            double b = sum(tempj.col(m).t()*Ss.slice(m).col(j)) - Psis(j,j,m)*Ss(j,j,m) + R::rgamma(r+1, 1/(abs(Psis(j,j,m)) + s));
            Psis(j,j,m) = (sqrt(b*b + 4*n(m)*Ss(j,j,m)) - b)/(2*Ss(j,j,m));
          }
        }
        
        
      }
      
    }
    // we have to update the last diagonal elements theta_{pp}^k seperately, because in the above for loop we are respecting the symmetry of the matrices
    // and doing so leaves out the last diagaonl element. Note j goes from 1 to p-1. 
    mat tempj = zeros(p,K);
    for(int mt=0; mt<M; mt++){
      vec templist = matlist(mt);
      for(int m=0; m<K; m++){
        if(any(templist == m)){
          tempj.col(m) += Psis.slice(mt).col(p-1);
        }
      }
    }
    
    for(int m=0; m <K; m++){
      double b = sum(tempj.col(m).t()*Ss.slice(m).col(p-1)) - Psis(p-1,p-1,m)*Ss(p-1,p-1,m) + R::rgamma(r+1, 1/(abs(Psis(p-1,p-1,m)) + s));
      Psis(p-1,p-1,m) = (sqrt(b*b + 4*n(m)*Ss(p-1,p-1,m)) - b)/(2*Ss(p-1,p-1,m));
    }
    
    // output the numebr of iterations once in 1000
    if(it - floor(it/1000)*1000 == 0){
      vec iteration = ones(1)*it;
      iteration.print();
    }
    
  }
  
  
  // reconstruct the adjacency matrices of all Psi^r in the decomposition model
  cube Classes = zeros(p, p, M);
  for(int j=0; j<(p-1); j++)
    for(int k=(j+1); k<p; k++){
      vec tempbs = Votes.tube(j,k);
      int indx = tempbs.index_max();
      if(indx != M){
        Classes(k,j,indx) +=1;
        Classes(j,k,indx) +=1;
      }
    }
    
    // reconstruct the adjacency matrices of all the K precision matrices 
    cube Omegas = zeros(p,p,K);
  for(int m=0; m<K;m++){
    Omegas.slice(m) += eye(p,p);
    for(int mt=0; mt < M; mt++){
      vec templist = matlist(mt);
      if(any(templist == m)){
        Omegas.slice(m) += Classes.slice(mt);
      }
    }
  }
  // output the adjacency matrices of all the matrices in the decomposition model and that of the precison matrices
  return(List::create( Named("Omegas") = Omegas));
}
