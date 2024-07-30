h2_gamma1_TBD <- c(0, 10^-2, 10^-1)

for( k_gamma1 in 1:length(h2_gamma1_TBD) ){
  # k_gamma1 <- 1
  h2_gamma1 <- h2_gamma1_TBD[k_gamma1]
  h2_delta  <- h2_gamma1*4/5
  h2_kappax <- h2_gamma1/5
  
  comparison <- function( seed, h2_delta, h2_kappax ){
    library(MASS)
    library(MendelianRandomization)
    library(MCMCpack)
    library(dlm)
    library(distr)
    library(progress)
    library(dplyr)
    library(PMVMR)
    
    set.seed(seed)
    genotype_continuous <- function( sample_size, nsnp, r_LD ){
      miu <- rep(0, nsnp)   ### Set mean variable
      Sigma <- matrix( 1:nsnp^2, ncol = nsnp )   #### Set initial covariance matrix
      for( i in 1:nsnp ){
        Sigma[i, ] <- r_LD^abs( i - (1:nsnp) )   ### Calculate the covariance matrix
      }  ### Elements for the i-th row are r^| i - n |
      return( mvrnorm( sample_size, miu, Sigma ) )
    }
    
    genotype_to_factor  <- function( gene_matrix, minor_allel_frequency = c() ){
      MAF_cut <- function( x ){
        MAF <- runif( 1, minor_allel_frequency[1], minor_allel_frequency[2] )
        qt_aa <- MAF^2
        qt_Aa <- MAF*(2 - MAF)   ### 0___aa(MAF*MAF)___Aa(MAF*(2-MAF))___1
        as.integer( cut(x, breaks = c(-Inf, quantile( x, c(qt_aa, qt_Aa) ), Inf), labels = c(0, 1, 2) ) ) - 1
      }
      G <- apply(gene_matrix, 2, MAF_cut)
      return(G)
    }
    
    calculateP <- function( mean, sd ){
      LCI = mean - 1.96*sd; UCI = mean + 1.96*sd;
      if( LCI >= 0 ){
        pvalue = pnorm( 0, mean, sd, T )
      }else if( UCI <= 0 ){
        pvalue = pnorm( 0, mean, sd, F )
      }else{
        pLower = pnorm( 0, mean, sd, T)
        pUpper = pnorm( 0, mean, sd, F)
        pvalue = 2*min( pLower, pUpper )
      }
      return( pvalue )
    }
    
    beta1_true <- ifelse( h2_delta == 0, 0, 0.01 )   ### beta1_true = 0, when and only when h2_delta = h2_kappax = 0
    
    data_generating <- function( 
      nc = 1000, nr = 20, nconf = 20,  # nc for common SNPs size; nr for rare SNPs; nr for confounders
      ns = 50000, ni = 50000,          # ns for sumamry-level data size; ni for individual-level
      r_kappa = 0.2, r_LD = 0.2,
      h2_delta = h2_delta, h2_kappax = h2_kappax, h2_gamma2 = 0.02, h2_kappay = 0.02, h2_theta = 0.02,
      beta1 = beta1_true, beta2 = 0.2 )
    {
      # nc = 1000; nr = 20; nconf = 20;  # nc for common SNPs size; nr for rare SNPs; nr for confounders
      # ns = 50000; ni = 50000;          # ns for sumamry-level data size; ni for individual-level
      # r_kappa = 0.5; r_LD = 0.2;
      # h2_delta = 0.02; h2_kappax = 0.02; h2_gamma2 = 0.02; h2_kappay = 0.05; h2_theta = 0.02;
      # beta1 = beta1_true; beta2 = 0.2;
      ########## Initialize parameters ##########
      nsnp <- nc + nr;  nsample <- ns + ni
      delta <- rnorm( nsnp, 0, 1 ); 
      Sigma_kappa <- matrix( c( 1, r_kappa, r_kappa, 1 ), nrow = 2 )
      kappa  <- mvrnorm( nsnp, mu = c(0, 0), Sigma = Sigma_kappa )
      kappax <- kappa[ ,1]; kappay <- kappa[ ,2]
      gamma2 <- rnorm( nsnp, 0, 1)
      theta <- rnorm( nsnp, 0, 1 )
      fai <- mvrnorm( nconf, mu = c(0, 0), Sigma = matrix( c(1, 0, 0, 1 ), nrow = 2 ) )
      faix <- fai[ ,1]; faiy <- fai[ ,2]
      
      ########## Generate genotype and confounders ##########
      Gc <- genotype_continuous( sample_size = nsample, nsnp = nc, r_LD = r_LD ) %>% 
        genotype_to_factor(., c(0.05, 0.5))
      Gr <- genotype_continuous( sample_size = nsample, nsnp = nr, r_LD = 0.5 ) %>% 
        genotype_to_factor(., c(0.005, 0.01))
      G  <- cbind(Gc, Gr)
      U  <- matrix( rnorm( nsample*nconf ), ncol = nconf )
      
      ########## Generate and scale data using structural model ##########
      k_scale <- ( beta1^2 + beta2^2 + 1 )/( 1 - ( h2_delta + h2_kappax + h2_gamma2 + h2_kappay + h2_theta ) )
      G_delta <- G %*% delta; G_kappax <- G %*% kappax;
      if( beta1 != 0 ){
        delta  <- delta *as.numeric( sqrt( h2_delta *k_scale/( beta1^2*var( G_delta  ) ) ) )
        kappax <- kappax*as.numeric( sqrt( h2_kappax*k_scale/( beta1^2*var( G_kappax ) ) ) )
        G_delta <- G %*% delta; G_kappax <- G %*% kappax
      }
      
      G_gamma2 <- G %*% gamma2
      if( beta2 != 0 ){
        gamma2 <- gamma2*as.numeric( sqrt( h2_gamma2*k_scale/( beta2^2*var( G_gamma2 ) ) ) )
        G_gamma2 <- G %*% gamma2
      }
      
      if( h2_kappay == 0 ){
        G_kappay <- G %*% rep( 0, nsnp )
      } else{
        G_kappay <- G %*% kappay;
        kappay <- kappay*as.numeric( sqrt( h2_kappay*k_scale/  var( G_kappay ) ) )
        G_kappay <- G %*% kappay
      }
      
      if( h2_theta == 0 ){
        G_theta  <- G %*% rep( 0, nsnp )
      } else{ 
        G_theta  <- G %*% theta
        theta  <- theta *as.numeric( sqrt( h2_theta *k_scale/  var( G_theta  ) ) )
        G_theta  <- G %*% theta
      }
      
      U_faix <- U %*% faix; U_faix <- U_faix / as.numeric( sqrt( var( U_faix )/0.8 ) ) ## 0.8
      U_faiy <- U %*% faiy; U_faiy <- U_faiy / as.numeric( sqrt( var( U_faiy )/0.8 ) ) 
      
      ##### Calculate coefficients in discovery data #####
      X1 <- G_delta + G_kappax + U_faix + rnorm( nsample )*as.numeric( sqrt( 1 - var( U_faix ) ) )
      X2 <- G_gamma2 + rnorm( nsample )
      model_s <- MR.CUE::fastSigLm( X1[1:ns], G[1:ns, ] )
      coef_r <- model_s$coef[(nc + 1):nsnp]
      PRS <- Gr %*% coef_r; PRS <- PRS / as.numeric( sqrt( var( PRS ) ) )
      PRS_max <- max(PRS); PRS_min <- min( PRS ); # var_U_faiy <- var( U_faiy )
      error_generate <- function( x, PRS_max, PRS_min ){
        rnorm( 1, 0, ( x - PRS_min )/ ( PRS_max - PRS_min ) + 0.05 )
      }
      error_y <- sapply( PRS, function(x) error_generate( x, PRS_max, PRS_min ) )
      Y <- beta1 * X1 + beta2 * X2 + G_kappay + G_theta + U_faiy + error_y*as.numeric( sqrt(1 - var(U_faiy) ) )
      # var( beta1*G_delta )/var(Y); var( beta1*G_kappax )/var(Y); var( G_kappay )/var(Y); var( G_theta )/var(Y); var( G_gamma2 )/var(Y)
      
      model_gamma1 <- MR.CUE::fastSigLm( X1[(ns+1):nsample],  G[(ns+1):nsample, ] )
      model_gamma2 <- MR.CUE::fastSigLm( X2[(ns+1):nsample],  G[(ns+1):nsample, ] )
      model_gammaP <- MR.CUE::fastSigLm( PRS[(ns+1):nsample], G[(ns+1):nsample, ] )
      model_Gamma  <- MR.CUE::fastSigLm( Y[(ns+1):nsample],   G[(ns+1):nsample, ] )
      
      mydata <- cbind( bx = cbind( model_gamma1$coef, model_gamma2$coef, model_gammaP$coef ), 
                       bxse = cbind( model_gamma1$std,  model_gamma2$std, model_gammaP$std ), 
                       by = model_Gamma$coef, 
                       byse = model_Gamma$std )
      colnames(mydata) <- c( paste0('bx.', 1:3), paste0('bxse.', 1:3), 'by', 'byse')
      # mydata[990:1020, ]
      return( mydata )
    }
    
    mydata <- data_generating( h2_delta = h2_delta, h2_kappax = h2_kappax, beta1 = beta1_true )
    
    ###### Calculate results  #####
    naive_mvinput <- mr_mvinput( bx   = mydata[ ,grep( 'bx\\.', colnames(mydata) )][, -1], 
                                 bxse = mydata[ ,grep( 'bxse\\.', colnames(mydata) )][, -1],
                                 by   = as.numeric( mydata[ , colnames(mydata) == 'by'] ), 
                                 byse = as.numeric( mydata[ , colnames(mydata) == 'byse'] ) )
    mvinput <- mr_mvinput( bx   = mydata[ ,grep( 'bx\\.', colnames(mydata) )], 
                           bxse = mydata[ ,grep( 'bxse\\.', colnames(mydata) )],
                           by   = as.numeric( mydata[ , colnames(mydata) == 'by'] ), 
                           byse = as.numeric( mydata[ , colnames(mydata) == 'byse'] ) )
    
    res_naiveGibbs <- adapt( mydata[ ,-c( 3, 6) ], 1000 )
    res_Gibbs <- adapt( mydata, 1000 )
    res_naiveivw <- mr_mvivw( naive_mvinput )
    res_ivw <- mr_mvivw( mvinput )
    res_naiveegger <- mr_mvegger( naive_mvinput )
    res_egger <- mr_mvegger( mvinput )
    res_naivelasso <- mr_mvlasso( naive_mvinput )
    res_lasso <- mr_mvlasso( mvinput )
    
    naiveGibbs <- c( res_naiveGibbs$beta1['Est'], res_naiveGibbs$beta1['LowerCI'], res_naiveGibbs$beta1['UpperCI'] )
    mvGibbs    <- c( res_Gibbs$beta1['Est'],      res_Gibbs$beta1['LowerCI'],      res_Gibbs$beta1['UpperCI'] )
    naivemvivw <- c( res_naiveivw@Estimate[1],    res_naiveivw@CILower[1],         res_naiveivw@CIUpper[1] )
    mvivw      <- c( res_ivw@Estimate[1],         res_ivw@CILower[1],              res_ivw@CIUpper[1] ) 
    naiveegger <- c( res_naiveegger@Estimate[1],  res_naiveegger@CILower.Est[1],   res_naiveegger@CIUpper.Est[1] )
    mvegger    <- c( res_egger@Estimate[1],       res_egger@CILower.Est[1],        res_egger@CIUpper.Est[1] )
    naivelasso <- c( res_naivelasso@Estimate[1],  res_naivelasso@CILower[1],       res_naivelasso@CIUpper[1] )
    mvlasso    <- c( res_lasso@Estimate[1],       res_lasso@CILower[1],            res_lasso@CIUpper[1] )
    
    
    ##### Type I error #####
    typeI_naiveGibbs  <- ifelse( naiveGibbs[2] > 0 | naiveGibbs[3] < 0, 'Reject', 'Accept' )
    typeI_Gibbs       <- ifelse( mvGibbs[2]    > 0 | mvGibbs[3]    < 0, 'Reject', 'Accept' )
    typeI_naiveivw    <- ifelse( naivemvivw[2] > 0 | naivemvivw[3] < 0, 'Reject', 'Accept' )
    typeI_ivw         <- ifelse( mvivw[2]      > 0 | mvivw[3]      < 0, 'Reject', 'Accept' )
    typeI_naiveegger  <- ifelse( naiveegger[2] > 0 | naiveegger[3] < 0, 'Reject', 'Accept' )
    typeI_egger       <- ifelse( mvegger[2]    > 0 | mvegger[3]    < 0, 'Reject', 'Accept' )
    typeI_naivelasso  <- ifelse( naivelasso[2] > 0 | naivelasso[3] < 0, 'Reject', 'Accept' )
    typeI_lasso       <- ifelse( mvlasso[2]    > 0 | mvlasso[3]    < 0, 'Reject', 'Accept' )
    typeIerror        <- c( typeI_naiveGibbs, typeI_Gibbs, typeI_naiveivw, typeI_ivw, typeI_naiveegger, typeI_egger, typeI_naivelasso, typeI_lasso )
    names(typeIerror) <- c('naive_Gibbs', 'Gibbs', 'naive_ivw', 'ivw', 'naive_egger', 'egger', 'naive_lasso', 'lasso')
    
    ##### Bias #####
    bias_naiveGibbs <- naiveGibbs[1] - beta1_true
    bias_Gibbs      <- mvGibbs[1]    - beta1_true
    bias_naiveivw   <- naivemvivw[1] - beta1_true
    bias_ivw        <- mvivw[1]      - beta1_true
    bias_naiveegger <- naiveegger[1] - beta1_true
    bias_egger      <- mvegger[1]    - beta1_true
    bias_naivelasso <- naivelasso[1] - beta1_true
    bias_lasso      <- mvlasso[1]    - beta1_true
    bias            <- c( bias_naiveGibbs, bias_Gibbs, bias_naiveivw, bias_ivw, bias_naiveegger, bias_egger, bias_naivelasso, bias_lasso )
    names(bias)     <- c('naive_Gibbs', 'Gibbs', 'naive_ivw', 'ivw', 'naive_egger', 'egger', 'naive_lasso', 'lasso')
    
    return( list( typeIerror = typeIerror, bias = bias )  )
    
  }
  
  
  ##### Parallel Computing #####
  library(parallel)
  mycores <- 30
  n_parallel <- 300
  cl <- makeCluster( mycores, setup_strategy = 'parallel' )
  
  temp <- parLapply( cl, 1:n_parallel, comparison, h2_delta, h2_kappax )
  
  save( temp, file = paste0( 'temp_', k_gamma1, '.RData' )  )
  
  stopCluster(cl)
  
  ########## Save results ##########
  typeI <- data.frame( matrix(1:8, ncol = 8) )
  bias  <- data.frame( matrix(1:8, ncol = 8) )
  for( i in 1:length( temp ) ){
    res <- matrix( unlist( temp[[i]] ), ncol = 8, byrow = T)
    typeI[i, ] <- res[1, ]
    bias[i, ]  <- res[2, ]
  }
  
  colnames(typeI) <- c('naive_Gibbs', 'Gibbs', 'naive_ivw', 'ivw', 'naive_egger', 'egger', 'naive_lasso', 'lasso')
  colnames(bias)  <- c('naive_Gibbs', 'Gibbs', 'naive_ivw', 'ivw', 'naive_egger', 'egger', 'naive_lasso', 'lasso')
  
  write.csv(typeI,  paste0('/home/user2/typeI_', k_gamma1, '_200.csv'), row.names = F)
  write.csv(bias,   paste0('/home/user2/bias_',  k_gamma1, '_200.csv'), row.names = F)
  
  gc()
  
}