#include <regex>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double calculateP( double mean, double sd );
List RARE( const arma::mat& mydata, int iter_times = 1000 );


// [[Rcpp::export]]
List RARE( const arma::mat& mydata, int iter_times ){
  
  int nsnps = mydata.n_rows;  int n_exposure = mydata.n_cols/2 - 1;
  
  arma::mat gamma( nsnps, n_exposure, arma::fill::value(0.1) );
  arma::vec s2g( n_exposure, arma::fill::value(0.1) );
  
  arma::mat gammah = mydata.cols( 0, n_exposure - 1 );
  arma::mat s2hg   = mydata.cols( n_exposure, 2*n_exposure - 1 ) % mydata.cols( n_exposure, 2*n_exposure - 1 );
  
  arma::vec Gamma( nsnps, arma::fill::value(0.1) );
  arma::vec Gammah = mydata.col( 2*n_exposure );
  arma::vec s2hG   = mydata.col( 2*n_exposure + 1 ) % mydata.col( 2*n_exposure + 1 );
  
  arma::vec kappay( nsnps, arma::fill::value(0.1) ); double sigma_kappay_sq = 0.1; 
  double sigma_theta_sq = 0.1;
  
  arma::vec beta( n_exposure, arma::fill::value(0.1) );
  arma::mat beta_new( iter_times, n_exposure );
  
  double a_shape, b_scale; a_shape = 1; b_scale = 1; 
  
  double sigmat_Gamma, miut_Gamma, sigmat_gamma, miut_gamma, sigmat_kappay, miut_kappay, 
  sigmat_beta, miut_beta, beta_ran, at_gamma, bt_gamma,
  at_kappay, bt_kappay, at_theta, bt_theta;

  for( int iter = 0; iter < iter_times; ++iter ){
 
    for( int snp = 0; snp < nsnps; ++snp ){
      ////////// Start inner circulation //////////
      arma::rowvec gammav = gamma.row( snp ); arma::rowvec gammahv = gammah.row( snp ); arma::rowvec s2hgv = s2hg.row( snp );
      /// Gamma ///
      sigmat_Gamma = 1/( 1/s2hG( snp ) + 1/sigma_theta_sq );
      miut_Gamma   = ( Gammah( snp )/s2hG( snp ) + ( arma::dot( gammav, beta ) + kappay( snp ) )/sigma_theta_sq ) * sigmat_Gamma;
      Gamma(snp) = R::rnorm( miut_Gamma, sqrt( sigmat_Gamma ) );
      
      for( int exposure = 0; exposure < n_exposure; ++exposure ){
        arma::rowvec gammav_copy = gammav; arma::vec beta_copy = beta;
        gammav_copy.shed_col( exposure ); beta_copy.shed_row( exposure );
        
        sigmat_gamma = 1/( 1/s2hgv( exposure ) + beta( exposure )*beta( exposure )/sigma_theta_sq + 1/s2g( exposure ) );
        miut_gamma   = ( gammahv( exposure )/s2hgv( exposure ) + beta( exposure )*( Gamma( snp ) - arma::dot( gammav_copy, beta_copy ) - kappay( snp ) )/sigma_theta_sq )*sigmat_gamma;
        gamma( snp, exposure ) = R::rnorm( miut_gamma, sqrt( sigmat_gamma ) );
      }
      
      /// kappay ///
      sigmat_kappay = 1/( 1/sigma_theta_sq + 1/sigma_kappay_sq );
      miut_kappay   =  ( Gamma( snp ) - arma::dot( gammav, beta ) )/sigma_theta_sq *sigmat_kappay; 
      kappay( snp ) = R::rnorm( miut_kappay, sqrt( sigmat_kappay ) );
      ////////// End inner circulation //////////
    } 
    
    /// Update beta ///
    for( int exposure = 0; exposure < n_exposure; ++exposure ){
      arma::mat gamma_copy = gamma; arma::vec beta_copy = beta;
      gamma_copy.shed_col( exposure ); beta_copy.shed_row( exposure);
      
      sigmat_beta = 1/( arma::dot( gamma.col(exposure), gamma.col(exposure) )/sigma_theta_sq );
      miut_beta   = arma::dot( ( Gamma - gamma_copy * beta_copy - kappay ), gamma.col( exposure ) )/sigma_theta_sq * sigmat_beta;  
      
      beta_ran    = R::rnorm( miut_beta, sqrt( sigmat_beta ) );
      beta( exposure ) = beta_ran; beta_new( iter, exposure ) = beta_ran;
      
      at_gamma  = a_shape + nsnps/2; bt_gamma = b_scale + arma::dot( gamma.col( exposure ), gamma.col( exposure ) )/2;
      s2g( exposure ) = 1/as<double>( Rcpp::rgamma( 1, at_gamma, 1.0/bt_gamma) ); 
    }
    
    at_kappay = a_shape + nsnps/2; bt_kappay = b_scale + arma::dot( kappay, kappay )/2;
    sigma_kappay_sq = 1/as<double>( Rcpp::rgamma( 1, at_kappay, 1.0/bt_kappay) );
    
    at_theta = a_shape + nsnps/2; bt_theta = b_scale + arma::dot( Gamma- gamma*beta - kappay, Gamma- gamma*beta - kappay )/2;
    sigma_theta_sq = 1/as<double>( Rcpp::rgamma( 1, at_theta, 1.0/bt_theta) );
  } 
  
  /// Output ///
  List results;
  for( int exposure = 0; exposure < n_exposure; ++exposure ){
    
    std::string varName;
    varName = "beta" + std::to_string( exposure + 1 );
    
    // if( exposure == ( n_exposure - 1 ) ){
    //   varName = "betaPRS";
    // }else{
    //   varName = "beta" + std::to_string( exposure + 1 );
    // }
    
    double beta_mean = arma::mean( beta_new.col( exposure ).subvec( iter_times/2, iter_times - 1 ) );
    double beta_sd   = arma::stddev( beta_new.col( exposure ).subvec( iter_times/2, iter_times - 1 ) );
    double beta_p    = calculateP( beta_mean, beta_sd );
    
    NumericVector  result = NumericVector::create( beta_mean, beta_sd, beta_mean - 1.96*beta_sd, beta_mean + 1.96*beta_sd, beta_p ); 
    result.names() = CharacterVector::create( "Est", "Std", "LowerCI", "UpperCI", "P value" );
    
    results[varName] = result;
    
  }
  return results;
  
} 


// [[Rcpp::export]]
double calculateP( double mean, double sd ){
  double LCI = mean - 1.96*sd; double UCI = mean + 1.96*sd;
  double pvalue, pLower, pUpper;
  if( LCI >= 0 ){
    pvalue = R::pnorm( 0, mean, sd, true,  false );
  } else if( UCI <= 0 ){
    pvalue = R::pnorm( 0, mean, sd, false, false );
  } else{
    pLower = R::pnorm( 0, mean, sd, true,  false );
    pUpper = R::pnorm( 0, mean, sd, false, false );
    pvalue = 2*std::min( pLower, pUpper );
  }
  return pvalue;
}