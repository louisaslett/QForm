#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/FFT>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/mpfr.hpp>

using namespace boost::multiprecision;

// [[Rcpp::depends(BH,RcppEigen,Rmpfr)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector newfft2(NumericVector re, NumericVector im) {
  NumericVector res(re.length());

  #if 0
  typedef double FFTFloat;
  #else
  //typedef number<cpp_dec_float<50> > cpp_dec_float_50_ljma;
  //typedef cpp_dec_float_50_ljma FFTFloat;
  typedef number<mpfr_float_backend<50, allocate_dynamic>, et_off> FFTFloat;
  #endif
  int L = re.length();               // Length of signal
  typedef std::complex<FFTFloat> C;
  std::vector<C> x(L);
  std::vector<C> y;
  for (int i = 0; i < L; i++) {
    x[i].real(re[i]);
    x[i].imag(im[i]);
  }
  Eigen::FFT<FFTFloat> fft;
  fft.fwd(y, x);
  for (int i = 0; i < y.size(); i++) {
    res[i] = y[i].real().convert_to<double>();
    //res[i] = y[i].real();
  }

  return res;
}
