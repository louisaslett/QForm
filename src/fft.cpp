#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/FFT>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/mpfr.hpp>

using namespace boost::multiprecision;

// [[Rcpp::depends(BH,RcppEigen,Rmpfr)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector newfft(NumericVector x) {
  NumericVector res(x.length());

  #if 0
  typedef double FFTFloat;
  #else
  //typedef number<cpp_dec_float<50> > cpp_dec_float_50_ljma;
  //typedef cpp_dec_float_50_ljma FFTFloat;
  typedef number<mpfr_float_backend<50, allocate_dynamic>, et_off> FFTFloat;
  #endif
  int L = x.length();               // Length of signal
  std::vector<FFTFloat> timebuf(L);
  typedef std::complex<FFTFloat> C;
  std::vector<C> freqbuf;
  for (int i = 0; i < L; i++) {
    timebuf[i] = x[i];
  }
  Eigen::FFT<FFTFloat> fft;
  fft.fwd(freqbuf, timebuf);
  for (int i = 0; i < freqbuf.size(); i++) {
    res[i] = freqbuf[i].real().convert_to<double>();
  }

  return res;
}
