#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/FFT>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/constants/constants.hpp>

using namespace boost::multiprecision;
using boost::math::constants::pi;

// [[Rcpp::depends(BH,RcppEigen,Rmpfr)]]
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector newfft3(NumericVector k, NumericVector lambda, NumericVector ncps, double a, double b, int n) {
  // Write wrapper with checks that inputs are sane
  // eg length(lambda) == length(ncps) etc
  NumericVector res(k.length());
  int nevals = lambda.length();

  #if 0
  typedef double FFTFloat;
  #else
  //typedef number<cpp_dec_float<50> > cpp_dec_float_50_ljma;
  //typedef cpp_dec_float_50_ljma FFTFloat;
  typedef number<mpfr_float_backend<50, allocate_dynamic>, et_off> FFTFloat;
  #endif

  int L = k.length();               // Length of signal
  typedef std::complex<FFTFloat> FFTFloatC;
  typedef std::vector<FFTFloat> FFTFloatV;
  typedef std::vector<FFTFloatC> FFTFloatCV;
  FFTFloatCV x(L);

Rcout << "Start consts\n";
  // Some constants
  FFTFloat c1, c2, c3, one(1.0), zero(0.0);
  c1 = pi<FFTFloat>()*FFTFloat(n)/(FFTFloat(b)-FFTFloat(a));
Rcout << "c1 = " << c1 << "\n";
  c2 = FFTFloat(2.0)/FFTFloat(n);
Rcout << "c2 = " << c2 << "\n";
  c3 = pi<FFTFloat>()*FFTFloat(a)*FFTFloat(n)/(FFTFloat(b)-FFTFloat(a));
Rcout << "c3 = " << c3 << "\n";

  FFTFloatV nl(nevals), l(nevals);
Rcout << "nl, l\n";
  for(int i = 0; i < nevals; i++) {
    nl[i] = FFTFloat(ncps[i])*FFTFloat(lambda[i]);
    l[i] = FFTFloat(-2.0*lambda[i]);
  }

  FFTFloatC p5(FFTFloat(0.5), zero);
Rcout << "p5 = " << p5 << "\n";
Rcout << "End consts\n";

  FFTFloatC log_phi;
  FFTFloat t, tl;
  for(int i = 0; i < L; i++) {
if(i<5 || i>65500) Rcout << "Start rho[" << i << "]\n";
    t = c1*(c2*FFTFloat(k[i]-1.0)-one);
if(i<5 || i>65500) Rcout << "t = " << t << "\n";

    log_phi.real(zero);
    log_phi.imag(zero);
    for(int j = 0; j < nevals; j++) {
      tl = l[j] * t;
      log_phi += FFTFloatC(zero, nl[j]*t)/FFTFloatC(one, tl) - p5 * log(FFTFloatC(one, tl));
    }
if(i<5 || i>65500) Rcout << "log_phi = " << log_phi << "\n";

    x[i] = exp( log_phi - FFTFloatC(zero, c3*(c2*FFTFloat(k[i]-1.0)-one)) );
if(i<5 || i>65500) Rcout << "rho[" << i << "] = " << x[i] << "\n";
  }

  FFTFloatCV y;
  Eigen::FFT<FFTFloat> fft;
  fft.fwd(y, x);
  for (int i = 0; i < y.size(); i++) {
    res[i] = y[i].real().convert_to<double>();
    //res[i] = y[i].real();
  }

  return res;
}
