#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/FFT>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/constants/constants.hpp>
#include <thread>
#include <algorithm>
#include <memory>

using namespace boost::multiprecision;
using boost::math::constants::pi;

// [[Rcpp::depends(BH,RcppEigen,Rmpfr)]]
using namespace Rcpp;

#if 0
typedef double FFTFloat;
#else
typedef number<cpp_dec_float<50>> FFTFloat;
//typedef number<mpfr_float_backend<50, allocate_dynamic>, et_on> FFTFloat;
#endif
typedef std::complex<FFTFloat> FFTFloatC;
typedef std::vector<FFTFloat> FFTFloatV;
typedef std::vector<FFTFloatC> FFTFloatCV;


void computeRho(FFTFloatCV& rho, int from_i, int to_i, NumericVector& k, FFTFloatV& nl, FFTFloatV& l,
                FFTFloat& c1, FFTFloat& c2, FFTFloat& c3, FFTFloat& zero, FFTFloat& one, FFTFloatC& p5,
                int nevals) {
  FFTFloatC log_phi;
  FFTFloat t, tl;
  for(int i = from_i; i <= to_i; i++) {
    // if(i<5 || i>65500) Rcout << "Start rho[" << i << "]\n";
    t = c1*(c2*FFTFloat(k[i]-1.0)-one);
    // if(i<5 || i>65500) Rcout << "t = " << t << "\n";

    log_phi.real(zero);
    log_phi.imag(zero);
    for(int j = 0; j < nevals; j++) {
      tl = l[j] * t;
      log_phi += FFTFloatC(zero, nl[j]*t)/FFTFloatC(one, tl) - p5 * log(FFTFloatC(one, tl));
    }
    // if(i<5 || i>65500) Rcout << "log_phi = " << log_phi << "\n";

    rho[i] = exp( log_phi - FFTFloatC(zero, c3*(c2*FFTFloat(k[i]-1.0)-one)) );
    // if(i<5 || i>65500) Rcout << "rho[" << i << "] = " << x[i] << "\n";
  }
}

// [[Rcpp::export]]
NumericVector newfft4(NumericVector k, NumericVector lambda, NumericVector ncps, double a, double b, int n, int nthreads) {
  // Write wrapper with checks that inputs are sane
  // eg length(lambda) == length(ncps) etc
  int L = k.length();               // Length of signal

  NumericVector res(k.length());
  int nevals = lambda.length();

  // FFT input
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

  std::vector<std::thread> threads;
  int by = ceil(((double) L)/((double) nthreads));
  for(int th = 0; th < std::min(nthreads, L); th++) {
    threads.push_back(std::thread(computeRho,
                                  std::ref(x), th*by, (th+1)*by-1, std::ref(k), std::ref(nl), std::ref(l),
                                  std::ref(c1), std::ref(c2), std::ref(c3), std::ref(zero), std::ref(one), std::ref(p5),
                                  nevals));
  }
  for(auto &th : threads) {
    th.join();
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
