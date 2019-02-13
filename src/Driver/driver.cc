// Comparison of the super-fermion approach with Landauer predictions

#include <iostream>
#include <vector>
#include <cmath>
#include <sys/time.h>

#include "mkl.h"
#include "../Utils/Utils.h"

double seconds()
{
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}

std::vector<double> linspace(double a, double b, size_t n) {
    double h = (b - a) / static_cast<double>(n - 1);
    std::vector<double> vec(n);
    typename std::vector<double>::iterator x;
    double val;
    for (x = vec.begin(), val = a; x != vec.end(); ++x, val += h)
        *x = val;
    return vec;
}

int main(int argc, char **argv)
{
  // Properties of the leads
  double pi = M_PI;
  double w = 5.0;
  double Emin = -1.0 * w;
  double Emax =  1.0 * w;
  MKL_INT N = 200; // Number of eigenmodes in each lead
  std::vector<double> en = linspace(Emin, Emax, N);
  double d_en = en[1] - en[0];
  double mu_l = 0.5; // chemical potential left
  double mu_r = -0.5; // chemical potential right
  double t_l = 0.1; // temperature left
  double t_r = 0.1; // temperature right 

  // Gamma from LB
  double Gamma = 1.0;
  // Hopping parameter between dot and eigenmodes
  double t = std::sqrt( (d_en * Gamma) / (2.0 * pi) );
  // Thermalisation rate
  double therm = 2.0 * d_en;

  double tic = seconds();

  // Fermi-Dirac Distributions
  std::vector<double> f_l(N);
  std::vector<double> f_r(N);
  for(MKL_INT i = 0; i < N; ++i){
    f_l[i] = 1.0 / (1.0 + (std::exp( (en[i] - mu_l) / t_l )));
    f_r[i] = 1.0 / (1.0 + (std::exp( (en[i] - mu_r) / t_r )));
  }
  // Thermalisation excitation and decay rates
  MKL_INT h_size = (2 * N) + 1;
  std::vector<double> Ge(h_size);
  std::vector<double> Gd(h_size);
  Ge[0] = 0.0; Gd[0] = 0.0;
  for(MKL_INT i = 0; i < N; ++i){
    Ge[i + 1] = therm * f_l[i];
    Ge[N + i + 1] = therm * f_r[i];
    Gd[i + 1] = therm * (1.0 - f_l[i]);
    Gd[N + i + 1] = therm * (1.0 - f_r[i]);
  }

  // Dot energy
  MKL_INT samp = 20;
  std::vector<double> epsilon = linspace(-1.8 * w, 1.8 * w, samp);

  std::vector<double> dens_sf(samp, 0.0);
  std::vector<double> curr_sf(samp, 0.0);
  std::vector<double> ener_sf_dot(samp, 0.0);
  std::vector<double> ener_sf_lead_left(samp, 0.0);
  std::vector<double> ener_sf_lead_right(samp, 0.0);
  for(MKL_INT sp = 0; sp < samp; ++sp){
    
    // Superfermion solution
    // Construct the entire system Hamiltonian: dot + left leaf + right lead
    std::vector<double> H( h_size * h_size, 0.0 );
    for(MKL_INT i = 0; i < N; ++i){
      H[(h_size * (i + 1)) + (i + 1)] = en[i]; 
      H[(h_size * (N + i + 1)) + (N + i + 1)] = en[i]; 
    }
    for(MKL_INT i = 1; i < h_size; ++i){
      H[i] = -1.0 * t; 
      H[(h_size * i)] = -1.0 * t; 
    }
    H[0] = epsilon[sp];
  
    // Superfermion solution to the correlation matrix
    std::vector< std::complex<double> > corr;
    corr = Utils::super_fermion_solve(H, Ge, Gd, h_size);

    // Dot density
    double n = real(corr[0]);

    // Current
    std::complex<double> im(0.0, 1.0);
    std::vector<double> curr(h_size * h_size);
    for(MKL_INT i = 0; i < h_size; ++i){
      for(MKL_INT j = 0; j < h_size; ++j){
        curr[(i * h_size) + j] = real(-1.0 * im * H[(i * h_size) + j]
          * (corr[(i * h_size) + j] - corr[(j * h_size) + i]));
      }
    }
    
    double j_left = 0.0;
    double j_right = 0.0;
    for(MKL_INT i = 1; i < (N + 1); ++i){
      j_left += curr[i];
      j_right += curr[N + i];
    }

    double ener_left_dot = epsilon[sp] * j_left;
    double ener_right_dot = epsilon[sp] * j_right;

    double ener_lead_left = 0.0;
    double ener_lead_right = 0.0;
    for(MKL_INT i = 1; i < (N + 1); ++i){
      ener_lead_left += curr[i] * en[i - 1];
      ener_lead_right += curr[N + i] * en[i - 1];
    }

    dens_sf[sp] = n;
    curr_sf[sp] = j_left;
    ener_sf_dot[sp] = ener_left_dot;
    ener_sf_lead_left[sp] = ener_lead_left;
    ener_sf_lead_right[sp] = ener_lead_right;
  }

  // Landauer-Buttiker solutions
  MKL_INT wsamp = 10000;
  std::vector<double> w_lb = linspace(-1.3 * w, 1.3 * w, wsamp);
  double dw = w_lb[1] - w_lb[0];
  std::vector<double> box(wsamp, 0.0);
  for(MKL_INT i = 0; i < wsamp; ++i){
    if( (w_lb[i] >= Emin) && (w_lb[i] <= Emax) ) 
      box[i] = Gamma;
  }

  // Re-evaluate Fermi-Dirac distros for left-right
  std::vector<double> f_l_lb(wsamp);
  std::vector<double> f_r_lb(wsamp);
  for(MKL_INT i = 0; i < wsamp; ++i){
    f_l_lb[i] = 1.0 / (1.0 + (std::exp( (w_lb[i] - mu_l) / t_l ))) * box[i];
    f_r_lb[i] = 1.0 / (1.0 + (std::exp( (w_lb[i] - mu_r) / t_r ))) * box[i];
  }

  // Analytic expressions for real and imaginary self-energy components
  std::vector<double> gamma(box);
  std::vector<double> lambd(wsamp, 0.0);
  for(MKL_INT i = 0 ; i < wsamp; ++i)
    lambd[i] = (gamma[i] / (2.0 * pi)) * std::log( std::fabs( (w_lb[i] - Emin) / (w_lb[i] - Emax) ) );
 
  std::vector<double> dens_lb(samp, 0.0);
  std::vector<double> curr_lb(samp, 0.0);
  std::vector<double> ener_lb(samp, 0.0);
  for(MKL_INT i = 0; i < samp; ++i){
 
    // Density and current
    for(MKL_INT j = 0; j < wsamp; ++j){
      dens_lb[i] += (1.0 / (2.0 * pi)) * dw * ( gamma[j] * ( f_l_lb[j] + f_r_lb[j] ) )
        / (std::pow(w_lb[j] - epsilon[i] - (2.0 * lambd[j]), 2) + std::pow(gamma[j], 2));
      curr_lb[i] += (1.0 / (2.0 * pi)) * dw * ( gamma[j] * gamma[j] * 
        ( f_l_lb[j] - f_r_lb[j] ) )
          / (std::pow(w_lb[j] - epsilon[i] - (2.0 * lambd[j]), 2) + std::pow(gamma[j], 2));
      ener_lb[i] += (1.0 / (2.0 * pi)) * dw * ( w_lb[j] * gamma[j] * gamma[j] * 
        ( f_l_lb[j] - f_r_lb[j] ) )
          / (std::pow(w_lb[j] - epsilon[i] - (2.0 * lambd[j]), 2) + std::pow(gamma[j], 2));
    }
  }

  double toc = seconds();

  // Print solutions
  std::cout << std::fixed;
  std::cout.precision(8);
  std::cout << "# N = " << N << std::endl;
  std::cout << "# Gamma = " << Gamma << std::endl;
  std::cout << "# d_en = " << d_en << std::endl;
  std::cout << "# therm = " << therm << std::endl;
  std::cout << "# t_k = " << t << std::endl;
  std::cout << "# mu_l = " << mu_l << " mu_r = " << mu_r << std::endl;
  std::cout << "# t_l = " << t_l << " t_r = " << t_r << std::endl;
  std::cout << "# E_min = " << Emin << " E_max = " << Emax << std::endl;
  std::cout << "# Time = " << toc - tic << std::endl;
  std::cout << "# Dot_En" << " " << "Dot_n_LB" << " " << "Dot_n_SF" << " " << 
    "J_left_LB" << " " << "J_left_SF" << " " << "J_ener_LB" << " " << "J_ener_SF_dot" <<
      " " << "J_ener_SF_lead_left" << " " << "J_ener_SF_lead_right" << std::endl;
  for(MKL_INT sp = 0; sp < samp; ++sp){
    std::cout << epsilon[sp] << " " << dens_lb[sp] << " " << dens_sf[sp] << " " <<
      curr_lb[sp] << " " << curr_sf[sp] << " " << ener_lb[sp] << " " << ener_sf_dot[sp] << 
        " " << ener_sf_lead_left[sp] << " " << ener_sf_lead_right[sp] << std::endl; 
  }

  return 0;
}
