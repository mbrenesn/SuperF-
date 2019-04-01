// Comparison of the super-fermion approach with Landauer predictions for a single dot

// Log-Linear lead discretisation 

#include <iostream>
#include <vector>
#include <algorithm>
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
  MKL_INT N_lin = 777;
  MKL_INT N_log = 777;
  double temperature = 0.777;
  double Gamma = 0.777;
  double mu = 0.777;
  double V = 0.777;

  if(argc != 13){
    std::cerr << "Usage: " << argv[0] << " --N_lin [lin discretised modes] --N_log [log discretised modes PER SIDE] --temperature [T] --Gamma [Big Gamma] --mu [mu] --V [V]" << std::endl;
    exit(1);
  }
  for(int i = 0; i < argc; ++i){
    std::string str = argv[i];
    if(str == "--N_lin") N_lin = atoi(argv[i + 1]);
    else if(str == "--N_log") N_log = atoi(argv[i + 1]);
    else if(str == "--temperature") temperature = atof(argv[i + 1]);
    else if(str == "--Gamma") Gamma = atof(argv[i + 1]);
    else if(str == "--mu") mu = atof(argv[i + 1]);
    else if(str == "--V") V = atof(argv[i + 1]);
    else continue;
  }
  if(N_lin == 777 || N_log == 777 || temperature == 0.777 || Gamma == 0.777 || mu == 0.777 || V == 0.777){
    std::cerr << "Error setting parameters" << std::endl;
    std::cerr << "Usage: " << argv[0] << " --N_lin [lin discretised modes] --N_log [log discretised modes PER SIDE] --temperature [T] --Gamma [Big Gamma] --mu [mu] --V [V]" << std::endl;
    exit(1);
  }

  double tic = seconds();
  // Properties of the leads
  double mu_l = mu + (V / 2.0); // chemical potential left
  double mu_r = mu - (V / 2.0); // chemical potential right
  double t_l = temperature; // temperature left
  double t_r = temperature; // temperature right 
  // Here, we use a Log-Linear discretisation of the leads, N_lin are discretised linearly
  // in the vicinity of epsilon and N_log are discretised logarithmically beyond N_lin
  MKL_INT N = N_lin + N_log + N_log;
  double base_log = 2.0;
  double pi = M_PI;
  double w = 8.0;
  double Emin = -1.0 * w;
  double Emax =  1.0 * w;
  double Estar_min = -0.5 * w;
  double Estar_max =  0.5 * w;

  // Setting thermalisation rate = spacing in the modes
  // Linearly spaced portion
  std::vector<double> linear_spacings = linspace(Estar_min, Estar_max, N_lin + 1);
  std::vector<double> en_lin(N_lin, 0.0);
  for(MKL_INT i = 0; i < N_lin; ++i){
    en_lin[i] = (linear_spacings[i] + linear_spacings[i + 1]) / 2.0;
  }
  std::vector<double> therm_lin(N_lin, 0.0);
  for(MKL_INT i = 0; i < N_lin; ++i){
    therm_lin[i] = std::abs(linear_spacings[i] - linear_spacings[i + 1]);
  }

  // Logarithmically spaced portion: Base = 2
  // Right side
  std::vector<double> log_spacings_right = linspace(std::log2(Estar_max), std::log2(Emax), N_log + 1);
  for(MKL_INT i = 0; i < (N_log + 1); ++i){
    log_spacings_right[i] = std::pow(2.0, log_spacings_right[i]);
  }
  std::vector<double> en_log_right(N_log, 0.0);
  for(MKL_INT i = 0; i < N_log; ++i){
    en_log_right[i] = (log_spacings_right[i] + log_spacings_right[i + 1]) / 2.0;
  }
  std::vector<double> therm_log_right(N_log, 0.0);
  for(MKL_INT i = 0; i < N_log; ++i){
    therm_log_right[i] = std::abs(log_spacings_right[i] - log_spacings_right[i + 1]);
  }
  // Left side: equal to right side but reverted / sign changed
  std::vector<double> en_log_left(N_log, 0.0);
  for(MKL_INT i = 0; i < N_log; ++i){
    en_log_left[i] = -1.0 * en_log_right[N_log - 1 - i];
  }
  std::vector<double> therm_log_left(N_log, 0.0);
  for(MKL_INT i = 0; i < N_log; ++i){
    therm_log_left[i] = std::abs(log_spacings_right[N_log - 1 - i] - log_spacings_right[N_log - i]);
  }

  std::vector<double> en;
  en.insert(en.end(), en_log_left.begin(), en_log_left.end());
  en.insert(en.end(), en_lin.begin(), en_lin.end());
  en.insert(en.end(), en_log_right.begin(), en_log_right.end());

  std::vector<double> therm;
  therm.insert(therm.end(), therm_log_left.begin(), therm_log_left.end());
  therm.insert(therm.end(), therm_lin.begin(), therm_lin.end());
  therm.insert(therm.end(), therm_log_right.begin(), therm_log_right.end());

  // Hopping parameter between dot and eigenmodes. (therm = d_en)
  std::vector<double> t(N, 0.0);
  for(MKL_INT i = 0; i < N; ++i){
    t[i] = -1.0 * std::sqrt( (therm[i] * Gamma) / (2.0 * pi) );
  }

  // Dot energies to sample
  MKL_INT samp = 20;
  std::vector<double> epsilon = linspace(-6.5, 6.5, samp);
  //std::vector<double> epsilon(1, 1.0);
  
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
    Ge[i + 1] = 0.5 * therm[i] * f_l[i];
    Ge[N + i + 1] = 0.5 * therm[i] * f_r[i];
    Gd[i + 1] = 0.5 * therm[i] * (1.0 - f_l[i]);
    Gd[N + i + 1] = 0.5 * therm[i] * (1.0 - f_r[i]);
  }

  std::vector<double> dens_sf(samp, 0.0);
  std::vector<double> curr_sf(samp, 0.0);
  std::vector<double> ener_sf(samp, 0.0);
  for(MKL_INT sp = 0; sp < samp; ++sp){  
    // Superfermion solution
    // Construct the entire system Hamiltonian: dot + left leaf + right lead
    std::vector<double> H( h_size * h_size, 0.0 );
    for(MKL_INT i = 0; i < N; ++i){
      H[(h_size * (i + 1)) + (i + 1)] = en[i]; 
      H[(h_size * (N + i + 1)) + (N + i + 1)] = en[i]; 
    }
    for(MKL_INT i = 1; i < (N + 1); ++i){
      H[i] = t[i - 1];
      H[i + N] = t[i - 1];
      H[(h_size * i)] = t[i - 1];
      H[(h_size * i) + (N * h_size)] = t[i - 1];
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
        curr[(i * h_size) + j] = std::real(-1.0 * im * H[(i * h_size) + j]
          * (corr[(i * h_size) + j] - corr[(j * h_size) + i]));
      }
    }
    
    double j_left = 0.0;
    double j_right = 0.0;
    for(MKL_INT i = 1; i < (N + 1); ++i){
      j_left += curr[i];
      j_right += curr[N + i];
    }

    // Energy
    std::vector<double> cross(h_size * h_size);
    for(MKL_INT i = 0; i < h_size; ++i){
      for(MKL_INT j = 0; j < h_size; ++j){
        cross[(i * h_size) + j] = std::real( -1.0 * H[(i * h_size) + j]
          * (corr[(i * h_size) + j] + corr[(j * h_size) + i]));
      }
    }

    double accum = 0.0;
    for(MKL_INT i = 1; i < (N + 1); ++i){
      double n_local = std::real( corr[(i * h_size) + i] );
      accum += ( cross[i] * therm[i - 1] * 0.25 ) + 
        ( 0.5 * therm[i - 1] * en[i - 1] * (f_l[i - 1] - n_local) ); 
    }

    dens_sf[sp] = n;
    curr_sf[sp] = j_left;
    ener_sf[sp] = accum;
  }

  // Landauer-Buttiker solutions
  MKL_INT wsamp = 10000;
  std::vector<double> w_lb = linspace(-1.3 * w, 1.3 * w, wsamp);
  double dw = w_lb[1] - w_lb[0];
  std::vector<double> box(wsamp, 0.0);
  for(MKL_INT i = 0; i < wsamp; ++i){
    if( (w_lb[i] >= Emin) && (w_lb[i] <= Emax) ) 
      box[i] = 1.0;
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
  for(MKL_INT i = 0 ; i < wsamp; ++i)
    gamma[i] = gamma[i] * Gamma;
  std::vector<double> lambd(wsamp, 0.0);
  for(MKL_INT i = 0 ; i < wsamp; ++i)
    lambd[i] = (gamma[i] / (2.0 * pi)) * std::log( std::fabs( (w_lb[i] - Emin) / (w_lb[i] - Emax) ) );
  
  std::vector<double> dens_lb(samp, 0.0);
  std::vector<double> curr_lb(samp, 0.0);
  std::vector<double> ener_lb(samp, 0.0);
  for(MKL_INT i = 0; i < samp; ++i){
 
    // Density and current
    for(MKL_INT j = 0; j < wsamp; ++j){
      double denom = std::pow( w_lb[j] - epsilon[i] - (2.0 * lambd[j]), 2.0) 
          + std::pow( gamma[j], 2.0);
      double nn = gamma[j] * ( f_l_lb[j] + f_r_lb[j] );
      double jj = gamma[j] * gamma[j] * ( f_l_lb[j] - f_r_lb[j] );
      double ee = w_lb [j] * gamma[j] * gamma[j] * ( f_l_lb[j] - f_r_lb[j] );

      dens_lb[i] += (1.0 / (2.0 * pi)) * ( nn / denom ) * dw;
      curr_lb[i] += (1.0 / (2.0 * pi)) * ( jj / denom ) * dw;
      ener_lb[i] += (1.0 / (2.0 * pi)) * ( ee / denom ) * dw;
    }
  }

  double toc = seconds();
  std::cout << std::fixed;
  std::cout.precision(8);

  //Power
  std::vector<double> power_lb(samp, 0.0);
  std::vector<double> power_sf(samp, 0.0);
  for(MKL_INT sp = 0; sp < samp; ++sp){
    power_lb[sp] = curr_lb[sp] * V;
    power_sf[sp] = curr_sf[sp] * V;
  }
  // Efficiency
  double carnot_ef = 1.0 - (t_r / t_l);
  std::vector<double> efficiency_lb(samp, 0.0);
  std::vector<double> efficiency_sf(samp, 0.0);
  for(MKL_INT sp = 0; sp < samp; ++sp){
    efficiency_lb[sp] = power_lb[sp] / ( ener_lb[sp] - (mu_l * curr_lb[sp]) );
    efficiency_sf[sp] = power_sf[sp] / ( ener_sf[sp] - (mu_l * curr_sf[sp]) );
    efficiency_lb[sp] = efficiency_lb[sp] / carnot_ef;
    efficiency_sf[sp] = efficiency_sf[sp] / carnot_ef;
    //std::cout << "Particle " << curr_lb[sp] << " " << curr_sf[sp] << std::endl;
    //std::cout << "Power " << power_lb[sp] << " " << power_sf[sp] << std::endl;
    //std::cout << "Energy " << ener_lb[sp] << " " << ener_sf_dot[sp] << std::endl;
  }
  
  //std::cout << "# mu_l = " << mu_l << " mu_r = " << mu_r << std::endl;
  //std::cout << "# t_l = " << t_l << " t_r = " << t_r << std::endl;
  //std::cout << "# Time = " << toc - tic << std::endl;
  //std::cout << "# mu" << " " << "V" << " " << "Power_LB" << " " << "Power_SF" << " " <<
  //  "Efficiency_LB" << " " << "Effciency_SF" << std::endl;
  for(MKL_INT sp = 0; sp < samp; ++sp){
    //std::cout << mu - epsilon[sp] << " " << V << " " << power_lb[sp] << " " <<
    //  power_sf[sp] << " " << efficiency_lb[sp] << " " << efficiency_sf[sp] << std::endl; 
    std::cout << epsilon[sp] << " " << dens_lb[sp] << " " << dens_sf[sp] << " " <<
      curr_lb[sp] << " " << curr_sf[sp] << " " << ener_lb[sp] << " " << ener_sf[sp] << std::endl; 
  }

  return 0;
}
