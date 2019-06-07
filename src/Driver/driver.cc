// Comparison of the super-fermion approach with Landauer predictions

// Multiple fermionic sites in the central system

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
  MKL_INT L = 777;
  double temperature = 0.777;
  double Gamma = 0.777;
  double mu = 0.777;
  double V = 0.777;
  double epsilon_i = 0.777; // self-energy in each site
  double t_s = 0.777; //hopping within system

  if(argc != 19){
    std::cerr << "Usage: " << argv[0] << " --N_lin [lin discretised modes] --N_log [log discretised modes PER SIDE] --L [sites in system] --eps [On-site energies] --ts [System hoppings] --temperature [T] --Gamma [Big Gamma] --mu [mu] --V [V]" << std::endl;
    exit(1);
  }
  for(int i = 0; i < argc; ++i){
    std::string str = argv[i];
    if(str == "--N_lin") N_lin = atoi(argv[i + 1]);
    else if(str == "--N_log") N_log = atoi(argv[i + 1]);
    else if(str == "--L") L = atoi(argv[i + 1]);
    else if(str == "--eps") epsilon_i = atof(argv[i + 1]);
    else if(str == "--ts") t_s = atof(argv[i + 1]);
    else if(str == "--temperature") temperature = atof(argv[i + 1]);
    else if(str == "--Gamma") Gamma = atof(argv[i + 1]);
    else if(str == "--mu") mu = atof(argv[i + 1]);
    else if(str == "--V") V = atof(argv[i + 1]);
    else continue;
  }
  if(N_lin == 777 || N_log == 777 || L == 777 || temperature == 0.777 || Gamma == 0.777 || mu == 0.777 || V == 0.777 || epsilon_i == 0.777 || t_s == 0.777){
    std::cerr << "Error setting parameters" << std::endl;
    std::cerr << "Usage: " << argv[0] << " --N_lin [lin discretised modes] --N_log [log discretised modes PER SIDE] --L [sites in system] --temperature [T] --Gamma [Big Gamma] --mu [mu] --V [V]" << std::endl;
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

  // Vector containing the discretised energies of the leads
  std::vector<double> en;
  en.insert(en.end(), en_log_left.begin(), en_log_left.end());
  en.insert(en.end(), en_lin.begin(), en_lin.end());
  en.insert(en.end(), en_log_right.begin(), en_log_right.end());

  // Vector containing the thermalising rate for each mode in the leads, differs depending on the 
  // spacings
  std::vector<double> therm;
  therm.insert(therm.end(), therm_log_left.begin(), therm_log_left.end());
  therm.insert(therm.end(), therm_lin.begin(), therm_lin.end());
  therm.insert(therm.end(), therm_log_right.begin(), therm_log_right.end());

  // Hopping parameter between dot and eigenmodes. (therm = d_en)
  std::vector<double> t(N, 0.0);
  for(MKL_INT i = 0; i < N; ++i){
    t[i] = -1.0 * std::sqrt( (therm[i] * Gamma) / (2.0 * pi) );
  }

  // Fermi-Dirac Distributions
  std::vector<double> f_l(N);
  std::vector<double> f_r(N);
  for(MKL_INT i = 0; i < N; ++i){
    f_l[i] = 1.0 / (1.0 + (std::exp( (en[i] - mu_l) / t_l )));
    f_r[i] = 1.0 / (1.0 + (std::exp( (en[i] - mu_r) / t_r )));
  }
  // Thermalisation excitation and decay rates
  MKL_INT h_size = (2 * N) + L;
  std::vector<double> Ge(h_size, 0.0);
  std::vector<double> Gd(h_size, 0.0);
  for(MKL_INT i = 0; i < N; ++i){
    Ge[i + L] = therm[i] * f_l[i];
    Ge[N + i + L] = therm[i] * f_r[i];
    Gd[i + L] = therm[i] * (1.0 - f_l[i]);
    Gd[N + i + L] = therm[i] * (1.0 - f_r[i]);
  }
  
  // System parameters and system Hamiltonian
  MKL_INT samp = 1;

  // Superfermion solution
  // Construct the entire system Hamiltonian: system + left lead + right lead
  std::vector<double> H( h_size * h_size, 0.0 );
  if( L == 0 ) { std::cerr << "L can't be zero" << std::endl; exit(1); }
  else if( L == 1 ) H[0] = epsilon_i;
  else{
    for(MKL_INT i = 1; i < (L - 1); ++i){
      H[(i * h_size) + i] = epsilon_i;
      H[(i * h_size) + i - 1] = -1.0 * t_s;
      H[(i * h_size) + i + 1] = -1.0 * t_s;
    }
    H[0] = epsilon_i;
    H[1] = -1.0 * t_s;
    H[((L - 1) * h_size) + L - 1] = epsilon_i;
    H[((L - 1) * h_size) + L - 2] = -1.0 * t_s;
  }
  for(MKL_INT i = 0; i < N; ++i){
    H[(h_size * (i + L)) + (i + L)] = en[i]; 
    H[(h_size * (N + i + L)) + (N + i + L)] = en[i]; 
  }
  for(MKL_INT i = 0; i < N; ++i){
    H[i + L] = t[i];
    H[i + ((L - 1) * h_size) + (L + N)] = t[i];
    H[(h_size * (i + L) )] = t[i];
    H[(h_size * (i + L) ) + (N * h_size) + (L - 1)] = t[i];
  }

  std::vector<double> curr_sf(samp, 0.0);
  std::vector<double> ener_sf(samp, 0.0);
  for(MKL_INT sp = 0; sp < samp; ++sp){  
  
    // Superfermion solution to the correlation matrix
    std::vector< std::complex<double> > corr;
    corr = Utils::super_fermion_solve(H, Ge, Gd, h_size);

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
    for(MKL_INT i = 0; i < N; ++i){
      j_left += curr[i + L];
      j_right += curr[i + L + N];
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
    for(MKL_INT i = L; i < (N + L); ++i){
      double n_local = std::real( corr[(i * h_size) + i] );
      accum += ( cross[i] * therm[i - L] * 0.5 ) + 
        ( 1.0 * therm[i - L] * en[i - L] * (f_l[i - L] - n_local) ); 
    }

    curr_sf[sp] = j_left;
    ener_sf[sp] = accum;
  }

  // Landauer-Buttiker solutions for tight-binding central system
  MKL_INT wsamp = 10000;
  std::vector<double> w_lb = linspace(-1.0 * w, 1.0 * w, wsamp);
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

  // Transmission function
  std::vector<double> transmission(wsamp, 0.0);
  for(MKL_INT i = 0; i < wsamp; ++i){

    // M(w) = Iw - H_s - \Sigma
    std::vector< std::complex<double> > m_diag(L, 0.0);
    for(MKL_INT j = 0; j < L; ++j){
      m_diag[j] = w_lb[i] - epsilon_i;
    }
    std::complex<double> ig(0.0, Gamma / 2.0);
    m_diag[0] += ig;
    m_diag[L - 1] += ig;

    // Determinant of M(w)
    std::complex<double> f_m1(0.0, 0.0);
    std::complex<double> f_0(1.0, 0.0);
    std::complex<double> f_n;
    for(MKL_INT j = 0; j < L; ++j){
      f_n = (m_diag[j] * f_0) - (t_s * t_s * f_m1);
      f_m1 = f_0;
      f_0 = f_n;
    }

    if( L == 1 ) transmission[i] = (Gamma * Gamma) / std::norm(f_n);
    else transmission[i] = (Gamma * Gamma * std::pow(t_s, 2 * (L - 1))) / std::pow( std::abs(f_n), 2 );
  }

  std::vector<double> curr_lb(samp, 0.0);
  std::vector<double> ener_lb(samp, 0.0);
  for(MKL_INT i = 0; i < samp; ++i){
 
    // Density and current
    for(MKL_INT j = 0; j < wsamp - 1; ++j){
      curr_lb[i] += (1.0 / (2.0 * pi)) * transmission[j] * (f_l_lb[j] - f_r_lb[j]) * dw;
      ener_lb[i] += (1.0 / (2.0 * pi)) * w_lb[j] * transmission[j] * (f_l_lb[j] - f_r_lb[j]) * dw;
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
  }

  std::cout << "# L = " << L << " N = " << N << std::endl;
  std::cout << "# T_L = " << t_l << " T_R = " << t_r << std::endl;
  std::cout << "# mu_L = " << mu_l << " mu_R = " << mu_r << std::endl;
  std::cout << "# epsilon = " << epsilon_i << " t_S = " << t_s << std::endl;
  std::cout << "# Gamma = " << Gamma << std::endl;

  for(MKL_INT sp = 0; sp < samp; ++sp){
    std::cout << L << " " << N << " " << curr_lb[sp] << " " << curr_sf[sp] << " " << ener_lb[sp] << " " << ener_sf[sp] << std::endl;
    //std::cout << mu - epsilon_i << " " << V << " " << power_lb[sp] << " " << power_sf[sp] << " " << efficiency_lb[sp] << " " << efficiency_sf[sp] << std::endl;
  }

  return 0;
}
