#include "Utils.h"

namespace Utils
{
  MZType super_fermion_solve(std::vector<double> &H,
                             std::vector<double> &Ge,
                             std::vector<double> &Gd,
                             MKL_INT &n)
  {
    MZType block_1(n * n, 0.0);
    MZType block_2(n * n, 0.0);
    MZType block_3(n * n, 0.0);
    MZType block_4(n * n, 0.0);
    for(MKL_INT i = 0; i < n; ++i){
      for(MKL_INT j = 0; j < n; ++j){
        if(i == j){ 
          block_1[(i * n) + j] = std::complex<double>(H[(i * n) + j], (Ge[i] - Gd[i]) / 2.0);
          block_2[(i * n) + j] = std::complex<double>(0.0, Ge[i]);
          block_3[(i * n) + j] = std::complex<double>(0.0, Gd[i]);
          block_4[(i * n) + j] = std::complex<double>(H[(i * n) + j], (Gd[i] - Ge[i]) / 2.0);
        }
        else{
          block_1[(i * n) + j] = std::complex<double>(H[(i * n) + j], 0.0);
          block_4[(i * n) + j] = std::complex<double>(H[(i * n) + j], 0.0);
        }
      }
    }

    MZType H_super((2 * n) * (2 * n), 0.0);
    for(MKL_INT i = 0; i < n; ++i){
      for(MKL_INT j = 0; j < n; ++j){
        H_super[(i * (2 * n)) + j] = block_1[(i * n) + j];
        H_super[( (i * (2 * n)) + n ) + j] = block_2[(i * n) + j];
        H_super[(i * (2 * n)) + j + (2 * n * n)] = block_3[(i * n) + j];
        H_super[( (i * (2 * n)) + n ) + j + (2 * n * n)] = block_4[(i * n) + j];
      }
    } 
    
    MZType eigvals(2 * n, 0.0);
    MZType V((2 * n) * (2 * n), 0.0);
    MKL_INT info;

    // Diag
    info = LAPACKE_zgeev(LAPACK_ROW_MAJOR,
                         'N',
                         'V',
                         2 * n,
                         &H_super[0],
                         2 * n,
                         &eigvals[0],
                         NULL,
                         2 * n,
                         &V[0],
                         2 * n);
    if(info > 0){
      std::cout << "Eigenproblem" << std::endl;
      exit(1);
    }

    // Inverse
    MZType Vinv((2 * n) * (2 * n), 0.0);
    std::vector<MKL_INT> ipiv( 2 * n );
    cblas_zcopy( (2 * n) * (2 * n),
                 &V[0],
                 1,
                 &Vinv[0],
                 1);
    info = LAPACKE_zgetrf( LAPACK_ROW_MAJOR,
                           2 * n,
                           2 * n,
                           &Vinv[0],
                           2 * n,
                           &ipiv[0]);
    if(info > 0){
      std::cout << "LUproblem" << std::endl;
      exit(1);
    }
    info = LAPACKE_zgetri( LAPACK_ROW_MAJOR,
                           2 * n,
                           &Vinv[0],
                           2 * n,
                           &ipiv[0]);
    if(info > 0){
      std::cout << "INVproblem" << std::endl;
      exit(1);
    }

    // Matrix of normal mode correlations for the NESS
    MZType D((2 * n) * (2 * n), 0.0);
    for(MKL_INT i = 0; i < (2 * n); ++i){
      if(imag(eigvals[i]) > 0)
        D[(i * 2 * n) + i] = std::complex<double>(1.0, 0.0);
      else
        D[(i * 2 * n) + i] = std::complex<double>(0.0, 0.0);
    }

    // Compute G-> Important: This is already computing the transpose
    MZType temp((2 * n) * (2 * n), 0.0);
    MZType G((2 * n) * (2 * n), 0.0);
    double alpha = 1.0;
    double beta = 0.0;
    cblas_zgemm( CblasRowMajor,
                 CblasNoTrans,
                 CblasTrans,
                 2 * n,
                 2 * n,
                 2 * n,
                 &alpha,
                 &D[0],
                 2 * n,
                 &V[0],
                 2 * n,
                 &beta,
                 &temp[0],
                 2 * n);
    cblas_zgemm( CblasRowMajor,
                 CblasTrans,
                 CblasNoTrans,
                 2 * n,
                 2 * n,
                 2 * n,
                 &alpha,
                 &Vinv[0],
                 2 * n,
                 &temp[0],
                 2 * n,
                 &beta,
                 &G[0],
                 2 * n);
    // Get info about the system only
    MZType corr( n * n );
    for(MKL_INT i = 0; i < n; ++i)
      for(MKL_INT j = 0; j < n; ++j)
        corr[(i * n) + j] = G[(i * 2 * n) + j];

    return corr;
  }
}
