#ifndef __UTILS_H
#define __UTILS_H

#include <complex>
#include <vector>
#include <iostream>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"

typedef std::vector< std::complex<double> > MZType;

namespace Utils
{
  MZType super_fermion_solve(std::vector<double> &H,
                             std::vector<double> &Ge,
                             std::vector<double> &Gd,
                             MKL_INT &n);
}
#endif
