/*
 * Copyright (c) 2018 Marzia Rivi
 *
 * This file is part of RadioLensfit.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see <http://www.gnu.org/licenses/>.
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "distributions.h"
#include "generate_random_values.h"

#ifdef __cplusplus
extern "C" {
#endif
   
void evaluate_CDF(double *F, double min_value, double max_value, double (*CDFunc)(double,double), double param)
{
   double h = (max_value-min_value)/N;
   double CFmin = (*CDFunc)(param,min_value);
   double CFrange = (*CDFunc)(param,max_value) - CFmin;

   F[0]=0.;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
   for (unsigned int i=1; i<=N; i++) F[i] = ((*CDFunc)(param,min_value+i*h) - CFmin)/CFrange;
}	


// Generate random data in the range [min_value,max_value] according to a distribution with a given Cumulative Distribution Function CDFunc(param, x)
double generate_random_data(double u, double min_value, double max_value, double *F)
{
  double h = (max_value-min_value)/N;
        
  unsigned int k=0;
  while ((k < N) && (u > F[k])) k++;
  // u is in the interval (F(k-1),F(k)), find F^{-1}(u) using a linear interpolation
  double y = min_value+h*(k-1)+h*(u-F[k-1])/(F[k]-F[k-1]);
  return y;

}
 

#ifdef __cplusplus
}
#endif

 
