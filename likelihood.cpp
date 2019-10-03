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
 *                  */


//  Computation of the likelihood function through galaxy model-fitting
//  It depends on 6 parameters: ellipcticy (e1,e2), scalelength (r), position (l,m), flux (S).

#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datatype.h"
#include "distributions.h"
#include "generate_random_values.h"
#include "likelihood.h"
#include "multinest.h"


// parameters: Cube = [l,m,flux,scalelength,e1,e2]
// return list of 6 values drawn randomly from each prior distribution
void LogLike(double* Cube, int &ndim, int &npars, double &lnew, void *context)
{
   likelihood_params *param = (likelihood_params *) context;

   // Position prior
   // Uniformly distributed positions between [-FOV/2 , FOV/2]
   Cube[0] = (Cube[0] - 0.5) * param->FOV;
   Cube[1] = (Cube[1] - 0.5) * param->FOV;

   // Flux prior 
   Cube[2] = generate_random_data(Cube[2], Smin, Smax, param->F_flux)*1e-6; 

   // Scalelength prior
   Cube[3] = generate_random_data(Cube[3], Rmin, Rmax, param->F_scale);

   // Ellipticity prior
   double e_mod = generate_random_data(Cube[4], 0., Emax, param->F_ellip);
   double theta = 2 * PI * Cube[5];
   Cube[4] = e_mod * cos(theta);
   Cube[5] = e_mod * sin(theta);
 
#ifdef USE_MPI
   MPI_Bcast(Cube, 6 , MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
   //compute loglikelihood
   double loc_chisq = compute_chisq(param->nchannels, param->spec, param->wavenumbers, param->band_factor, param->acc_time, Cube[4], Cube[5], Cube[3], Cube[2], Cube[0], Cube[1], param->ncoords, param->uu, param->vv, param->ww, param->sigma2, param->flag, param->data);

   double chisq;
#ifdef USE_MPI
   MPI_Reduce(&loc_chisq, &chisq, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#else
   chisq = loc_chisq;
#endif

   lnew = -chisq*0.5;
}


// Compute galaxy visibilities (sersic index 1) and chi-squared
double compute_chisq(unsigned int nchannels, double *spec, double *wavenumbers, double band_factor, double acc_time,
		     double e1, double e2, double scalelength, double flux, double l, double m, 
		     unsigned long int num_coords, double* uu_metres, double* vv_metres, double* ww_metres, float* sigma2, bool* flag, complexd* visData)
{
   double detA = 1.-e1*e1-e2*e2;
   double scale = scalelength*ARCS2RAD;  // scale in rad
   double scale_factor = (scale*scale)/(detA*detA);
   double n = sqrt(1.-l*l-m*m) - 1.;
   
   double chisq = 0.;

#ifdef _OPENMP
#pragma omp parallel for
#endif
   for (unsigned long int ch=0; ch < nchannels; ch++)
   {
     double spectra = spec[ch];
     double wavenumber = wavenumbers[ch];
     double wavenumber2 = wavenumber*wavenumber;

     unsigned long int ch_vis = ch*num_coords;
     double chi_par = 0.;

     for (unsigned long int i = 0; i < num_coords; ++i)
     {
       if (!flag[ch_vis])
       {
         double u = uu_metres[i];
         double v = vv_metres[i];
         double w = ww_metres[i];
                                                                                                                                   
         double phase = u*l+v*m+w*n;
         double smear = 1.;
         if (phase !=0.)
         {
           smear = band_factor*phase;
           smear = sin(smear)/smear;
         }
         phase = wavenumber*phase;
                                                             
         double k1 = (1.+e1)*u + e2*v;
         double k2 = e2*u + (1.-e1)*v;
         double den = 1. + scale_factor*wavenumber2*(k1*k1+k2*k2);
         double shape = spectra*flux/(den*sqrt(den));  //primary beam effect already included in the flux value
            
         double visGal_real = shape*cos(phase)*smear;
         double visGal_imag = shape*sin(phase)*smear;

         double temp = (visData[ch_vis].real - visGal_real)*(visData[ch_vis].real - visGal_real);
         temp += (visData[ch_vis].imag - visGal_imag)*(visData[ch_vis].imag - visGal_imag);
         chi_par += temp/sigma2[i];
       }
       ch_vis++;
     }
#ifdef _OPENMP
#pragma omp critical
#endif
     chisq = chisq + chi_par;
   }

   return chisq;
}



/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
// Arguments:
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value from the default (non-INS) mode
// INSlogZ						= log evidence value from the INS mode
// logZerr						= error on log evidence value
// context						void pointer, any additional information


void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr,
	       	double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context)
{
   // convert the 2D Fortran arrays to C++ arrays
    		
   // the posterior distribution
   // postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
    					
   int i, j;
    							
   double postdist[nSamples][nPar + 2];
   for( i = 0; i < nPar + 2; i++ )
   	for( j = 0; j < nSamples; j++ )
           postdist[j][i] = posterior[0][i * nSamples + j];

   // last set of live points
   // pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
   double pLivePts[nlive][nPar + 1];
   for( i = 0; i < nPar + 1; i++ )
   	for( j = 0; j < nlive; j++ )
   	  pLivePts[j][i] = physLive[0][i * nlive + j];
}

