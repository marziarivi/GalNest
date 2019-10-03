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

#ifndef ____likelihood__
#define ____likelihood__

void LogLike(double* Cube, int &ndim, int &npars, double &lnew, void *context);

double compute_chisq(unsigned int nchannels, double *spec, double *wavenumbers, double band_factor, 
		double acc_time, double e1, double e2, double scalelength, double flux, double l, double m, 
		unsigned long int num_coords, double* uu_metres, double* vv_metres, double* ww_metres, 
		float* sigma2, bool* flag, complexd* visData);
void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr,
	       	double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context);


#endif /* defined(____likelihood__) */
