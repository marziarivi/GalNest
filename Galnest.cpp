/*
 * Copyright (c) 2019 Marzia Rivi
 *
 * This file is part of GalNest.
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


//  GalNest.cpp
/*
    Detect and measure star forming galaxy parameters from a radio observation.
    A single model approach is adopted using MultNest.

    Data must be provided in a Measurement Set.

    Command line input parameters:
    argv[1]  filename Measurement Set
    argv[2]  number of live points
*/

#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <new>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "datatype.h"
#include "utils.h"
#include "measurement_set.h"
#include "likelihood.h"
#include "generate_random_values.h"
#include "distributions.h"
#include "multinest.h"


using namespace std;

int main(int argc, char *argv[])
{
    int nprocs, rank, num_threads=1;
#ifdef USE_MPI
    MPI_Init(&argc, &argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs) ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double start_tot = MPI_Wtime();
#else
    nprocs=1;
    rank=0;
    
    long long start_tot;
    start_tot = current_timestamp();
#endif
#ifdef _OPENMP
#pragma omp parallel
    num_threads = omp_get_num_threads();
    if (rank==0) cout << "Number of OpenMP threads = " << num_threads << endl;
#endif
    
    if (argc < 4)
    {
        cout << "ERROR: parameter missing!" << endl;
        cout << "usage: GalNest <filename MS> <num_live_points> <FoV [arcsec]> <resume flag 1/0>" << endl;
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD,rank);
#else
        exit(EXIT_FAILURE);
#endif
    }

#ifdef USE_MPI
    double data_time = 0.;
    double fitting_time = 0.;
    double start_data,end_data,start_fitting,end_fitting;
    start_data = MPI_Wtime();
#else
    double data_time = 0;
    double fitting_time = 0;
    long long start_data,end_data,start_fitting,end_fitting;
    start_data = current_timestamp();
#endif

    // Read Measurement Set --------------------------------------------------------------------------------------------------------------------------------------------------------------------
    char filename[50];
    sprintf(filename,"%s%d.ms", argv[1],7+rank);
    cout << "rank " << rank << ": " << filename << endl;
    RL_MeasurementSet* ms = ms_open(filename);

    //double RA = ms_phase_centre_ra_rad(ms);                 // Phase Centre coordinates
    //double Dec = ms_phase_centre_dec_rad(ms);   
    unsigned int num_stations = ms_num_stations(ms);        // Number of stations
    unsigned int num_channels = ms_num_channels(ms);        // Number of frequency channels
    unsigned int num_rows = ms_num_rows(ms);                // Number of rows 
    double freq_start_hz = ms_freq_start_hz(ms);            // Start Frequency, in Hz
    double channel_bandwidth_hz = ms_freq_inc_hz(ms);       // Frequency channel bandwidth, in Hz
    double full_bandwidth_hz = channel_bandwidth_hz * num_channels;  // Frequency total bandwidth, in Hz
    int time_acc = ms_time_inc_sec(ms);                     // accumulation time (sec)

    //double efficiency = 0.9;     // system efficiency
    //double SEFD = 320;    // System Equivalent Flux Density (in Jy) of JVLA antenna 
    double ref_frequency_hz = 5.5e+9;  //Reference frequency in Hz at which fluxes are measured
    
    unsigned int num_baselines = num_stations * (num_stations - 1) / 2;
    if (rank==0)
    {
        cout << "Number of channels: " << num_channels << endl;
        cout << "Channels bandwidth (Hz): " << channel_bandwidth_hz << endl;
        cout << "Reference frequency (Hz): " << ref_frequency_hz << endl;
        cout << "Accumulation time (sec): " << time_acc << endl;
    }
    cout << "rank " << rank << ": Starting frequency (Hz): " << freq_start_hz << endl;
    cout << "rank " << rank << ": Number of rows: " << num_rows << endl;
    
    double sizeGbytes, totGbytes = 0.;
    
    // Allocate and read uv coordinates 
    unsigned long int num_coords = ms_num_rows(ms);
    double* uu_metres = new double[num_coords];
    double* vv_metres = new double[num_coords];
    double* ww_metres = new double[num_coords];
    sizeGbytes = 3*num_coords*sizeof(double)/((double)(1024*1024*1024));
    cout << "rank " << rank << ": allocated original coordinates: " << sizeGbytes  << " GB" << endl;
    totGbytes += sizeGbytes;
    
    int status = 0;
    double len = ms_read_coords(ms,0,num_coords,uu_metres,vv_metres,ww_metres,&status);
    
    // Allocate and read noise variance
    float *sigma2;
    try
    {
      sigma2 = new float[num_coords];
      sizeGbytes = num_coords*sizeof(float)/((double)(1024*1024*1024));
      cout << "rank " << rank << ": allocated SIGMA column, size = " << sizeGbytes << " GB" << endl;
      totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
      cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << endl;
    }
    ms_read_sigma(ms,0,num_coords,sigma2,&status);

    // Allocate and read Data visibilities
    unsigned long int num_vis  = (unsigned long int) num_channels * num_coords;
    complexd *visData;
    try
    {
        visData = new complexd[num_vis];
        sizeGbytes = num_vis*sizeof(complexd)/((double)(1024*1024*1024));
        cout << "rank " << rank << ": allocated original data visibilities: " << num_vis << ", size = " << sizeGbytes  << " GB" << endl;
        totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
        cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << '\n';
    }

    ms_read_vis(ms, 0, 0, num_channels, num_rows, "DATA", visData, &status);

    bool *Flag;
    try
    {
	Flag = new bool[num_vis];
	sizeGbytes = num_vis*sizeof(bool)/((double)(1024*1024*1024));
	cout << "rank " << rank << ": allocated FLAG column, size = " << sizeGbytes << " GB" << endl;
	totGbytes += sizeGbytes;
    }
    catch (bad_alloc& ba)
    {
	cerr << "rank " << rank << ": bad_alloc caught: " << ba.what() << endl;
    }    
    ms_read_Flag(ms, 0, 0, num_channels, num_rows, "FLAG", Flag, &status);

    cout << "rank " << rank << ": Finish reading MS..." << endl;
    ms_close(ms);
    if (status)
    {
        cout << "rank " << rank << ": MS ERROR! " << status << endl;
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD,rank);
#else
        exit(EXIT_FAILURE);
#endif
    } 

#ifdef USE_MPI
    end_data = MPI_Wtime();
    data_time = end_data - start_data;
#else
    end_data = current_timestamp();
    data_time = (double)(end_data - start_data)/1000.;
#endif
    
    // Pre-compute wavenumber and spectral factor for each channel 
    // They corresponds to the central frequency of each channel
    double *wavenumbers = new double[num_channels];
    double ch_freq = freq_start_hz + 0.5*channel_bandwidth_hz;
    double *spec = new double[num_channels];

    for (unsigned int ch = 0; ch < num_channels; ch++)
    {
        wavenumbers[ch] = 2.0 * PI * ch_freq / C0;
        spec[ch] = pow(ch_freq/ref_frequency_hz,-0.7);
        ch_freq += channel_bandwidth_hz;
    }

    // Evaluate CDF of priors
    double *F_flux;
    double *F_scale;
    double *F_ellip;

    if (rank == 0)
    {
      F_flux = new double[N+1];
      F_scale = new double[N+1];
      F_ellip = new double[N+1]; 
      evaluate_CDF(F_flux, Smin, Smax, flux_CDF,beta);
      evaluate_CDF(F_scale, Rmin, Rmax, scale_CDF, mu);
      double h = Emax/N;
      F_ellip[0] = 0.;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
      for (unsigned int i=1; i<N; i++) F_ellip[i] = CDF(e_pdf,i*h);
    }

    // Set likelihood computation parameters
    likelihood_params par;
    par.FOV = atof(argv[3])*ARCS2RAD;
    par.nchannels = num_channels;
    par.band_factor = channel_bandwidth_hz*PI/C0;
    par.acc_time = time_acc;
    par.spec = spec;
    par.wavenumbers = wavenumbers; // wavenumbers for the model
    par.sigma2 = sigma2;
    //par.sigma2 = (SEFD*SEFD)/(2.*time_acc*channel_bandwidth_hz*efficiency*efficiency); // visibility noise variance
    par.ncoords = num_coords;
    par.uu = uu_metres;
    par.vv = vv_metres;
    par.ww = ww_metres;
    par.flag = Flag;
    par.data = visData;
    if (rank == 0)
    {
      par.F_flux = F_flux;
      par.F_scale = F_scale;
      par.F_ellip = F_ellip;
    //  cout << "sigma_vis = " << sqrt(par.sigma2) << " Jy" << endl;
    }
#ifdef USE_MPI
    start_fitting = MPI_Wtime();
#else
    start_fitting = current_timestamp();
#endif
    
    if (rank == 0)
    { 
      // Multinest Call  -----------------------------------------------------------------------------------------------------------------------------------
      // Set multinest sampling parameters 
      int IS = 0;				// do Nested Importance Sampling?   	
      int mmodal = 1;				// do mode separation?	
      int ceff = 0;				// run in constant efficiency mode?	
      int nlive = atoi(argv[2]);		// number of live points
      double efr = 0.8;				// set the required efficiency	
      double tol = 0.1;				// tol, defines the stopping criteria	
      int ndims = 6;				// dimensionality (no. of free parameters)
      int nPar = 6;				// total no. of parameters including free & derived parameters	
      int nClsPar = 2;				// no. of parameters to do mode separation on	
      int updInt = 100; 			// after how many iterations feedback is required & the output files should be updated
						// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
      double Ztol = -1E90;			// all the modes with logZ < Ztol are ignored
      int maxModes = 100;			// expected max no. of modes (used only for memory allocation)
      int pWrap[ndims];				// which parameters to have periodic boundary conditions?
      for(int i = 0; i < ndims; i++) pWrap[i] = 0;
      char root[100] = "/share/data1/mrivi/JVLA-5GHz/JVLA-2ch-500-5GHZ-";	// root for output files
      int seed = 1;				// random no. generator seed, if < 0 then take the seed from system clock
      int fb = 1;				// need feedback on standard output?
      int resume = atoi(argv[4]);		// resume from a previous job?
      int outfile = 1;				// write output files?
      int initMPI = 0;				// initialize MPI routines?, relevant only if compiling with MPI
                                                // set it to F if you want your main program to handle MPI initialization
      double logZero = -1E90;                   // points with loglike < logZero will be ignored by MultiNest
      int maxiter = -1;			// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it
                                                // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
      void *context = (void *) &par;		// not required by MultiNest, any additional information user wants to pass

      // calling MultiNest
      nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, 
		    pWrap, fb, resume, outfile, initMPI,logZero, maxiter, LogLike, dumper, context);	
    }
#ifdef USE_MPI
    else
    {
       while (1)
       {
         double Cube[6];
         MPI_Bcast(Cube, 6 , MPI_DOUBLE, 0, MPI_COMM_WORLD);

         //compute loglikelihood
         double loc_chisq = compute_chisq(par.nchannels, par.spec, par.wavenumbers, par.band_factor, par.acc_time, Cube[4], Cube[5], Cube[3], Cube[2], Cube[0], Cube[1], par.ncoords, par.uu, par.vv, par.ww, par.sigma2, par.flag, par.data);
         
         double chisq;
         MPI_Reduce(&loc_chisq, &chisq, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
        }
      }
#endif

#ifdef USE_MPI    
    end_fitting = MPI_Wtime();
    fitting_time = end_fitting - start_fitting;
    double end_tot = MPI_Wtime();
    double total_time = end_tot - start_tot;
#else
    end_fitting = current_timestamp();
    fitting_time = (double)(end_fitting - start_fitting)/1000.;
    long long end_tot = current_timestamp();
    double total_time = (double)(end_tot - start_tot)/1000.;
#endif

    cout << "Set up time (sec): " << total_time - data_time - fitting_time << endl;
    cout << "Data reading time (sec): " << data_time << endl;
    cout << "Data fitting computation time (min): " << fitting_time/60. << endl;
    cout << "Total time (min): " << total_time/60. << endl;
    
    // free memory ----------------------------------------------------------------------------------------------------------------
    delete[] visData;
    delete[] Flag;
    delete[] sigma2;
    delete[] uu_metres;
    delete[] vv_metres;
    delete[] ww_metres;
    delete[] wavenumbers;
    delete[] spec;
    delete[] F_flux;
    delete[] F_scale;
    delete[] F_ellip;

#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    return 0;
}
