/*
 * Copyright (c) 2018 Marzia Rivi
 *
 * This file is part of RadioLensfit2.
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

#ifndef DATATYPE_H_
#define DATATYPE_H_

/**
 * @file datatype.h
 */

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef C0
#define C0 299792458.0
#endif

#ifndef DEG2RAD
#define DEG2RAD PI/180.0
#endif

#ifndef ARCS2RAD
#define ARCS2RAD PI/648000.0
#endif


typedef struct{
    double real;
    double imag;
} complexd;

typedef struct{
    unsigned int ncoords;
    unsigned int nchannels;
    double band_factor;
    double acc_time;
    double *spec;
    double *wavenumbers;
    float  *sigma2;
    double *uu;
    double *vv;
    double *ww;
    bool   *flag;
    complexd *data;
    double *F_flux;
    double *F_scale;
    double *F_ellip;
    double FOV;
} likelihood_params;


#endif
