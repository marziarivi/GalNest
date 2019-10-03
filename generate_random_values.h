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

#ifndef _generate_values_h
#define _generate_values_h

#ifdef __cplusplus
extern "C" {
#endif
    
void evaluate_CDF(double *F, double min_value, double max_value, double (*CDFunc)(double,double), double param);
double generate_random_data(double u, double min_value, double max_value, double *F);

#ifdef __cplusplus
}
#endif

#endif
