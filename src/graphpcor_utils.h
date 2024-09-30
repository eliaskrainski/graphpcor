
/* graphpcor_utils.h
 *
 * Copyright (C) 2023 Elias T Krainski
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * The author's contact information:
 *
 *        Elias T Krainski
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 */


#include <stdio.h>

void exchangeableU(int n, double r, double *cc);
void dl2fullQ(int n, double *d, double *l, double *q);
void dl2Qu(int n, double *d, double *l, double *qu);
void cov2cor(int verbose, int n, double *cc) ;
double cov2kld(int verbose, int n, double *C0, double *C1);
void covariance_parent_children(int verbose, int np, int N, int niiv, int *iiv, int *jjv, int *ipar, int *itop, double *sch, double *v2, double *CC) ;
void correlation_parent_children(int verbose, int np, int N, int niiv, int *iiv, int *jjv, int *ipar, int *itop, double *sch, double *v2, double *CC) ;
void theta_parent_children_kldh(int verbose, int np, int N, int niiv, int *iiv, int *jjv, int *ipar, int *itop, double *sch, double *theta, double hs, double *kld, double *kldd) ;
