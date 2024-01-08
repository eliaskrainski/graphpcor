
/* corgraphs_utils.h
 *
 * Copyright (C) 2023 Elias Krainski
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
 *        Elias Krainski
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#if __GNUC__ > 7
typedef size_t fortran_charlen_t;
#else
typedef int fortran_charlen_t;
#endif
#define F_ONE ((fortran_charlen_t)1)

void dgemm_(char *transa, char *transb,
            int *m, int *n, int *k, double *alpha,
            double *a, int *lda, double *b, int *ldb,
            double *beta, double *c, int *ldc, fortran_charlen_t);

void dgesv_(int *N, int *NRHS, double *A, int *LDA, int *IPIV,
            double *B, int *LDB, int *INLFO, fortran_charlen_t);

void dpotrf_(char *uplo, int *N, double *A, int *LDA,
             int *INLFO, fortran_charlen_t);

void dposv_(char *uplo, int *N, int *NRHS, double *A, int *LDA,
            double *B, int *LDB, int *INLFO, fortran_charlen_t);

void cov2cor(int verbose, int n, double *cc) ;
double cov2kld(int verbose, int n, double *C0, double *C1);
void covariance_parent_children(int verbose, int np, int N, int niiv, int *iiv, int *jjv, int *ipar, int *itop, double *sch, double *v2, double *CC) ;
void correlation_parent_children(int verbose, int np, int N, int niiv, int *iiv, int *jjv, int *ipar, int *itop, double *sch, double *v2, double *CC) ;
void theta_parent_children_kldh(int verbose, int np, int N, int niiv, int *iiv, int *jjv, int *ipar, int *itop, double *sch, double *theta, double hs, double *kld, double *kldd) ;

