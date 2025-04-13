
/* cgeneric_pc_prec_correl.c
 *
 * Copyright (C) 2025 Elias Krainski
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

#include "graphpcor.h"
#include "graphpcor_utils.h"

double *inla_cgeneric_pc_prec_correl(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{

	// This is the cgeneric implementatin of the
	// PC-prior for a Precision matrix Q with dimension N,
	// given (scalar) parameter 'lambda'.
	// Q is parametrized as follows
	// step 1: | 1 |
	// | \theta[1] 1 SYMMETRIC |
	// Q0 = | \theta[2] \theta[n] |
	// | : ...  ...  |
	// | \theta[n-1] ...  \theta[m] 1 |
	// step 2: V = Q0^{-1}
	// step 3: S = diag(V)^{-1/2}
	// step 4: C = S %*% V %*% S : this is C(Q)
	// step 5: Q = C^{-1}
	// p(Q|lambda) = p_d( d(\theta[1:m]) | lambda) |J(C(Q))|
	// where p_d() is the prior for the distance.
	// d(\theta) = (\theta - \theta_0)' I(\theta_0) (\theta - \theta_0)
	// I(\theta_0) is the Fisher information at \theta_0
	// see section 6.2 of Simpson et. al. (2017)
	// It returns for if 'cmd' is
	// 'graph': i,j index set for the upper triangle of Q;
	// 'Q': the inverse of C;
	// 'mu': 0.0 (zero);
	// 'initial': theta[k] = 3.0
	// 'log_prior': the PC-prior

	double *ret = NULL;
	int i, j, k, N, M, nth;

	// the size of the model
	assert(data->n_ints > 1);
	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	N = data->ints[0]->ints[0];
	assert(N > 0);
	M = (int) ((double) N * ((double) (N + 1)) / 2.0);
	nth = (int) ((double) N * ((double) (N - 1)) / 2.0);

	assert(!strcasecmp(data->ints[1]->name, "debug"));     // this will always be the case
//	int debug = data->ints[1]->ints[0];

//  printf("debug = %d\n", debug);

	assert(!strcasecmp(data->doubles[0]->name, "lambda"));
	double lambda = data->doubles[0]->doubles[0];
	assert(lambda > 0);

	assert(!strcasecmp(data->doubles[1]->name, "lconst"));
	double lconst = data->doubles[1]->doubles[0];
/*
	if (debug > 999) {
		printf("N=%d, nth=%d, M=%d, lambda=%f, lc=%f\n", N, nth, M, lambda, lconst);
	}
*/
	double daux;
	double param[nth];

	if (theta) {
		param[nth - 1] = atan2(theta[nth - 1], theta[nth - 2]);
		if (param[nth - 1] < 0) {
			param[nth - 1] += 2.0 * M_PI;
		}
		daux = SQR(theta[nth - 1]) + SQR(theta[nth - 2]);
		if (nth > 2) {
			for (i = nth - 2; i > 0; i--) {
				// phi[i]
				param[i] = atan2(sqrt(daux), theta[i - 1]);
				daux += SQR(theta[i - 1]);
			}
		}
		// r = sqrt(sum_i(theta[i]^2))
		param[0] = sqrt(daux);
/*
		if (debug > 999) {
			printMat(param, 1, nth, "param:\n");
		}
*/
	} else {
		for (i = 0; i < nth; i++) {
			param[i] = NAN;
		}
	}

	switch (cmd) {
	case INLA_CGENERIC_GRAPH:
	{
		k = 2;
		ret = Calloc(k + 2 * M, double);
		ret[0] = N;				       /* dimension */
		ret[1] = M;				       /* number of (i <= j) */

		for (i = 0; i < N; i++) {
			for (j = i; j < N; j++) {
				ret[M + k] = j;
				ret[k++] = i;
			}
		}
/*
		if (debug > 999) {
			printf("ii: ");
			for (k = 0; k < M; k++) {
				printf("%f ", ret[2 + k]);
			}
			printf("\njj: ");
			for (k = 0; k < M; k++) {
				printf("%f ", ret[2 + M + k]);
			}
			printf("\n");
		}
*/

	}
		break;
	case INLA_CGENERIC_Q:
	{
		// Q = (CC)^{-1}
		// with C = LL' parametrized as in Section 6.2 of the PC-prior paper
		int offset = 2;

		ret = Calloc(offset + M, double);
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;				       /* REQUIRED */

		// compute the (upper) correlation matrix
		int std = 1;
		theta2Qcorrel(N, std, &theta[0], &ret[offset]);

		int info;
		char uplo = 'L';
		dpptrf_(&uplo, &N, &ret[offset], &info, F_ONE);	// chol
		dpptri_(&uplo, &N, &ret[offset], &info, F_ONE);	// chol2inv

	}
		break;
	case INLA_CGENERIC_MU:
	{
		// return (N, mu). if N==0 then mu is not needed as its taken to be mu[]==0
		ret = Calloc(1, double);
		ret[0] = 0;
	}
		break;

	case INLA_CGENERIC_INITIAL:
	{
		// return c(P, initials)
		// where P is the number of hyperparameters
		ret = Calloc(nth + 1, double);
		ret[0] = nth;
		for (i = 0; i < nth; i++) {
			ret[1 + i] = 3.0;
		}
	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		ret = Calloc(1, double);
		double ldJacobian;

		ldJacobian = ((double) (nth - 1)) * log(param[0]);
		if (nth > 2) {
			for (i = 1; i < (nth - 1); i++) {
				ldJacobian += ((double) (nth - 1 - i)) * log(sin(param[i]));
			}
		}
/*
		if (debug > 999) {
			printf("log det Jacobian = %2.7f\n", ldJacobian);
		}
*/
		// the log prior:
		// lconst should be equal to
		// log(lambda) -(m-1)*log(pi)-log(2)-log(|H|)
		ret[0] = lconst - lambda * param[0] + ldJacobian;

	}
		break;

	case INLA_CGENERIC_VOID:
	case INLA_CGENERIC_LOG_NORM_CONST:
	case INLA_CGENERIC_QUIT:
	default:
		break;
	}

	return (ret);
}
