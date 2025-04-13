
/* cgeneric_LKJ.c
 *
 * Copyright (C) 2024 Elias Krainski
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

double *inla_cgeneric_LKJ(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{

	// Lewandowski-Kurowicka-Joe - LKJ prior
	// This is the cgeneric implementatin of the
	// LKJ prior for a correlation matrix C with dimension N,
	// given (scalar) parameter 'eta'.
	// Parametrized from a hypershere decomposition,
	// see Rapisarda, Brigo and Mercurio (2007).
	// It returns for if 'cmd' is
	// 'graph': i,j index set for the upper triangle of Q;
	// 'Q': the inverse of C;
	// 'mu': 0.0 (zero);
	// 'initial': theta[k] = 3.0
	// 'log_prior': the LKJ prior for c_k, k = 1,...,M

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

	assert(!strcasecmp(data->doubles[0]->name, "eta"));
	double eta = data->doubles[0]->doubles[0];
	assert(eta > 0);

	assert(!strcasecmp(data->doubles[1]->name, "lc"));
	double lc = data->doubles[1]->doubles[0];

/*
	if (debug > 999) {
		printf("N=%d, nth=%d, M=%d, eta=%f, lc=%f\n", N, nth, M, eta, lc);
	}
*/

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

	}
		break;
	case INLA_CGENERIC_Q:
	{
		// Q = (CC)^{-1}
		int offset = 2;

		ret = Calloc(offset + M, double);
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;				       /* REQUIRED */

		double hld;
		// build L that factorize C giving C = LL'
		theta2gamma2Lcorr(N, &hld, &theta[0], &ret[2]);

/*
		if (debug > 999) {
			printf("L:\n");
			k = 0;
			for (i = 0; i < N; i++) {
				for (j = i; j < N; j++) {
					printf("%2.3f ", ret[offset + k++]);
				}
				printf("\n");
			}
		}
*/
		// chol2inv: Q = C^{-1} computed from L
		int info;
		char uplo = 'L';
		dpptri_(&uplo, &N, &ret[offset], &info, F_ONE);

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

		// log determinant is computed in the theta2gamma2Lcorrel
		double lhdetC, ll[N * (N + 1) / 2], cc1[N * (N - 1) / 2], cc2[N * (N - 1) / 2];
		theta2gamma2Lcorr(N, &lhdetC, &theta[0], &ll[0]);

/*
		if (debug > 999) {
			printf("L[upper]:\n");
			k = 0;
			for (i = 0; i < N; i++) {
				for (j = i; j < N; j++) {
					printf("%2.3f ", ll[k++]);
				}
				printf("\n");
			}
			printL(ll, N, N, "Lcorr\n");
			printf("|R| = %2.5f, lc+(eta-1)|R| = %2.5f\n", 2.0 * lhdetC, lc + 2.0 * (eta - 1) * lhdetC);
		}
*/

		L2Cupper(N, &ll[0], &cc1[0]);


/*
 		if (debug > 999) {
			printf("C:\n");
			k = 0;
			for (i = 0; i < (N - 1); i++) {
				for (j = i; j < (N - 1); j++) {
					printf("%2.5f ", cc1[k++]);
				}
				printf("\n");
			}
		}
*/

		// compute the Jacobian[nth*nth]
		int kj, k2 = 0;
		double ldtmp, daux, h = 0.0005;
		double h2 = 2.0 * h;
		double mJacobian[nth * nth];
		for (kj = 0; kj < nth; kj++) {
			daux = theta[kj];
			theta[kj] = daux - h;
			theta2gamma2Lcorr(N, &ldtmp, &theta[0], &ll[0]);
			L2Cupper(N, &ll[0], &cc1[0]);
			theta[kj] = daux + h;
			theta2gamma2Lcorr(N, &ldtmp, &theta[0], &ll[0]);
			L2Cupper(N, &ll[0], &cc2[0]);
			theta[kj] = daux;
			// store derivatives
			for (k = 0; k < nth; k++) {
				daux = (cc2[k] - cc1[k]) / h2;
				mJacobian[k2++] = daux;
			}
		}

		int info, pivot[nth], lwork = 2 * nth + nth + 1;
		double work[lwork], tau[nth];
		dgeqp3_(&nth, &nth, &mJacobian[0], &nth, &pivot[0], &tau[0], &work[0], &lwork, &info, F_ONE);
		double ldJacobian = 0.0;
		for (i = 0; i < nth; i++) {
			ldJacobian += log(fabs(mJacobian[i * nth + i]));
		}
/*
		if (debug > 99) {
			printf("\ndet Jacobian = %f\n", ldJacobian);
		}
*/
		if (ldJacobian < 0) {
			ldJacobian *= -1.0;
		}
		// store the log-prior
		// 'lc' is a pre-computed constant
		ret[0] = lc + 2.0 * lhdetC * (eta - 1.0);
		ret[0] += ldJacobian;

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
