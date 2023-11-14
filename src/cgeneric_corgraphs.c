
/* cgeneric_corgraphs.c
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

#include "cgeneric_defs.h"

double *inla_cgeneric_corgraphs(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
{

	double *ret = NULL;
	int i, j, k, N, M, ith, nth;

	// the size of the model
	assert(data->n_ints > 1);
	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

	assert(!strcasecmp(data->ints[1]->name, "debug"));     // this will always be the case
	int debug = data->ints[1]->ints[0];		       // this will always be the case
	assert(debug >= 0);				       // just to 'find an use for "debug" ...'
	if (debug>0) debug = 1;

	assert(!strcasecmp(data->ints[2]->name, "NC"));
	int NC = data->ints[2]->ints[0];
	assert(NC > 0);
	double sigmas[NC];
	int ifix[NC];
	M = NC + ( ( NC*(NC-1) ) / 2 ) ;

        assert(!strcasecmp(data->ints[3]->name, "NP"));
        int NP = data->ints[3]->ints[0];
        assert(NP >= 0);

	assert(!strcasecmp(data->doubles[0]->name, "lambda"));
	double lambda = data->doubles[0]->doubles[0];
	assert(lambda>0);

	assert(!strcasecmp(data->doubles[1]->name, "slambdas"));
        inla_cgeneric_vec_tp *slambdas = data->doubles[1];
	assert(slambdas->length>0);
	assert(slambdas->length>=NC);

        assert(!strcasecmp(data->doubles[2]->name, "plambdas"));
        inla_cgeneric_vec_tp *plambdas = data->doubles[2];
        assert(plambdas->length>0);
	assert(plambdas->length>=NC);

        nth = 0;
	for(i = 0; i<NC; i++) {
		if (iszero(plambdas->doubles[i])) {
			ifix[i] = 1;
		} else {
			ifix[i] = 0;
			nth++;
		}
	}

	switch (cmd) {
	case INLA_CGENERIC_GRAPH:
	{
		k = 0;
		ret = Calloc(2 + 2 * M, double);
		ret[0] = N;				       /* dimension */
		ret[1] = M;				       /* number of (i <= j) */
      for (int i = 0; i < N; i++) {
        for(int j = i; j < N; j++) {
          ret[k++] = i;
        }
      }
      for (int i = 0; i < N; i++) {
        for(int j = i; j < N; j++) {
          ret[k++] = j;
        }
      }
	}
		break;

	case INLA_CGENERIC_Q:
	{
		int offset = 2;
		ret = Calloc(offset + M, double);

		assert(!strcasecmp(data->mats[0]->name, "xx"));
		inla_cgeneric_mat_tp *xx = data->mats[0];
		assert(xx->nrow == N);
		assert(xx->ncol == N);

		int ipiv[N], info=0;

		dgesv_(&N, &N, &xx->x[0], &N, ipiv, &ret[offset], &N, &info, F_ONE);

		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;				       /* REQUIRED */
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
		ith = 0;
		ret[ith++] = (double) nth;
		for(i = 0; i < NC; i++) {
			if (ifix[i] == 0) {
				ret[ith++] = 0.0;
			}
		}
	}

		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		ret = Calloc(1, double);
		// PC-priors
		ret[0] = 0.0;
		ith = 0;
		double lam = 0;
		for(i = 0; i < nth; i++) {
			if (ifix[i] == 0) {
				lam = -log(plambdas->doubles[i]) / slambdas->doubles[i];
				ret[0] += log(lam) + theta[ith] - lam * exp(theta[ith]);
				ith++;
			}
		}
		assert(ith == nth);
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
