
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
	int i, j, k, N, M, n2, np, N2;

	// the size of the model
	assert(data->n_ints > 1);
	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	N = data->ints[0]->ints[0];
	assert(N > 0);
	M = N * N;

	assert(!strcasecmp(data->ints[1]->name, "p"));     // this will always be the case
	np = data->ints[1]->ints[0];
	assert(np > 0);
	n2 = N + np;
	N2 = (n2 * n2);

	assert(!strcasecmp(data->ints[2]->name, "i2th"));     // this will always be the case
	inla_cgeneric_vec_tp *i2th = data->ints[2];
	assert(i2th->len == np);

	assert(!strcasecmp(data->ints[3]->name, "iq2th"));     // this will always be the case
	inla_cgeneric_vec_tp *iq2th = data->ints[3];
	assert(i2th->len == iq2th->len);

	assert(!strcasecmp(data->ints[4]->name, "i1th"));     // this will always be the case
	inla_cgeneric_vec_tp *i1th = data->ints[4];

	assert(!strcasecmp(data->ints[5]->name, "iq1th"));     // this will always be the case
	inla_cgeneric_vec_tp *iq1th = data->ints[5];
	assert(i1th->len == iq1th->len);

	assert(!strcasecmp(data->doubles[0]->name, "sth"));     // this will always be the case
	inla_cgeneric_vec_tp *sth = data->doubles[0];
	assert(i1th->len == sth->len);

	assert(!strcasecmp(data->doubles[1]->name, "q"));
	inla_cgeneric_vec_tp *q = data->doubles[1];
	assert(q->len == M);

	assert(!strcasecmp(data->doubles[2]->name, "lambda"));
	double lambda = data->doubles[1]->doubles[2];
	assert(lambda>0);

	assert(!strcasecmp(data->doubles[3]->name, "slambdas"));
	inla_cgeneric_vec_tp *slambdas = data->doubles[3];
	assert(slambdas->len>0);
	assert(slambdas->len == N);

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

		int ipiv[N], info=0;

		double vv[N2], qq[N2], vc[M];
		for(i=0; i<M; i++) {
		  qq[i] = q->doubles[i];
		}
		for(i=0; i<N; i++) {
		  qq[iq2th->ints[i]] += exp(-2 * theta[i2th->ints[i]]);
		}
		k = 0;
		for(i=0; i < i1th->len; i++) {
		  qq[iq1th->ints[i]] = -sth->doubles[i] * exp( -2 * theta[i1th->ints[i]]);
		}

		dgesv_(&n2, &n2, qq, &n2, ipiv, vv, &n2, &info, F_ONE);

		k=0;
		int k2 = 0;
		for(j=0; j<N; j++) {
		  for(i=0; i<n2; i++) {
		    if((i<N) & (j<N)) {
		      vc[k2] = vv[k];
		      k2++;
		    }
		    k++;
		  }
		}

		dgesv_(&N, &N, vc, &N, ipiv, &ret[offset], &N, &info, F_ONE);

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
		ret = Calloc(n2 + 1, double);
		for(i = 0; i < N2; i++) {
			ret[1+i] = 1.0;
		}
	}

		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		ret = Calloc(1, double);
		// PC-priors
		ret[0] = 0.0;
		/*
		double lam = 0;
		for(i = 0; i < N; i++) {
			lam = slambdas->doubles[i];
		  ret[0] += log(lam) + theta[i] - lam * exp(theta[i]);
		} */
		for(i = 0; i < np; i++) {
		  ret[0] += -0.5 * SQR(theta[i]);
		}
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
