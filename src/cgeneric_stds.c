
/* cgeneric_stds.c
 *
 * Copyright (C) 2026 Elias Krainski
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

double *inla_cgeneric_stds_dev(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{
  return &inla_cgeneric_stds(cmd, &theta[0], &data[0])[0];
}

double *inla_cgeneric_stds(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{

	// theta : vector of unknown parameters
	// 0 <= length(theta) <= n
	// actualtheta : vector of n model parameters
	// sigmas[i] = exp(actualtheta[i])

	double *ret = NULL;
	int i,  k;

	// the size of the model
	assert(data->n_ints > 1);
	assert(!strcasecmp(data->ints[0]->name, "n"));    // this will always be the case
	int N = data->ints[0]->ints[0];
	assert(N > 0);

//	printf("N = %d \n", N);

	assert(!strcasecmp(data->ints[1]->name, "debug")); // this will always be the case
//	int debug = data->ints[1]->ints[0];
	//assert(debug>=0);

	assert(!strcasecmp(data->ints[2]->name, "iparams"));
	assert(N == data->ints[2]->len);

	assert(!strcasecmp(data->ints[3]->name, "ifixed"));
	int npars = data->ints[3]->len;
	assert(npars>0);
	assert(npars<=N);
	double th0[npars];

//	printf("npars = %d \n", npars);

	int nUnk = npars;
	for(i=0; i<npars; i++) {
	  nUnk -= data->ints[3]->ints[i];
	}
	assert(nUnk<=npars);

//	printf("N=%d, npars=%d, nUnk=%d\n", N, npars, nUnk);

	assert(!strcasecmp(data->doubles[0]->name, "sigmaref"));
	assert(data->doubles[0]->len == npars);
	for (i=0; i<npars; i++) {
	  assert(data->doubles[0]->doubles[i] > 0);
	}

	assert(!strcasecmp(data->doubles[1]->name, "sigmaprob"));
	assert(npars == data->doubles[1]->len);

	double actualtheta[N];

	  k=0;
	  for(i=0; i<npars; i++) {
	    if(data->ints[3]->ints[i]) {
	      th0[i] = log(data->doubles[0]->doubles[i]);
	    } else {
	      if(theta) {
	        th0[i] = theta[k++]; // it can have no theta
	      } else {
	        th0[i] = NAN;
	      }
	    }
	  }
	  assert(nUnk==k);

	  for(i=0; i<N; i++) {
	    actualtheta[i] = th0[data->ints[2]->ints[i]];
	  }

	switch (cmd) {
	case INLA_CGENERIC_GRAPH:
	{
	  k = 2;
	  ret = Calloc(k + 2*N, double);
		ret[0] = N;   /* dimension */
		ret[1] = N;   /* (nnzaq/2 */
		for (i = 0; i < N; i++) {
			ret[k] = i;
			ret[N + k] = i;
			k++;
		}

	}
		break;
	case INLA_CGENERIC_Q:
	{
		k = 2;
		ret = Calloc(k + N, double);
		// memset(ret + k, 0, M * sizeof(double));
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = N;				       /* REQUIRED */

		for (i = 0; i < N; i++) {
			ret[k + i] = exp(-actualtheta[i] * 2);
		}

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
		ret = Calloc(nUnk + 1, double);
		ret[0] = nUnk;
		for (i = 1; i <= nUnk; i++) {
			ret[i] = -2.0;
		}
	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		ret = Calloc(1, double);
	  ret[0] = 0.0;

	  // PC prior for UNKNOWN sigma[i]
	  if(nUnk>0) {
	    double lam;
	    k=0;
	    for (i = 0; i < npars; i++) {
	      if (data->ints[3]->ints[i]==0) {
	        lam = -log(data->doubles[1]->doubles[i]);
	        lam /= data->doubles[0]->doubles[i];
	        ret[0] += gpc_pc_logsigma(theta[k++], lam);
	      }
	    }
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
