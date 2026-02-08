
/* cgeneric_pc_correl.c
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

double *inla_cgeneric_pc_correl(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{

	// This is the cgeneric implementatin of the
	// PC-prior for a correlation matrix C with dimension N,
	// given (scalar) parameter 'lambda'.
	// The correlation matrix is parametrized using the
	//  1. Canonical Partial Correlation, Lewandowski-Kurowicka-Joe (2009).
	//  2. hypershere decomposition, Rapisarda, Brigo and Mercurio (2007).
	// It returns for if 'cmd' is
	// 'graph': i,j index set for the upper triangle of Q;
	// 'Q': the inverse of C;
	// 'mu': 0.0 (zero);
	// 'initial': theta[k] = 3.0
	// 'log_prior': the PC-prior

	double *ret = NULL;
	int i, j, k, N, M;

	// the size of the model
	assert(data->n_ints > 1);
	assert(!strcasecmp(data->ints[0]->name, "n"));	    // this will always be the case
	N = data->ints[0]->ints[0];
	assert(N > 0);
	M = (int) ((double) N * ((double) (N + 1)) / 2.0);

	assert(!strcasecmp(data->ints[1]->name, "debug"));  // this will always be the case
//	int debug = data->ints[1]->ints[0];

	assert(!strcasecmp(data->doubles[1]->name, "sigmaref"));
	int np1 = data->doubles[1]->len;
	for (i = 0; i < np1; i++) {
	  assert(data->doubles[1]->doubles[i] > 0);
	}

	assert(!strcasecmp(data->ints[3]->name, "ifixed"));
	int npars = data->ints[3]->len;
	int np2 = npars-np1;
	double th0[npars];

	int nunk1=0, nunk2=0;
	for(i=0; i<np1; i++) {
	  nunk1 += (data->ints[3]->ints[i]==0);
	}
	for(i=0; i<np2; i++) {
	  nunk2 += (data->ints[3]->ints[np1+i]==0);
	}
	int nUnk=nunk1+nunk2;
//	printf("np1 %d, np2 %d, npars %d, nu1 %d, nu2 %d \n",
  //      np1, np2, npars, nunk1, nunk2);

	double actualtheta[npars], actualsigmas[N];
	if (theta) {

	  assert(!strcasecmp(data->ints[4]->name, "iparams"));
	  assert(!strcasecmp(data->doubles[2]->name, "sigmaprob"));
	  assert(!strcasecmp(data->doubles[3]->name, "lconst"));
	  assert(!strcasecmp(data->doubles[4]->name, "thetabase"));

	  k=0;
	  for(i=0; i<np1; i++) {
	    if(data->ints[3]->ints[i]) {
	      th0[i] = log(data->doubles[1]->doubles[i]);
	    } else {
	      th0[i] = theta[k++];
	    }
	  }
	  for(i=0; i<np2; i++) {
	    j = np1 + i;
	    if(data->ints[3]->ints[j]) {
	      th0[j] = data->doubles[4]->doubles[i];
	    } else {
	      th0[j] = theta[k++];
	    }
	  }
	  assert(nUnk==k);

	  for(i=0; i<npars; i++) {
	    actualtheta[i] = th0[data->ints[4]->ints[i]];
	  }

	  for(i=0; i<N; i++) {
	    actualsigmas[i] = exp(actualtheta[i]);
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

	}
		break;
	case INLA_CGENERIC_Q:
	{
		// Q = (CC)^{-1}, with C = LL'
		int offset = 2;

		ret = Calloc(offset + M, double);
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;				       /* REQUIRED */

// Cholesky of the correlation matrix
//    if(parametrization==1) { // cpc parametrization
      double ldet, aJac;
      assert(!strcasecmp(data->ints[2]->name, "iLtheta"));

      double ltheta[M-N];
      for(i=0; i<(M-N); i++) {
        ltheta[i] = 0.0;
      }
      k = N;
      for(i=0; i<data->ints[2]->len; i++) {
        ltheta[data->ints[2]->ints[i]] = actualtheta[k++];
      }
      cpcCholesky(&N, &ltheta[0], &ret[offset], &ldet, &aJac);
	//  } else { // old parametrization
    //  double hld;
      //theta2gamma2Lcorr(N, &hld, &theta[0], &ret[offset]);
    //}

/*
		if (debug > 999) {
			printf("L:\n");
			for (i = 0; i < N; i++) {
				k = i;
				for (j = 0; j <= i; j++) {
					printf("%2.3f ", ret[offset + k]);
					k += (N - j - 1);
				}
				printf("\n");
			}
		}
*/

  // include sigmas in the Cholesky
  k = 2;
  for(i=0; i<N; i++) {
    for(j=i; j<N; j++) {
      ret[k] *= actualsigmas[j];
      k++;
    }
  }

  // chol2inv
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
		ret = Calloc(nUnk + 1, double);
	  ret[0] = nUnk;
	  for (i = 0; i < nUnk; i++) {
	    ret[1 + i] = 0.0;
	  }

	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		ret = Calloc(1, double);
	  ret[0] = 0.0;

	  // PC prior for sigma[i]
	  if(nunk1>0) {
	    double lam;
	    k=0;
	    for (i = 0; i < np1; i++) {
	      if (data->ints[3]->ints[i]==0) {
	        lam = -log(data->doubles[2]->doubles[i]);
	        lam /= data->doubles[1]->doubles[i];
	        ret[0] += pclogsigma(theta[k++], lam);
	      }
	    }
	  }

	  // p(theta|lambda) = p(xi|lambda) |det(I(theta0))|
	  // lconst = |det(I)^{1/2}|
	  if(nunk2>0) {
	    assert(!strcasecmp(data->doubles[0]->name, "lambda"));
	    assert(!strcasecmp(data->doubles[3]->name, "lconst"));
	    assert(!strcasecmp(data->doubles[4]->name, "thetabase"));
	    assert(!strcasecmp(data->mats[0]->name, "Ihalf"));
	    double thb[nunk2];
	    k=0;
	    for(i=0; i<np2; i++) {
	      j = np1 + i;
	      if(data->ints[3]->ints[j]==0) {
	        thb[k++] = data->doubles[4]->doubles[i];
	      }
	    }
//	    printMat(thb,1,nunk2,"thb\n");
	    ret[0] += pcmultivar(
	      nunk2, data->doubles[0]->doubles[0], &thb[0],
        &data->mats[0]->x[0], &data->doubles[3]->doubles[0], &theta[nunk1]);
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
