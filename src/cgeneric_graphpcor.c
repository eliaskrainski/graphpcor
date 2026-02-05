
/* cgeneric_graphpcor.c
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

double *inla_cgeneric_graphpcor(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{

	// theta : vector of unknown parameters
	// 0 <= length(theta) <= n + m
	// actualtheta : vector of n+m model parameters
	// actualtheta = { log(sigmas), lowtheta }
	// sigmas[i] = exp(actualtheta[i])
	// actualtheta[n+1:m] = H^{1/2}(actualtheta[n+1:m] - base[1:m])
	// = H^{1/2} actualtheta[n+1:m] - thetabasescaled

	double *ret = NULL;
	int i, j, k;

	// the size of the model
	assert(data->n_ints > 1);
	assert(!strcasecmp(data->ints[0]->name, "n"));    // this will always be the case
	int N = data->ints[0]->ints[0];
	assert(N > 0);

	assert(!strcasecmp(data->ints[1]->name, "debug")); // this will always be the case
//	int debug = data->ints[1]->ints[0];
	//assert(debug>=0);

	assert(!strcasecmp(data->ints[2]->name, "ne"));
	int ne = data->ints[2]->ints[0];
	assert(ne > 0);
	int M = N + ne;

	assert(!strcasecmp(data->ints[3]->name, "nfi"));
	int nfi = data->ints[3]->ints[0];

	assert(!strcasecmp(data->ints[4]->name, "ii"));
	inla_cgeneric_vec_tp *ii = data->ints[4];
	assert(M == ii->len);

	assert(!strcasecmp(data->ints[5]->name, "jj"));
	inla_cgeneric_vec_tp *jj = data->ints[5];
	assert(M == jj->len);

	assert(!strcasecmp(data->ints[6]->name, "iuq"));
	inla_cgeneric_vec_tp *iuq = data->ints[6];
	assert(M == iuq->len);

	assert(!strcasecmp(data->ints[7]->name, "iuqpac"));
	inla_cgeneric_vec_tp *iuqpac = data->ints[7];
	assert(M == iuqpac->len);

	assert(!strcasecmp(data->ints[8]->name, "ifi"));
	inla_cgeneric_vec_tp *ifi = data->ints[8];
	assert(nfi == ifi->len);

	assert(!strcasecmp(data->ints[9]->name, "jfi"));
	inla_cgeneric_vec_tp *jfi = data->ints[9];
	assert(nfi == jfi->len);

	assert(!strcasecmp(data->doubles[1]->name, "sigmaref"));
	int np1 = data->doubles[1]->len;

	assert(!strcasecmp(data->ints[11]->name, "itheta"));
	int npars = data->ints[11]->len;
	int np2 = npars-np1;

  assert(!strcasecmp(data->ints[10]->name, "ifixed"));
  int nunk1=0, nunk2=0;
  for(i=0; i<npars; i++) {
    if(i<np1) {
      nunk1 += (data->ints[10]->ints[i]==0);
    } else {
      nunk2 += (data->ints[10]->ints[i]==0);
    }
  }
  int nUnk=nunk1+nunk2;

//  printf("np1 %d, np2 %d, npars %d, nu1 %d, nu2 %d \n",
  //       np1, np2, npars, nunk1, nunk2);

	double actualtheta[M], th0[npars];
	double actualsigmas[N];

	for(i=0; i<np1; i++) {
	  th0[i] = log(data->doubles[1]->doubles[i]);
	}
	for(i=0; i<np2; i++) {
	  th0[np1+i] = data->doubles[4]->doubles[i];
	}
	//printMat(th0,1,M, "theta0\n");

	if (theta) {

	  assert(!strcasecmp(data->doubles[2]->name, "sigmaprob"));
	  assert(!strcasecmp(data->doubles[3]->name, "lconst"));
	  assert(!strcasecmp(data->doubles[4]->name, "thetabase"));

	  aethetafn(npars, &theta[0], &th0[0], &data->ints[10]->ints[0],
             &data->ints[11]->ints[0], &actualtheta[0]);

	  //printMat(actualtheta,1,M, "actualtheta\n");

	  for(i=0; i<N; i++) {
	    actualsigmas[i] = exp(actualtheta[i]);
	  }

	  /*
		 if (debug > 9999) {
		   printMat(theta, 1, k, "theta\n");
	//	   printMat(&data->doubles[4]->doubles[0], 1,
    //          data->doubles[4]->len, "theta base\n");
		   printMat(actualtheta, 1, M, "actual theta\n");
//		   printMat(data->mats[0]->x, ne, ne, "I.5\n");
		 }
*/

	} else {
		for (i = 0; i < N; i++) {
			actualsigmas[i] = NAN;
		}
		for (i = 0; i < M; i++) {
			actualtheta[i] = NAN;
		}
	}

	switch (cmd) {
	case INLA_CGENERIC_GRAPH:
	{
		k = 2;
		ret = Calloc(k + 2 * M, double);
		ret[0] = N;				       /* dimension */
		ret[1] = M;
		/*
		 * number of (i <= j)
		 */
		for (i = 0; i < M; i++) {
			ret[k] = ii->ints[i];
			ret[M + k] = jj->ints[i];
			k++;
		}

	}
		break;
	case INLA_CGENERIC_Q:
	{
		int offset = 2;
		ret = Calloc(offset + M, double);
		// memset(ret + offset, 0, M * sizeof(double));
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;				       /* REQUIRED */

		int m2 = N * (N + 1) / 2;
		double ll[N * N], qtemp[m2];

		double d0[N];
		if(data->n_doubles>5) {
		  for(i = 0; i<N; i++)
		    d0[i] = data->doubles[5]->doubles[i];
		} else {
		  for(i = 0; i<N; i++)
		    d0[i] = (double)(N-i);
		}
		//printMat(d0,1,N,"d0:\n");

		// star L with diag, off-diag are zero
		k = 0;
		for (i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				if (i == j) {
					ll[k] = d0[i];
				} else {
					ll[k] = 0.0;
				}
				k++;
			}
		}

		/*
		 if (debug > 9999) {
			printMat(ll, N, N, "L[i,j]:\n");
		 }*/

		// add low theta to L
		k = 0;
		for (i = 0; i < M; i++) {
			if (ii->ints[i] != jj->ints[i]) {
				ll[iuq->ints[i]] = actualtheta[N+k++];
			}
		}
		//printMat(ll, N, N, "L\n");

/*
		 if (debug > 9999) {
			printMat(ll, N, N, "L[i,j]:\n");
		 } */

		if (nfi > 0) {

		  /*
	    if (debug > 9999) {
				printf("filling %d entries\n", nfi);
			}*/

			fillL(&N, &nfi, &ifi->ints[0], &jfi->ints[0], &ll[0]);

		  /*
		   if (debug > 9999) {
				printMat(ll, N, N, "L[i,j]:\n");
		   }*/


		}
		//printMat(ll, N, N, "filled L\n");

		// copy L to q (to be worked later in-place)
		k = 0;
		for (i = 0; i < N; i++) {
			for (j = i; j < N; j++) {
				qtemp[k++] = ll[N * i + j];
			}
		}

		/*
		 if (debug > 9999) {
			printf("L0 (upper)\n");
			k = 0;
			for (i = 0; i < N; i++) {
				for (j = i; j < N; j++) {
					printf("%2.3f ", qtemp[k++]);
				}
				printf("\n");
			}
		 } */


		// chol2inv: to compute V0 = Q_0^{-1}
		int info;
		char uplo = 'L';
		dpptri_(&uplo, &N, &qtemp[0], &info, F_ONE);

		/*
		 if (debug > 9999) {
			printf("V0 (upper)\n");
			k = 0;
			for (i = 0; i < N; i++) {
				for (j = i; j < N; j++) {
					printf("%2.3f ", qtemp[k++]);
				}
				printf("\n");
			}
		}
*/

		// si = diag(V0)^{1/2}
		// C = diag(1/si) V0 diag(1/si)
		double si[N];
		k = 0;
		for (i = 0; i < N; i++) {
			si[i] = actualsigmas[i] / sqrt(qtemp[k]);
			k += (N - i);
		}

		/*
		 if (debug > 999) {
			printMat(si, 1, N, "si:\n");
		 } */

		k = 0;
		for (i = 0; i < N; i++) {
			for (j = i; j < N; j++) {
				qtemp[k++] *= (si[i] * si[j]);
			}
		}

		/*
		 if (debug > 9999) {
			printf("V (upper)\n");
			k = 0;
			for (i = 0; i < N; i++) {
				for (j = i; j < N; j++) {
					printf("%2.3f ", qtemp[k++]);
				}
				printf("\n");
			}
		 } */

		// chol(V)
		dpptrf_(&uplo, &N, &qtemp[0], &info, F_ONE);

		/*
		if (debug > 9999) {
			printf("chol(V) (upper)\n");
			k = 0;
			for (i = 0; i < N; i++) {
				for (j = i; j < N; j++) {
					printf("%2.3f ", qtemp[k++]);
				}
				printf("\n");
			}
		} */


		// Q = chol2inv(chol(V))
		dpptri_(&uplo, &N, &qtemp[0], &info, F_ONE);

		/*
		 if (debug > 9999) {
			printf("Q (upper)\n");
			k = 0;
			for (i = 0; i < N; i++) {
				for (j = i; j < N; j++) {
					printf("%2.3f ", qtemp[k++]);
				}
				printf("\n");
			}
		 } */


		// copy the non-zero to return
		for (i = 0; i < M; i++) {
			ret[offset + i] = qtemp[iuqpac->ints[i]];
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
		for (i = 0; i < nUnk; i++) {
			ret[1 + i] = 0.0;
		}
	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		ret = Calloc(1, double);
	  ret[0] = 0.0;

	  // PC prior for UNKNOWN sigma[i]
	  if(nunk1>0) {
	    double lam;
	    for (i = 0; i < nunk1; i++) {
	      if (data->ints[10]->ints[i]) {
	        k = data->ints[11]->ints[i];
	        lam = -log(data->doubles[1]->doubles[k]);
	        lam /= data->doubles[2]->doubles[k];
	        ret[0] += pclogsigma(theta[i], lam);
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
		  ret[0] += pcmultivar(
		    nunk2, data->doubles[0]->doubles[0], &data->doubles[4]->doubles[0],
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
