
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
	int i, j, k, N, M, M2, n2, np, m2;

	// the size of the model
	assert(data->n_ints > 1);
	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	N = data->ints[0]->ints[0];
	assert(N > 0);
	M2 = N * N;
	M = N + (int)((double)N * ((double)(N-1)) / 2.0);

	assert(!strcasecmp(data->ints[1]->name, "debug"));	       // this will always be the case
	int debug = data->ints[0]->ints[0];

	assert(!strcasecmp(data->ints[2]->name, "p"));     // this will always be the case
	np = data->ints[2]->ints[0];
	assert(np > 0);
	n2 = N + np;
	m2 = (n2 * n2);
	if(debug) {
	  printf("(N = %d, M = %d, uM = %d, n2 = %d, m2 = %d)\n", N, M, M2, n2, m2) ;
	}
	double q2[np];

	assert(!strcasecmp(data->ints[3]->name, "i2th"));     // this will always be the case
	inla_cgeneric_vec_tp *i2th = data->ints[3];
	assert(i2th->len == np);

	assert(!strcasecmp(data->ints[4]->name, "iq2th"));     // this will always be the case
	inla_cgeneric_vec_tp *iq2th = data->ints[4];
	assert(i2th->len == iq2th->len);

	assert(!strcasecmp(data->ints[5]->name, "i1th"));     // this will always be the case
	inla_cgeneric_vec_tp *i1th = data->ints[5];

	assert(!strcasecmp(data->ints[6]->name, "iq1th"));     // this will always be the case
	inla_cgeneric_vec_tp *iq1th = data->ints[6];
	assert(i1th->len == iq1th->len);

	assert(!strcasecmp(data->ints[7]->name, "iq1ch"));     // this will always be the case
	inla_cgeneric_vec_tp *iq1ch = data->ints[7];

	assert(!strcasecmp(data->ints[8]->name, "ii"));     // this will always be the case
	inla_cgeneric_vec_tp *ii = data->ints[8];

	assert(!strcasecmp(data->ints[9]->name, "jj"));     // this will always be the case
	inla_cgeneric_vec_tp *jj = data->ints[9];

	assert(!strcasecmp(data->doubles[0]->name, "sch"));     // this will always be the case
	inla_cgeneric_vec_tp *sch = data->doubles[0];
	assert(iq1ch->len == sch->len);

	assert(!strcasecmp(data->doubles[1]->name, "sth"));     // this will always be the case
	inla_cgeneric_vec_tp *sth = data->doubles[1];
	assert(i1th->len == sth->len);

	assert(!strcasecmp(data->doubles[2]->name, "q"));
	inla_cgeneric_vec_tp *q = data->doubles[2];
	assert(q->len == m2);

	assert(!strcasecmp(data->doubles[3]->name, "lambda"));
	double lambda = data->doubles[1]->doubles[3];
	assert(lambda>0);

	assert(!strcasecmp(data->doubles[4]->name, "slambdas"));
	inla_cgeneric_vec_tp *slambdas = data->doubles[4];
	assert(slambdas->len>0);
	assert(slambdas->len == N);

	if(theta) {
	  if(debug>1) {
	    for (i=0; i<N; i++)
	      printf("theta[%d] = %2.1f \n", i, theta[i]);
	  }
	  for(i=0; i<np; i++) {
	    q2[i] = exp(-2 * theta[N +i]);
	    if(debug>1) {
	      printf("%2.1f %2.1f\n", theta[N+i], q2[i]);
	    }
	  }
	} else {
	  for(i=0; i<np; i++) {
	    q2[i] = NAN;
	  }
	}

	switch (cmd) {
	case INLA_CGENERIC_GRAPH:
	{
		k = 2;
	  ret = Calloc(k + 2 * M, double);
	  ret[0] = N;                                    /* dimension */
    ret[1] = M;                                    /* number of (i <= j) */
    for (i = 0; i < M; i++) {
      assert(ii->ints[i] <= jj->ints[i]);
      ret[k++] = ii->ints[i];
    }
    for (i = 0; i < M; i++) {
      ret[k++] = jj->ints[i];
      if(debug>1){
        printf("%d %2.1f %2.1f\n", i, ret[k-M], ret[k]);
      }
    }

	}
	  break;
	case INLA_CGENERIC_Q:
	{
		int offset = 2;
		ret = Calloc(offset + M, double);
//		memset(ret + offset, 0, M * sizeof(double));
		ret[0] = -1;				       /* REQUIRED */
    ret[1] = M;				       /* REQUIRED */

		int ipiv[N], ipiv2[n2], info=0;

		double qq[m2], vv[m2], vc[M];

		for(i=0; i<m2; i++) {
		  qq[i] = q->doubles[i];
		}
		for(i=0; i<sch->len; i++){
		  qq[iq1ch->ints[i]] *= sch->doubles[i];
		}
		for(i=0; i<iq2th->len; i++) {
		  qq[iq2th->ints[i]] += q2[i2th->ints[i]];
		}
		k = 0;
		for(i=0; i<i1th->len; i++) {
		  qq[iq1th->ints[i]] = -sth->doubles[i] * q2[i1th->ints[i]];
		}
		if(debug>1){
		  printf("Q[i,j]:\n");
		  k=0;
		  for(i=0; i<n2; i++) {
		    for(j=0; j<n2; j++) {
		      printf("%2.1f ", qq[k]);
  		    k++;
	  	  }
		    if(debug>1){
		      printf("\n");
		    }
		  }
		  printf("vv[i,j]:\n");
		}
		k=0;
		for(i=0; i<n2; i++) {
		  for(j=0; j<n2; j++) {
		    if(i==j) {
		      vv[k] = 1.0;
		    } else {
		      vv[k] = 0;
		    }
		    if(debug>1){
		      printf("%2.3f ", vv[k]);
		    }
		    k++;
		  }
		  if(debug>1){
		    printf("\n");
		  }
		}
		assert(k==m2);

		dgesv_(&n2, &n2, qq, &n2, ipiv2, vv, &n2, &info, F_ONE);

		if(debug>1){
		  printf("vv[i,j]:\n");
		  k=0;
		  for(i=0; i<n2; i++) {
		    for(j=0; j<n2; j++) {
		      printf("%2.1f ", vv[k]);
		      k++;
		    }
		    if(debug>1){
		      printf("\n");
		    }
		  }
		}

		k=0;
		int kc = 0;
		double ss[N];
		for(i=0; i<N; i++) {
		  for(j=0; j<n2; j++) {
		    if(j<N) {
		      vc[kc] = vv[k];
		      if(i==j) ss[i] = sqrt(vc[kc]);
		      kc++;
		    }
		    k++;
		  }
		  if(debug>1)
		    printf("ss[%d] = %2.1f\n", i, ss[i]) ;
		}
		if(debug>1)
		  printf("\n") ;
		assert(kc==M2);

		if(debug>1){
		  printf("vc[i,j]:\n");
		  k=0;
		  for(i=0; i<N; i++) {
		    for(j=0; j<N; j++) {
		      printf("%2.1f ", vc[k]);
		      k++;
		    }
		    if(debug>1){
		      printf("\n");
		    }
		  }
		}

		k=0;
    for(i=0; i<N; i++) {
      for(j=0; j<N; j++) {
        if(i==j) {
          qq[k] = exp(theta[i]);
        } else {
          qq[k] = 0.0;
        }
        vc[k] = exp(theta[i]+theta[j]) * vc[k] / (ss[i] * ss[j]);
		    k++;
		  }
		}
		dgesv_(&N, &N, vc, &N, ipiv, qq, &N, &info, F_ONE);

    k = 0;
    kc = 0;
    for(i=0; i<N; i++) {
      for(j=0; j<N; j++) {
        if(j>=i) {
          ret[offset + kc] = qq[k];
          kc++;
        }
        k++;
      }
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
		ret = Calloc(n2 + 1, double);
	  ret[0] = n2;
	  for(i = 0; i < n2; i++) {
			ret[1+i] = 1.0;
		}
	}
    break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		ret = Calloc(1, double);
		// PC-priors
		ret[0] = 0.0;
		double lam = 0;
		for(i = 0; i < N; i++) {
			lam = slambdas->doubles[i];
		  ret[0] += log(lam) + theta[i] - lam * exp(theta[i]);
		}
		for(i = 0; i < np; i++) {
		  ret[0] += -0.5 * SQR(theta[N+i]);
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
