
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

double *inla_cgeneric_corgraphs_cholQ(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
{

  double *ret = NULL;
  int i, j, k, N, M, nparam;

  // the size of the model
  assert(data->n_ints > 1);
  assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
  N = data->ints[0]->ints[0];
  assert(N > 0);

  assert(!strcasecmp(data->ints[1]->name, "debug"));	       // this will always be the case
  int debug = data->ints[1]->ints[0];

  assert(!strcasecmp(data->ints[2]->name, "ne"));     // this will always be the case
  int ne = data->ints[2]->ints[0];
  assert(ne > 0);
  nparam = N + ne;
  M = N + ne;

  assert(!strcasecmp(data->ints[3]->name, "ii"));     // this will always be the case
  inla_cgeneric_vec_tp *ii = data->ints[3];
  assert(ne == ii->len);

  assert(!strcasecmp(data->ints[4]->name, "jj"));     // this will always be the case
  inla_cgeneric_vec_tp *jj = data->ints[4];
  assert(ii->len == jj->len);

  assert(!strcasecmp(data->ints[5]->name, "iiq"));     // this will always be the case
  inla_cgeneric_vec_tp *iiq = data->ints[5];
  assert(ii->len == iiq->len);

  assert(!strcasecmp(data->ints[6]->name, "iifi"));     // this will always be the case
  inla_cgeneric_vec_tp *iifi = data->ints[6];
  int nfi = iifi->len;

  assert(!strcasecmp(data->ints[7]->name, "ifi"));
  inla_cgeneric_vec_tp *ifi = data->ints[7];
  assert(nfi = ifi->len);

  assert(!strcasecmp(data->ints[8]->name, "jfi"));
  inla_cgeneric_vec_tp *jfi = data->ints[8];
  assert(nfi = jfi->len);

  assert(!strcasecmp(data->doubles[0]->name, "lambda"));
  double lambda = data->doubles[0]->doubles[0];
  assert(lambda>0);

  assert(!strcasecmp(data->doubles[1]->name, "slambdas"));
  inla_cgeneric_vec_tp *slambdas = data->doubles[1];
  assert(slambdas->len>0);
  assert(slambdas->len == N);

  int iprint = 0;
  if((debug>9) & (iprint<1)) {
    iprint++;
    printf("(N = %d, M = %d, ne = %d)\n", N, M, ne) ;

    printf("\n ii:\n");
    for(i=0; i<ii->len; i++) {
      printf("%d ", ii->ints[i]);
    }
    printf("\n jj:\n");
    for(i=0; i<jj->len; i++) {
      printf("%d ", jj->ints[i]);
    }

    printf("\n iiq:\n");
    for(i=0; i<iiq->len; i++) {
      printf("%d ", iiq->ints[i]);
    }
    printf("\n ifi:\n");
    for(i=0; i<ifi->len; i++) {
      printf("%d ", ifi->ints[i]);
    }

  }

  switch (cmd) {
  case INLA_CGENERIC_GRAPH:
  {
    k = 2;
    ret = Calloc(k + 2 * M, double);
    ret[0] = N;                                    /* dimension */
    ret[1] = M;
    /* number of (i <= j) */
    for (i=0; i<ne; i++) {
        ret[k] = ii->ints[i];
        ret[ne+k] = jj->ints[i];
        k++;
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

    int N2 = N*N;
    char transa = 'T';
    char transb = 'N';
    double qq[N2], ll[N2];

    k=0;
    for(i=0; i<N; i++) {
      ll[k] = exp(theta[i]);
      k += N;
    }
    for(i=0; i<iiq->len; i++) {
      ll[iiq->ints[i]] = theta[N+i];
    }

    if(nfi>0) {

      int l, ki, kj;
      for(l=0; l<iifi->len; l++) {
        i = ifi->ints[l];
        j = jfi->ints[l];
        ki = 1;
      }

    }

    double alpha = 1.0, beta = 0.0;
    dgemm_(&transa, &transb, &N, &N, &N, &alpha,
           &ll[0], &N, &ll[0], &N, &beta, &qq[0], &N, F_ONE);

    if(debug>99) {
      printf("Q[i,j]:\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=0; j<N; j++) {
          printf("%2.3f ", qq[k]);
          k++;
        }
        printf("\n");
      }
    }

    k=0;
    for(i=0; i<N; i++) {
      ret[offset+k] = qq[k];
      k += N;
    }
    for(i=0; i<iiq->len; i++) {
      ret[offset+iiq->ints[i]] = qq[iiq->ints[i]];
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
    ret = Calloc(nparam + 1, double);
    ret[0] = nparam;
    for(i = 0; i < M; i++) {
      ret[1+i] = 0.0;
    }
  }
    break;

  case INLA_CGENERIC_LOG_PRIOR:
  {
    ret = Calloc(1, double);
    ret[0] = 0.0;

    // temporary: N(0, 1/lambda_i)
    double lam;
    for(i = 0; i < N; i++) {
      lam = slambdas->doubles[i];
      ret[0] += -0.5 * pow2(theta[i]) / lam;
    }
    for(i=N; i<M; i++) {
      ret[0] += -0.5 * pow2(theta[i]) / lambda;
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
