
/* cgeneric_dag_cholQ.c
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
 *        Elias T Krainski
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 */

#include "corgraphs.h"

void fillL(int *d, int *m, int *ii, int *jj, double *x) {
  // fill lower triangle
  int i, j, l, k, p, n = *d, nij = *m;
  for(l=0; l<nij; l++) {
    i = ii[l];
    j = jj[l];
    //  printf("%d: ", l);
    if(j>0) {
      p = j*n + i;
      for(k=0; k<j; k++) {
        //        printf("%d: i %d, j %d, k %d\n", p, i, j, k);
        x[p] -= x[k*n+i] * x[k*n+j] / x[j*n+j];
      }
    }
  }
}

double *inla_cgeneric_dag_cholQ(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
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

  assert(!strcasecmp(data->ints[3]->name, "nnz"));	       // this will always be the case
  int nnz = data->ints[3]->ints[0];
  assert(M == nnz);

  assert(!strcasecmp(data->ints[4]->name, "nfi"));	       // this will always be the case
  int nfi = data->ints[4]->ints[0];

  assert(!strcasecmp(data->ints[5]->name, "ii"));     // this will always be the case
  inla_cgeneric_vec_tp *ii = data->ints[5];
  assert(M == ii->len);

  assert(!strcasecmp(data->ints[6]->name, "jj"));     // this will always be the case
  inla_cgeneric_vec_tp *jj = data->ints[6];
  assert(M == jj->len);

  assert(!strcasecmp(data->ints[7]->name, "ilq"));     // this will always be the case
  inla_cgeneric_vec_tp *ilq = data->ints[7];
  assert(M == ilq->len);

  assert(!strcasecmp(data->ints[8]->name, "iuq"));     // this will always be the case
  inla_cgeneric_vec_tp *iuq = data->ints[8];
  assert(M == iuq->len);

  assert(!strcasecmp(data->ints[9]->name, "ifi"));     // this will always be the case
  inla_cgeneric_vec_tp *ifi = data->ints[9];
  assert(nfi == ifi->len);

  assert(!strcasecmp(data->ints[10]->name, "jfi"));     // this will always be the case
  inla_cgeneric_vec_tp *jfi = data->ints[10];
  assert(nfi == jfi->len);

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
    printf("\n ilq:\n");
    for(i=0; i<ilq->len; i++) {
      printf("%d ", ilq->ints[i]);
    }
    printf("\n iuq:\n");
    for(i=0; i<iuq->len; i++) {
      printf("%d ", iuq->ints[i]);
    }

    printf("\n ifi:\n");
    for(i=0; i<ifi->len; i++) {
      printf("%d ", ifi->ints[i]);
    }
    printf("\n jfi:\n");
    for(i=0; i<jfi->len; i++) {
      printf("%d ", jfi->ints[i]);
    }

  }

  double param[M];
  if(theta) {
    k = 0;
    for(i=0; i<M; i++) {
      if(ii->ints[i]==jj->ints[i]) {
        param[i] = exp(theta[k++]);
      }
    }
    if(ne>0) {
      for(i=0; i<M; i++) {
        if(ii->ints[i]!=jj->ints[i]) {
          param[i] = theta[k++];
        }
      }
    }
  } else {
    for(i=0; i<M; i++) {
      param[i] = NAN;
    }
  }

  if(debug>99){
    printf("\n params:\n");
    for(i=0; i<M; i++) {
      printf("%2.2f ", param[i]);
    }
    printf("\n");

  }

  switch (cmd) {
  case INLA_CGENERIC_GRAPH:
  {
    k = 2;
    ret = Calloc(k + 2 * M, double);
    ret[0] = N;                                    /* dimension */
    ret[1] = M;
    /* number of (i <= j) */
    for (i=0; i<M; i++) {
        ret[k] = ii->ints[i];
        ret[M+k] = jj->ints[i];
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
    double qq[N2], ll[N2];
    k=0;
    for(i=0; i<N; i++) {
      for(j=0; j<N; j++) {
        ll[k] = 0.0;
        qq[k++] = 0.0;
      }
    }

    for(i=0; i<M; i++) {
      ll[ilq->ints[i]] = param[i];
    }
    if(debug>99) {
      printf("L[i,j]:\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=0; j<N; j++) {
          printf("%2.3f ", ll[k]);
          k++;
        }
        printf("\n");
      }
    }

    if(nfi>0) {

      fillL(&N, &nfi, &ifi->ints[0], &jfi->ints[0], &ll[0]) ;

    }
    if(debug>99) {
      printf("L[i,j]:\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=0; j<N; j++) {
          printf("%2.3f ", ll[k]);
          k++;
        }
        printf("\n");
      }
    }

    char transa = 'N';
    char transb = 'T';
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

    k=2;
    for(i=0; i<M; i++) {
      ret[k] = qq[iuq->ints[i]];
      k++;
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
