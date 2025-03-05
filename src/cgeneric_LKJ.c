
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

double *inla_cgeneric_LKJ(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
{

  // Lewandowski-Kurowicka-Joe - LKJ prior
  // This is the cgeneric implementatin of the
  // LKJ prior for a correlation matrix C with dimension N,
  //  given (scalar) parameter 'eta'.
  // Parametrized from a hypershere decomposition,
  //  see Rapisarda, Brigo and Mercurio (2007).
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
  M = (int)((double)N * ((double)(N+1)) / 2.0);
  nth = (int)((double)N * ((double)(N-1)) / 2.0);

  assert(!strcasecmp(data->ints[1]->name, "debug"));	       // this will always be the case
  int debug = data->ints[1]->ints[0];

  assert(!strcasecmp(data->doubles[0]->name, "eta"));
  double eta = data->doubles[0]->doubles[0];
  assert(eta>0);

  assert(!strcasecmp(data->doubles[1]->name, "lc"));
  double lc = data->doubles[1]->doubles[0];

  if(debug>999) {
    printf("N=%d, nth=%d, M=%d, eta=%f, lc=%f\n", N, nth, M, eta, lc);
  }

  switch (cmd) {
  case INLA_CGENERIC_GRAPH:
  {
    k = 2;
    ret = Calloc(k + 2 * M, double);
    ret[0] = N;                         /* dimension */
    ret[1] = M;                         /* number of (i <= j) */

      for (i = 0; i < N; i++) {
        for(j = i; j<N; j++) {
          ret[M+k] = j;
          ret[k++] = i;
        }
      }

      if(debug>999) {
        printf("ii: ");
        for(k=0; k<M; k++) {
          printf("%f ", ret[2+k]);
        }
        printf("\njj: ");
        for(k=0; k<M; k++) {
          printf("%f ", ret[2+M+k]);
        }
        printf("\n");
      }


  }
    break;
  case INLA_CGENERIC_Q:
  {
    // Q = (CC)^{-1}
    int offset = 2;

    ret = Calloc(offset + M, double);
    ret[0] = -1;			       /* REQUIRED */
    ret[1] = M;				       /* REQUIRED */

    double hld;
    theta2gamma2Lcorr(N, &hld, &theta[0], &ret[2]);

    if(debug>999) {
      printf("L:\n");
      for(i=0; i<N; i++) {
        k = i;
        for(j=0; j<=i; j++) {
          printf("%2.3f ", ret[offset+k]);
          k += (N-j-1);
        }
        printf("\n");
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
    ret = Calloc(nth + 1, double);
    ret[0] = nth;
    for(i=0; i<nth; i++) {
        ret[1+i] = 3.0;
    }
  }
    break;

  case INLA_CGENERIC_LOG_PRIOR:
  {
    ret = Calloc(1, double);

    // precomp constant
    ret[0] = lc;

    double daux;
    // add the Jacobian of theta[i] -> gamma[i]
    for(i=0; i<nth; i++) {
      daux = exp(-theta[i]);
      if(debug>9999)
        printf("theta[%d] = %2.5f, J = %2.5f\n",
               i, theta[i], M_PI * daux / SQR(1.0 + daux));
      ret[0] += M_PI * daux / SQR(1.0 + daux);
    }
    if(debug>99) {
      printf("Jacobian = %2.3f, %2.3f ",
             ret[0]-lc, ret[0]);
    }
    // log determinant
    double ll[M];
    theta2gamma2Lcorr(N, &daux, &theta[0], &ll[0]);
    if(daux<0) {
      daux *= -2.0;
    } else {
      daux *= 2.0;
    }
    if(debug>99) {
      printf("eta-1 = %2.3f, eta-1.0 = %2.3f, |R| = %2.5f, lc+(eta-1)|R| = %2.5f, %2.5f\n",
             eta-1, eta-1.0, daux, lc+(eta-1)*daux, daux * (eta-1.0));
    }
    ret[0] += daux * (eta-1.0);

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
