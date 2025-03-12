
/* cgeneric_corgraph.c
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

double *inla_cgeneric_corgraph(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data)
{

  double *ret = NULL;
  int i, j, k, N, M;

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
  M = N + ne;

  assert(!strcasecmp(data->ints[3]->name, "nfi"));	       // this will always be the case
  int nfi = data->ints[3]->ints[0];

  assert(!strcasecmp(data->ints[4]->name, "ii"));     // this will always be the case
  inla_cgeneric_vec_tp *ii = data->ints[4];
  assert(M == ii->len);

  assert(!strcasecmp(data->ints[5]->name, "jj"));     // this will always be the case
  inla_cgeneric_vec_tp *jj = data->ints[5];
  assert(M == jj->len);

  assert(!strcasecmp(data->ints[6]->name, "iuq"));     // this will always be the case
  inla_cgeneric_vec_tp *iuq = data->ints[6];
  assert(M == iuq->len);

  assert(!strcasecmp(data->ints[7]->name, "iuqpac"));     // this will always be the case
  inla_cgeneric_vec_tp *iuqpac = data->ints[7];
  assert(M == iuqpac->len);

  assert(!strcasecmp(data->ints[8]->name, "ifi"));     // this will always be the case
  inla_cgeneric_vec_tp *ifi = data->ints[8];
  assert(nfi == ifi->len);

  assert(!strcasecmp(data->ints[9]->name, "jfi"));     // this will always be the case
  inla_cgeneric_vec_tp *jfi = data->ints[9];
  assert(nfi == jfi->len);

  assert(!strcasecmp(data->ints[10]->name, "itheta"));     // this will always be the case
  inla_cgeneric_vec_tp *itheta = data->ints[10];
  assert(M == itheta->len);
  int nparams[3];
  nparams[0] = itheta->ints[N-1]+1;
  nparams[2] = itheta->ints[M-1]+1;
  nparams[1] = nparams[2]-nparams[0];

  assert(!strcasecmp(data->ints[11]->name, "sfixed"));     // this will always be the case
  assert(N == data->ints[11]->len);
  int sfixed[N];
  for(i=0; i<N; i++) {
    sfixed[i] = data->ints[11]->ints[i];
  }

  assert(!strcasecmp(data->doubles[0]->name, "lambda"));
  double lambda = data->doubles[0]->doubles[0];
  assert(lambda>0.0);

  assert(!strcasecmp(data->doubles[1]->name, "sigmaref"));
  inla_cgeneric_vec_tp *sigmaref = data->doubles[1];
  assert(sigmaref->len>0);
  assert(sigmaref->len == N);

  assert(!strcasecmp(data->doubles[2]->name, "sigmaprob"));
  inla_cgeneric_vec_tp *sigmaprob = data->doubles[2];
  assert(sigmaprob->len>0);
  assert(sigmaprob->len == N);

  assert(!strcasecmp(data->doubles[3]->name, "lconst"));
  double lconst = data->doubles[3]->doubles[0];

  assert(!strcasecmp(data->doubles[4]->name, "thetabasescaled"));
  assert(data->doubles[4]->length==ne);

  assert(!strcasecmp(data->mats[0]->name, "hHneg"));
  assert(data->mats[0]->nrow==ne);
  assert(data->mats[0]->ncol==ne);

  double lowtheta[ne], thetat[ne];
  double sigmas[N];
  int nsu = 0;
  if(theta) {
    for(i=0; i<N; i++) {
      assert(sigmaref->doubles[i]>0);
      if(sfixed[i]) {
        sigmas[i] = sigmaref->doubles[itheta->ints[i]];
      } else {
        sigmas[i] = exp(theta[itheta->ints[i]]);
        nsu++;
      }
    }
    if(debug>99) {
      printf("n sigma to be estimated = %d", nsu);
      printMat(sfixed,1,N,"sfixed\n");
      printMat(sigmas,1,N,"sigmas\n");
    }
    // this is I(\theta_0)^{-0.5} \theta_0
    for(i=0; i<ne; i++) {
      thetat[i] = data->doubles[4]->doubles[i];
      lowtheta[i] = theta[itheta->ints[N+i]];
    }
    if(debug>99){
      printMat(thetat, 1, ne, "thetat:\n");
      printMat(lowtheta, 1, ne, "lowtheta:\n");
      printMat(data->mats[0]->x, ne, ne, "hHneg\n");
    }
    int one=1;
    char trans = 'N';
    double alpha = 1.0, beta = -1.0;
    // (complete) thetat: I(\theta_0)^{-0.5} ( \theta - \theta_0)
    dgemv_(&trans, &ne, &ne, &alpha,
           &data->mats[0]->x[0], &ne,
           &lowtheta[0], &one,
           &beta, &thetat[0], &one, F_ONE);
    if(debug>99) {
      printMat(thetat, 1, ne, "actual thetat:\n");
    }
  } else {
    for(i=0; i<N; i++) {
      sigmas[i] = NAN;
    }
    for(i=0; i<ne; i++) {
      lowtheta[i] = NAN;
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

    int m2 = N*(N+1)/2;
    double ll[N * N], qtemp[m2];

    k=0;
    for(i=0; i<N; i++) {
      for(j=0; j<N; j++) {
        if(i==j) {
          ll[k] = 1.0;
        } else {
          ll[k] = 0.0;
        }
        k++;
      }
    }

    if(debug>99) {
      printMat(ll, N, N, "L[i,j]:\n");
    }

    k=0;
    for(i=0; i<M; i++) {
      if(ii->ints[i]!=jj->ints[i]) {
        ll[iuq->ints[i]] = lowtheta[k++];
      }
    }
    if(debug>99) {
      printMat(ll, N, N, "L[i,j]:\n");
    }

    if(nfi>0) {

      if(debug>99) {
        printf("filling %d entries\n", nfi);
      }

      fillL(&N, &nfi, &ifi->ints[0], &jfi->ints[0], &ll[0]) ;

      if(debug>99) {
        printMat(ll, N, N, "L[i,j]:\n");
      }

    }

    // copy L to q (to be worked later in-place)
    k=0;
    for(i=0; i<N; i++) {
      for(j=i; j<N; j++) {
        qtemp[k++] = ll[N*i+j];
      }
    }

    if(debug>99){
      printf("Q0 (upper)\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=i; j<N; j++) {
          printf("%2.3f ", qtemp[k++]);
        }
        printf("\n");
      }
    }

    // chol2inv: to compute V0 = Q_0^{-1}
    int info;
    char uplo = 'L';
    dpptri_(&uplo, &N, &qtemp[0], &info, F_ONE);

    if(debug>99){
      printf("V0 (upper)\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=i; j<N; j++) {
          printf("%2.3f ", qtemp[k++]);
        }
        printf("\n");
      }
    }

    // si = diag(V0)^{1/2}
    // C = diag(1/si) V0 diag(1/si)
    double si[N];
    k=0;
    for(i=0; i<N; i++) {
      si[i] = sigmas[i]/sqrt(qtemp[k]);
      k += (N-i);
    }

    if(debug>99) {
      printMat(si, 1, N, "si:\n");
    }
    k=0;
    for(i=0; i<N; i++) {
      for(j=i; j<N; j++) {
        qtemp[k++] *= (si[i] * si[j]);
      }
    }

    if(debug>99){
      printf("V (upper)\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=i; j<N; j++) {
          printf("%2.3f ", qtemp[k++]);
        }
        printf("\n");
      }
    }

    // chol(V)
    dpptrf_(&uplo, &N, &qtemp[0], &info, F_ONE);

    if(debug>99){
      printf("chol(V) (upper)\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=i; j<N; j++) {
          printf("%2.3f ", qtemp[k++]);
        }
        printf("\n");
      }
    }

    // Q = chol2inv(chol(V))
    dpptri_(&uplo, &N, &qtemp[0], &info, F_ONE);

    if(debug>99){
      printf("Q (upper)\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=i; j<N; j++) {
          printf("%2.3f ", qtemp[k++]);
        }
        printf("\n");
      }
    }

    // copy the non-zero to return
    for(i=0; i<M; i++) {
      ret[offset+i] = qtemp[iuqpac->ints[i]];
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
    ret = Calloc(nparams[2] + 1, double);
    ret[0] = nparams[2];
    for(i = 0; i < nparams[2]; i++) {
      ret[1+i] = 0.0;
    }
  }
    break;

  case INLA_CGENERIC_LOG_PRIOR:
  {
    ret = Calloc(1, double);

    // the log prior:
    // lconst should be equal to
    //    log(lambda) -(m-1)*log(pi)-log(2)-log(|H|)
    ret[0] = lconst;

    // PC prior for sigma[i]
    double lam, val, pparams[ne];
    for(i = 0; i < nparams[0]; i++) {
      if(sfixed[i]>0) {
        lam = -log(sigmaprob->doubles[i]) / sigmaref->doubles[i];
        ret[0] += pclogsigma(log(sigmas[i]), lam);
        if(debug>999) {
          printf("lamb[%d] = %2.3f, p %2.3f \n", i, lam, ret[0]);
        }
      }
    }

    // TO BE FIXED when nparams[2]<ne?
    pparams[ne-1] = atan2(thetat[ne-1], theta[ne-2]);
    if(pparams[ne-1]<0) {
      pparams[ne-1] += 2.0*M_PI;
    }
    val = SQR(thetat[ne-1]) + SQR(thetat[ne-2]);
    for(i=(ne-2); i>=0; i--) {
      pparams[i] = atan2(sqrt(val), thetat[i-1]);
      val += SQR(thetat[i-1]);
    }
    // (\theta-\theta_0)I(\theta_0)(\theta -\theta_0)
    pparams[0] = sqrt(val); // the approx. KLD

    double ldJacobian;
    ldJacobian = ((double)(ne-1)) * log(pparams[0]);
    if(ne>2) {
      for(i=1; i<(ne-1); i++) { // not the last one
        ldJacobian += ((double)(ne-1-i)) * log(sin(pparams[i]));
      }
    }
    if(debug>999) {
      printf("log det Jacobian = %2.7f\n", ldJacobian);
    }

    ret[0] += ldJacobian - lambda * pparams[0];

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
