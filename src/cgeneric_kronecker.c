
/* cgeneric_kronecker.c
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

#include <ltdl.h>
#include "cgeneric_defs.h"

double *inla_cgeneric_kronecker(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp * data) {

  double *ret1 = NULL; // to store output from model1.
  double *ret2 = NULL; // to store output from model2.
  double *ret = NULL; // to return;

  //  Q = Q1 (x) Q2

  // cmd is length 1 string, however,
  // theta contains theta[1:nth1-1, nth1:(nth1+nth2-1)]
  //   nth1 is #theta for M1, nth2 is #theta for M2
  // data contains data for model M1 and model M2.
  // data->ints[0]->ints contain
  //     {n, M,
  //      n1, ni1, nd1, nc1, nm1, nsm1,
  //      n2, ni2, nd2, nc2, nm2, nsm2}
  // n1 is n for mode1 and n2 is n for model2
  // ni1 is the number of ints for model1
  // ni2 is the number of ints for model2
  // ...
  // so that
  //  data->ints[1+1:(ni1-2)] contain ints for M1
  //  data->ints[ni1+1:(ni2-2)-1] contain ints for M2
  //  data->ints[ni1+ni2-1] contain index for the kronecker model
  //  data->doubles[1:nd1-1] contain doubles for M1
  //  data->doubles[nd1+1:nd2-1] contain doubles for M2
  //  data->chars[1+1:nc1] contain chars for M1
  //  data->chars[1+nc1+1:nc2] contain chars for M2
  //  data->mats[1:nm1-1] contain mats for M1
  //  data->mats[nm1+1:nm2-1] contain mats for M2
  //  data->smats[1:nsm1-1] contain smats for M1
  //  data->smats[nsm1+1:nsm2-1] contain smats for M2

  int N = data->ints[0]->ints[0];
  int debug = data->ints[1]->ints[0];
  int i, j, n1, n2, iaux;
  int ni1, nd1, nc1, nm1, nsm1;
  int ni2, nd2, nc2, nm2, nsm2;
  int M, M1, M2, k;

  M = data->ints[0]->ints[1];
  n1 = data->ints[0]->ints[2];
  M1 = data->ints[0]->ints[3];
  ni1 = data->ints[0]->ints[4];
  nd1 = data->ints[0]->ints[5];
  nc1 = data->ints[0]->ints[6];
  nm1 = data->ints[0]->ints[7];
  nsm1 = data->ints[0]->ints[8];

  n2 = data->ints[0]->ints[9];
  M2 = data->ints[0]->ints[10];
  ni2 = data->ints[0]->ints[11];
  nd2 = data->ints[0]->ints[12];
  nc2 = data->ints[0]->ints[13];
  nm2 = data->ints[0]->ints[14];
  nsm2 = data->ints[0]->ints[15];

  assert(ni1>1);
  assert(ni2>1);
  assert(nc1>1);
  assert(nc2>1);
  assert(data->n_chars>5);

  if(debug) {
    printf(//stderr,
      "n %d, M %d\nn1 %d, M1 %d, ni %d, nd %d, nc %d, nm %d, nsm %d \n",
      N, M, n1, M1, ni1, nd1, nc1, nm1, nsm1);
    printf(//stderr,
      "n2 %d, M2 %d, ni %d, nd %d, nc %d, nm %d, nsm %d \n",
      n2, M2, ni2, nd2, nc2, nm2, nsm2);
  }

  N = n1 * n2;
  assert(N == data->ints[0]->ints[0]);

  inla_cgeneric_data_tp *dataM1 = Calloc(1, inla_cgeneric_data_tp);
  inla_cgeneric_data_tp *dataM2 = Calloc(1, inla_cgeneric_data_tp);

  dataM1->n_ints = ni1;
  dataM2->n_ints = ni2;
  dataM1->n_doubles = nd1;
  dataM2->n_doubles = nd2;
  dataM1->n_chars = nc1;
  dataM2->n_chars = nc2;
  dataM1->n_mats = nm1;
  dataM2->n_mats = nm2;
  dataM1->n_smats = nsm1;
  dataM2->n_smats = nsm2;

  dataM1->ints = Calloc(dataM1->n_ints, inla_cgeneric_vec_tp *);
  dataM2->ints = Calloc(dataM2->n_ints, inla_cgeneric_vec_tp *);
  dataM1->doubles = Calloc(dataM1->n_doubles, inla_cgeneric_vec_tp *);
  dataM2->doubles = Calloc(dataM2->n_doubles, inla_cgeneric_vec_tp *);
  dataM1->chars = Calloc(dataM1->n_chars, inla_cgeneric_vec_tp *);
  dataM2->chars = Calloc(dataM2->n_chars, inla_cgeneric_vec_tp *);
  dataM1->mats = Calloc(dataM1->n_mats, inla_cgeneric_mat_tp *);
  dataM2->mats = Calloc(dataM2->n_mats, inla_cgeneric_mat_tp *);
  dataM1->smats = Calloc(dataM1->n_smats, inla_cgeneric_smat_tp *);
  dataM2->smats = Calloc(dataM2->n_smats, inla_cgeneric_smat_tp *);

  for(i=0; i<6; i++) {
    printf("%d %s \n", i, &data->chars[i]->chars[0]);
  }

  // copy ints for M1
  if(debug)
    printf("copy ints for M1 ...");
  for(i=0; i<2; i++) {
    if(debug>1)
      printf(" %d, ", i);
    dataM1->ints[i] = Calloc(1, inla_cgeneric_vec_tp);
    dataM1->ints[i]->name = data->ints[i]->name;
    dataM1->ints[i]->len = 1;
    dataM1->ints[i]->ints = Calloc(1, int);
  }
  dataM1->ints[0]->ints[0] = n1;
  dataM1->ints[1]->ints[0] = data->ints[1]->ints[0];
  if(ni1>2) {
    for(i=2; i<dataM1->n_ints; i++) {
      if(debug>1)
        printf(" %d, ", i);
      dataM1->ints[i] = Calloc(1, inla_cgeneric_vec_tp);
      dataM1->ints[i]->name = data->ints[i]->name;
      dataM1->ints[i]->len = data->ints[i]->len;
      dataM1->ints[i]->ints = Calloc(dataM1->ints[i]->len, int);
      for(j=0; j<dataM1->ints[i]->len; j++) {
        dataM1->ints[i]->ints[j] = data->ints[i]->ints[j];
      }
    }
  }
  if(debug)
    printf("done \n");

  // copy ints for M2
  if(debug)
    printf("copy ints for M2 ...");
  for(i=0; i<2; i++) {
    dataM2->ints[i] = Calloc(1, inla_cgeneric_vec_tp);
    dataM2->ints[i]->name = data->ints[i]->name;
    dataM2->ints[i]->len = 1;
    dataM2->ints[i]->ints = Calloc(1, int);
    printf(" %d, ", i);
  }
  dataM2->ints[0]->ints[0] = n2;
  dataM2->ints[1]->ints[0] = data->ints[1]->ints[0];
  if(ni2>2) {
    for(i=2; i<dataM2->n_ints; i++) {
      if(debug>1)
        printf(" %d, ", i);
      dataM2->ints[i] = Calloc(1, inla_cgeneric_vec_tp);
      dataM2->ints[i]->name = data->ints[ni1-2+i]->name;
      dataM2->ints[i]->len = data->ints[ni1-2+i]->len;
      dataM2->ints[i]->ints = Calloc(dataM2->ints[i]->len, int);
      for(j=0; j<dataM2->ints[i]->len; j++) {
        dataM2->ints[i]->ints[j] = data->ints[ni1-2+i]->ints[j];
      }
    }
  }
  if(debug)
    printf("done \n");

  if(nd1>0) { // copy doubles for M1
    if(debug)
      printf("copy doubles for M1 ...");
    for(i=0; i<dataM1->n_doubles; i++) {
      if(debug>1)
        printf(" %d, ", i);
      dataM1->doubles[i] = Calloc(1, inla_cgeneric_vec_tp);
      dataM1->doubles[i]->name = data->doubles[i]->name;
      dataM1->doubles[i]->len = data->doubles[i]->len;
      dataM1->doubles[i]->doubles = Calloc(dataM1->doubles[i]->len, double);
      for(j=0; j<dataM1->doubles[i]->len; j++) {
        dataM1->doubles[i]->doubles[j] = data->doubles[i]->doubles[j];
      }
    }
    if(debug)
      printf("done \n");
  }

  if(nd2>0) { // copy doubles for M2
    if(debug)
      printf("copy doubles for M2 ...");
    for(i=0; i<dataM2->n_doubles; i++) {
      dataM2->doubles[i] = Calloc(1, inla_cgeneric_vec_tp);
      dataM2->doubles[i]->name = data->doubles[nd1+i]->name;
      dataM2->doubles[i]->len = data->doubles[nd1+i]->len;
      dataM2->doubles[i]->doubles = Calloc(dataM2->doubles[i]->len, double);
      for(j=0; j<dataM2->doubles[i]->len; j++) {
        dataM2->doubles[i]->doubles[j] = data->doubles[nd1+i]->doubles[j];
      }
    }
    if(debug)
      printf("done \n");
  }

  if(nc1>0) { // copy chars for M1
    if(debug)
      printf("copy chars for M1 ...");
    for(i=0; i<dataM1->n_chars; i++) {
      if(debug>1)
        printf(" %d, ", i);
      dataM1->chars[i] = Calloc(1, inla_cgeneric_vec_tp);
      dataM1->chars[i]->name = data->chars[2+i]->name;
      dataM1->chars[i]->len = data->chars[2+i]->len;
      dataM1->chars[i]->chars = Calloc(dataM1->chars[i]->len, char);
      for(j=0; j<dataM1->chars[i]->len; j++) {
        dataM1->chars[i]->chars[j] = data->chars[2+i]->chars[j];
      }
    }
    if(debug)
      printf("done \n");
  }

  if(nc2>0) { // copy chars for M2
    if(debug)
      printf("copy chars for M2 ...");
    for(i=0; i<dataM2->n_chars; i++) {
      if(debug>1)
        printf(" %d, ", i);
      dataM2->chars[i] = Calloc(1, inla_cgeneric_vec_tp);
      dataM2->chars[i]->name = data->chars[2+nc1+i]->name;
      dataM2->chars[i]->len = data->chars[2+nc1+i]->len;
      dataM2->chars[i]->chars = Calloc(dataM2->chars[i]->len, char);
      for(j=0; j<dataM2->chars[i]->len; j++) {
        dataM2->chars[i]->chars[j] = data->chars[2+nc1+i]->chars[j];
      }
    }
    if(debug)
      printf("done \n");
  }

  if(nm1>0) { // copy mats for M1
    if(debug)
      printf("copy mats for M1 ...");
    for(i=0; i<dataM1->n_mats; i++) {
      if(debug>1)
        printf(" %d, ", i);
      dataM1->mats[i] = Calloc(1, inla_cgeneric_mat_tp);
      dataM1->mats[i]->name = data->mats[i]->name;
      dataM1->mats[i]->nrow = data->mats[i]->nrow;
      dataM1->mats[i]->ncol = data->mats[i]->ncol;
      iaux = dataM1->mats[i]->nrow * dataM1->mats[i]->ncol;
      dataM1->mats[i]->x = Calloc(iaux, double);
      for(j=0; j<iaux; j++) {
        dataM1->mats[i]->x[j] = data->mats[i]->x[j];
      }
    }
    if(debug)
      printf("done \n");
  }

  if(nm2>0) { // copy mats for M2
    if(debug)
      printf("copy mats for M2 ...");
    for(i=0; i<dataM2->n_mats; i++) {
      if(debug>1)
        printf(" %d, ", i);
      dataM2->mats[i] = Calloc(1, inla_cgeneric_mat_tp);
      dataM2->mats[i]->name = data->mats[nm1+i]->name;
      dataM2->mats[i]->nrow = data->mats[nm1+i]->nrow;
      dataM2->mats[i]->ncol = data->mats[nm1+i]->ncol;
      iaux = dataM2->mats[i]->nrow * dataM2->mats[i]->ncol;
      dataM2->mats[i]->x = Calloc(iaux, double);
      for(j=0; j<iaux; j++) {
        dataM2->mats[i]->x[j] = data->mats[nm1+i]->x[j];
      }
    }
    if(debug)
      printf("done \n");
  }

  if(nsm1>0) { // copy smats for M1
    if(debug)
      printf("copy smats for M1 ...");
    for(i=0; i<dataM1->n_smats; i++) {
      if(debug>1)
        printf(" %d, ", i);
      dataM1->smats[i] = Calloc(1, inla_cgeneric_smat_tp);
      dataM1->smats[i]->name = data->smats[i]->name;
      dataM1->smats[i]->nrow = data->smats[i]->nrow;
      dataM1->smats[i]->ncol = data->smats[i]->ncol;
      iaux = data->smats[i]->n;
      dataM1->smats[i]->n = iaux;
      dataM1->smats[i]->i = Calloc(iaux, int);
      dataM1->smats[i]->j = Calloc(iaux, int);
      dataM1->smats[i]->x = Calloc(iaux, double);
      for(j=0; j<iaux; j++) {
        dataM1->smats[i]->x[j] = data->smats[i]->x[j];
      }
    }
    if(debug)
      printf("done \n");
  }

  if(nsm2>0) { // copy smats for M2
    if(debug)
      printf("copy smats for M2 ...");
    for(i=0; i<dataM2->n_smats; i++) {
      if(debug>1)
        printf(" %d, ", i);
      dataM2->smats[i] = Calloc(1, inla_cgeneric_smat_tp);
      dataM2->smats[i]->name = data->smats[nsm1+i]->name;
      dataM2->smats[i]->nrow = data->smats[nsm1+i]->nrow;
      dataM2->smats[i]->ncol = data->smats[nsm1+i]->ncol;
      iaux = data->smats[nsm1+i]->n;
      dataM2->smats[i]->n = iaux;
      dataM2->smats[i]->i = Calloc(iaux, int);
      dataM2->smats[i]->j = Calloc(iaux, int);
      dataM2->smats[i]->x = Calloc(iaux, double);
      for(j=0; j<iaux; j++) {
        dataM2->smats[i]->x[j] = data->smats[nsm1+i]->x[j];
      }
    }
    if(debug)
      printf("done \n");
  }

  // load libs
  static int ltdl_init = 1;
  if (ltdl_init) {
    lt_dlinit();
    lt_dlerror();
  }

  lt_dlhandle handle1;
  lt_dlhandle handle2;

  if(debug)
    printf("open %s \n", &data->chars[3]->chars[0]);
  handle1 = lt_dlopen(&data->chars[3]->chars[0]);
  if(debug)
    printf("open %s \n", &data->chars[5]->chars[0]);
  handle2 = lt_dlopen(&data->chars[5]->chars[0]);

  // model functions
  inla_cgeneric_func_tp *model1_func = NULL;
  inla_cgeneric_func_tp *model2_func = NULL;
  model1_func = (inla_cgeneric_func_tp *) lt_dlsym(handle1, &data->chars[2]->chars[0]);
  model2_func = (inla_cgeneric_func_tp *) lt_dlsym(handle2, &data->chars[4]->chars[0]);

  int nth1 = (int)model1_func(INLA_CGENERIC_INITIAL, NULL, dataM1)[0];
  double *theta1 = &theta[0];
  double *theta2 = &theta[nth1];

  switch (cmd) {

  case INLA_CGENERIC_VOID:
  {
    assert(!(cmd == INLA_CGENERIC_VOID));
    break;
  }

  case INLA_CGENERIC_GRAPH:
  {

    assert(M == data->smats[nsm1+nsm2]->n);

    ret = Calloc(2 + 2 * M, double);
    ret[0] = N;
    ret[1] = M;

    // collect i
    for(i = 0; i<M; i++) {
      ret[2+i] = data->smats[nsm1+nsm2]->i[i];
    }
    // collect j
    for(i = 0; i<M; i++) {
      ret[2+M+i] = data->smats[nsm1+nsm2]->j[i];
    }

    break;
  }

  case INLA_CGENERIC_Q:
  {
    ret = Calloc(2 + M, double);
    ret[0] = -1;		/* REQUIRED */
    ret[1] = M;

    ret1 = model1_func(INLA_CGENERIC_Q, theta1, dataM1);
    ret2 = model2_func(INLA_CGENERIC_Q, theta2, dataM2);

    int m2e = M2-n2;
    int M2e = M2 + m2e;
    if(debug){
      printf("%d, m2e %d, M2e %d\n", data->ints[ni1+ni2-2]->len, m2e, M2e);
    }
    assert(m2e == data->ints[ni1+ni2-1]->len);
    assert((M1 * M2e) == data->smats[nsm1+nsm2]->n);
    double ret2e[M2e];
    double retE[M1 * M2e];
    double daux;

    for(k = 0; k<M2; k++) {
      ret2e[k] = ret2[2+k];
    }
    if(m2e>0){
      for(i=0; i<m2e; i++) {
        ret2e[k++] = ret2[2+data->ints[ni1+ni2-2]->ints[i]];
      }
    }

    k=0;
    for(i=0; i<M1; i++) {
      daux = ret1[2+i];
      for(j=0; j<M2e; j++) {
        retE[k++] = daux * ret2e[j];
      }
    }

    int ox;
    for(k=0; k<data->smats[nsm1+nsm2]->n; k++) {
      ox = (int)data->smats[nsm1+nsm2]->x[k];
      ret[2+k] = retE[ox];
    }

    break;
  }

  case INLA_CGENERIC_MU:
  {
    // return (N, mu)
    // if N==0 then mu is not needed as its taken to be mu[]==0
    ret = Calloc(1, double);
    ret[0] = 0;
    break;
  }

  case INLA_CGENERIC_INITIAL:
  {
    // return c(M, initials)
    // where M is the number of hyperparameters

    ret1 = model1_func(INLA_CGENERIC_INITIAL, NULL, dataM1);
    ret2 = model2_func(INLA_CGENERIC_INITIAL, NULL, dataM2);

    int nth1 = (int)ret1[0], nth2 = (int)ret2[0];

    ret = Calloc(1 +  nth1 + nth2, double);
    ret[0] = nth1 + nth2;

    for(i=0; i<nth1; i++) {
      ret[1 + i] = ret1[1 + i];
    }

    for(i=0; i<nth2; i++) {
      ret[1 + nth1 + i] = ret2[1 + i];
    }

    break;
  }

  case INLA_CGENERIC_LOG_NORM_CONST:
  {
    break;
  }

  case INLA_CGENERIC_LOG_PRIOR:
  {
    // return c(LOG_PRIOR)
    ret1 = model1_func(INLA_CGENERIC_LOG_PRIOR, theta1, dataM1);
    ret2 = model2_func(INLA_CGENERIC_LOG_PRIOR, theta2, dataM2);

    ret = Calloc(1, double);
    ret[0] = ret1[0] + ret2[0];
    break;
  }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  free(ret1);
  free(ret2);

  return (ret);
}
