
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

  double *ret = NULL;

  //  Q = Q1 (x) Q2

  // cmd is length 1 string, however,
  // theta contains theta[1:nth1-1, nth1:(nth1+nth2-1)]
  //   nth1 is #theta for M1, nth2 is #theta for M2
  // data contains data for model M1 and model M2.
  // data->ints[0]->ints contain
  //     {n,
  //      n1, ni1, nd1, nc1, nm1, nsm1,
  //      n2, ni2, nd2, nc2, nm2, nsm2}
  // so that
  //  data->ints[1:ni1-1] contain ints for M1
  //  data->ints[ni1+1:ni2-1] contain ints for M2
  //  data->doubles[1:nd1] contain doubles for M1
  //  data->doubles[nd1+1:nd2-1] contain doubles for M2
  //  data->chars[1:nc1-1] contain chars for M1
  //  data->chars[nc1+1:nc2-1] contain chars for M2
  //  data->mats[1:nm1-1] contain mats for M1
  //  data->mats[nm1+1:nm2-1] contain mats for M2
  //  data->smats[1:nsm1-1] contain smats for M1
  //  data->smats[nsm1+1:nsm2-1] contain smats for M2

  int i, j, n1, n2, iaux;
  int ni1, nd1, nc1, nm1, nsm1;
  int ni2, nd2, nc2, nm2, nsm2;
  int N, M, M1, M2, k;

  n1 = data->ints[0]->ints[1];
  ni1 = data->ints[0]->ints[2];
  nd1 = data->ints[0]->ints[3];
  nc1 = data->ints[0]->ints[4];
  nm1 = data->ints[0]->ints[5];
  nsm1 = data->ints[0]->ints[6];

  n2 = data->ints[0]->ints[7];
  ni2 = data->ints[0]->ints[8];
  nd2 = data->ints[0]->ints[9];
  nc2 = data->ints[0]->ints[10];
  nm2 = data->ints[0]->ints[11];
  nsm2 = data->ints[0]->ints[12];

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
  dataM1->mats = Calloc(dataM1->n_mats, inla_cgeneric_mat_tp *);
  dataM2->mats = Calloc(dataM2->n_mats, inla_cgeneric_mat_tp *);
  dataM1->smats = Calloc(dataM1->n_smats, inla_cgeneric_smat_tp *);
  dataM2->smats = Calloc(dataM2->n_smats, inla_cgeneric_smat_tp *);

  // copy ints for M1
  for(i=0; i<dataM1->n_ints; i++) {
    dataM1->ints[i] = Calloc(1, inla_cgeneric_vec_tp);
    dataM1->ints[i]->name = data->ints[i]->name;
    dataM1->ints[i]->len = data->ints[i]->len;
    dataM1->ints[i]->ints = Calloc(dataM1->ints[i]->len, int);
    for(j=0; j<dataM1->ints[i]->len; j++) {
      dataM1->ints[i]->ints[j] = data->ints[i]->ints[j];
    }
  }

  // copy ints for M2
  for(i=0; i<dataM2->n_ints; i++) {
    dataM2->ints[i] = Calloc(1, inla_cgeneric_vec_tp);
    dataM2->ints[i]->name = data->ints[ni1+i]->name;
    dataM2->ints[i]->len = data->ints[ni1+i]->len;
    dataM2->ints[i]->ints = Calloc(dataM2->ints[i]->len, int);
    for(j=0; j<dataM2->ints[i]->len; j++) {
      dataM2->ints[i]->ints[j] = data->ints[ni1+i]->ints[j];
    }
  }

  if(nd1>0) { // copy doubles for M1
    for(i=0; i<dataM1->n_doubles; i++) {
      dataM1->doubles[i] = Calloc(1, inla_cgeneric_vec_tp);
      dataM1->doubles[i]->name = data->doubles[i]->name;
      dataM1->doubles[i]->len = data->doubles[i]->len;
      dataM1->doubles[i]->doubles = Calloc(dataM1->doubles[i]->len, double);
      for(j=0; j<dataM1->doubles[i]->len; j++) {
        dataM1->doubles[i]->doubles[j] = data->doubles[i]->doubles[j];
      }
    }
  }

  if(nd2>0) { // copy doubles for M2
    for(i=0; i<dataM2->n_doubles; i++) {
      dataM2->doubles[i] = Calloc(1, inla_cgeneric_vec_tp);
      dataM2->doubles[i]->name = data->doubles[nd1+i]->name;
      dataM2->doubles[i]->len = data->doubles[nd2+i]->len;
      dataM2->doubles[i]->doubles = Calloc(dataM2->doubles[i]->len, double);
      for(j=0; j<dataM2->doubles[i]->len; j++) {
        dataM2->doubles[i]->doubles[j] = data->doubles[nd1+i]->doubles[j];
      }
    }
  }

  if(nc1>0) { // copy chars for M1
    for(i=0; i<dataM1->n_chars; i++) {
      dataM1->chars[i] = Calloc(1, inla_cgeneric_vec_tp);
      dataM1->chars[i]->name = data->chars[i]->name;
      dataM1->chars[i]->len = data->chars[i]->len;
      dataM1->chars[i]->chars = Calloc(dataM1->chars[i]->len, char);
      for(j=0; j<dataM1->chars[i]->len; j++) {
        dataM1->chars[i]->chars[j] = data->chars[i]->chars[j];
      }
    }
  }

  if(nc2>0) { // copy chars for M2
    for(i=0; i<dataM2->n_chars; i++) {
      dataM2->chars[i] = Calloc(1, inla_cgeneric_vec_tp);
      dataM2->chars[i]->name = data->chars[nd1+i]->name;
      dataM2->chars[i]->len = data->chars[nd2+i]->len;
      dataM2->chars[i]->chars = Calloc(dataM2->chars[i]->len, char);
      for(j=0; j<dataM2->chars[i]->len; j++) {
        dataM2->chars[i]->chars[j] = data->chars[nd1+i]->chars[j];
      }
    }
  }

  if(nm1>0) { // copy mats for M1
    for(i=0; i<dataM1->n_mats; i++) {
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
  }

  if(nm2>0) { // copy mats for M2
    for(i=0; i<dataM2->n_mats; i++) {
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
  }

  if(nsm1>0) { // copy smats for M1
    for(i=0; i<dataM1->n_smats; i++) {
      dataM1->smats[i] = Calloc(1, inla_cgeneric_smat_tp);
      dataM1->smats[i]->name = data->smats[i]->name;
      dataM1->smats[i]->nrow = data->smats[i]->nrow;
      dataM1->smats[i]->ncol = data->smats[i]->ncol;
      iaux = dataM1->smats[i]->nrow * dataM1->smats[i]->ncol;
      dataM1->smats[i]->n = iaux;
      dataM1->smats[i]->i = Calloc(iaux, int);
      dataM1->smats[i]->j = Calloc(iaux, int);
      dataM1->smats[i]->x = Calloc(iaux, double);
      for(j=0; j<iaux; j++) {
        dataM1->smats[i]->x[j] = data->smats[i]->x[j];
      }
    }
  }

  if(nsm2>0) { // copy smats for M2
    for(i=0; i<dataM2->n_smats; i++) {
      dataM2->smats[i] = Calloc(1, inla_cgeneric_smat_tp);
      dataM2->smats[i]->name = data->smats[nsm1+i]->name;
      dataM2->smats[i]->nrow = data->smats[nsm1+i]->nrow;
      dataM2->smats[i]->ncol = data->smats[nsm1+i]->ncol;
      iaux = dataM2->smats[i]->nrow * dataM2->smats[i]->ncol;
      dataM2->smats[i]->n = iaux;
      dataM2->smats[i]->i = Calloc(iaux, int);
      dataM2->smats[i]->j = Calloc(iaux, int);
      dataM2->smats[i]->x = Calloc(iaux, double);
      for(j=0; j<iaux; j++) {
        dataM2->smats[i]->x[j] = data->smats[nsm1+i]->x[j];
      }
    }
  }

  // load lib
  static int ltdl_init = 1;
  if (ltdl_init) {
    lt_dlinit();
    lt_dlerror();
  }
  char *shlib1, *shlib2;
  char *model1, *model2;
  lt_dlhandle handle1;
  lt_dlhandle handle2;
  model1 = dataM1->chars[0]->chars;
  model2 = dataM2->chars[0]->chars;
  shlib1 = dataM1->chars[1]->chars;
  shlib2 = dataM2->chars[1]->chars;
  handle1 = lt_dlopen(shlib1);
  inla_cgeneric_func_tp *model1_func = NULL;
  model1_func = (inla_cgeneric_func_tp *) lt_dlsym(handle1, model1);
  handle2 = lt_dlopen(shlib2);
  inla_cgeneric_func_tp *model2_func = NULL;
  model2_func = (inla_cgeneric_func_tp *) lt_dlsym(handle2, model2);

  double *out1, *ret1 = NULL;
  double *out2, *ret2 = NULL;

  out1 = model1_func(INLA_CGENERIC_GRAPH, NULL, dataM1);
  out2 = model2_func(INLA_CGENERIC_GRAPH, NULL, dataM2);
  M1 = (int)out1[1];
  M2 = (int)out2[1];
  M = M1 * M2;

  switch (cmd) {

  case INLA_CGENERIC_VOID:
  {
    assert(!(cmd == INLA_CGENERIC_VOID));
    break;
  }

  case INLA_CGENERIC_GRAPH:
  {

    ret = Calloc(2 + 2 * M, double);
    ret[0] = N;
    ret[1] = M;

    k = 2;
    int k1 = 2, k2, m2;
    // collect i
    for(i = 0; i<M1; i++) {
      k2 = 2;
      m2 = M1*i;
      for(j = 0; j<M2; j++) {
        ret[k++] = out1[k1] * m2 + out2[k2++];
      }
      k1++;
    }
    // collect j
    k1 = 2 + M1;
    for(i = 0; i<M1; i++) {
      k2 = M2 + 2;
      m2 = M1*i;
      for(j = 0; j<M2; j++) {
        ret[k++] = out1[k1] * m2 + out2[k2++];
      }
      k1++;
    }

    break;
  }

  case INLA_CGENERIC_Q:
  {
    ret = Calloc(2 + M, double);
    ret[0] = -1;		/* REQUIRED */
    ret[1] = M;

    ret1 = model1_func(INLA_CGENERIC_Q, NULL, dataM1);
    ret2 = model2_func(INLA_CGENERIC_Q, NULL, dataM2);

    k = 2;
    int k1 = 2, k2;
    // collect x_1 (x) x_2
    for(i = 0; i<M1; i++) {
      k2 = 2;
      for(j = 0; j<M2; j++) {
        ret[k++] = ret1[k1] * ret2[k2++];
      }
      k1++;
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
    ret1 = model1_func(INLA_CGENERIC_LOG_PRIOR, NULL, dataM1);
    ret2 = model2_func(INLA_CGENERIC_LOG_PRIOR, NULL, dataM2);

    ret = Calloc(1, double);
    ret[0] = ret1[0] + ret2[0];
    break;
  }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return (ret);
}
