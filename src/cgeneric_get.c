
/* cgeneric_get.c
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
 *        Elias Krainski
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 */

#include <ltdl.h>
#include <R.h>
#include <R_ext/Utils.h> // needed to allow user interrupts
#include <Rdefines.h>
#include <Rinternals.h>
#include "cgeneric_defs.h"

SEXP cgeneric_element_get(SEXP Rcmd, SEXP Stheta, SEXP ints, SEXP doubles, SEXP chars) {

  // initial check
  if(!isNewList(ints))
    error("'ints' must be a list");
  if(!isNewList(doubles))
    error("'doubles' must be a list");
  if(!isNewList(chars))
    error("'chars' must be a list");

  // get initial info
  char *CMD = (char*)CHAR(STRING_ELT(Rcmd, 0));
  int i, n, debug;
  n = asInteger(VECTOR_ELT(ints, 0));
  debug = asInteger(VECTOR_ELT(ints, 1));
  if(debug>0) {
    Rprintf("Rcmd is %s\n", CMD);
    Rprintf("n = %d, debug = %d\n", n, debug);
  }
  int ni = length(ints);
  int nd = length(doubles);
  int nc = length(chars);
  if(debug>0) {
    Rprintf("ni = %d, ", ni);
    Rprintf("nd = %d, ", nd);
    Rprintf("nc = %d \n", nc);
  }
  assert(ni>0);

  double *theta = NULL;
  if(!isNull(Stheta))
    theta = REAL(Stheta);
  if(debug) {
    Rprintf("theta: ");
    for(i=0; i<length(Stheta); i++) {
      Rprintf("%f ", theta[i]);
    }
    Rprintf("\n");
  }

  int j, ilen[ni], dlen[nd], clen[nc];

  // collect data lengths and names
  const char *caux;
  SEXP inames = getAttrib(ints, R_NamesSymbol);
  SEXP dnames = getAttrib(doubles, R_NamesSymbol);
  SEXP cnames = getAttrib(chars, R_NamesSymbol);
  for(i=0; i<ni; i++) {
    ilen[i] = length(VECTOR_ELT(ints, i));
    if(debug>0) {
      caux = CHAR(STRING_ELT(inames, i));
      Rprintf("length(ints[[%d]]), %s, is %d\n",
              i+1, caux, ilen[i]);
    }
  }
  if(nd>0) {
    for(i=0; i<nd; i++) {
      dlen[i] = length(VECTOR_ELT(doubles, i));
      if(debug>0) {
        caux = CHAR(STRING_ELT(dnames, i));
        Rprintf("length(doubles[[%d]]), %s, is %d\n", i+1, caux, dlen[i]);
      }
    }
  }
  if(nc>0) {
    for(i=0; i<nc; i++) {
      clen[i] = length(VECTOR_ELT(chars, i));
      caux = CHAR(STRING_ELT(cnames, i));
      if(debug>0) {
        Rprintf("length(chars[[%d]]), %s, is %d\n", i+1, caux, clen[i]);
      }
      assert(clen[i]==1);
    }
  }

// define objects
  lt_dlhandle handle;
  inla_cgeneric_func_tp *model_func = NULL;
  inla_cgeneric_data_tp *cgeneric_data = Calloc(1, inla_cgeneric_data_tp);
  char *cgeneric_shlib;
  char *cgeneric_model;
  int naux, *iaux;
  double *daux, *ret = NULL;
  char *pcaux;

  // allocate and collect ints
  cgeneric_data->n_ints = ni;
  cgeneric_data->ints = Calloc(ni, inla_cgeneric_vec_tp *);
  for(i=0; i<ni; i++) {
    cgeneric_data->ints[i] = Calloc(1, inla_cgeneric_vec_tp);
    pcaux = (char*) CHAR(STRING_ELT(inames,i));
    cgeneric_data->ints[i]->name = pcaux;
    cgeneric_data->ints[i]->len = ilen[i];
    cgeneric_data->ints[i]->ints = Calloc(ilen[i], int);
    if(debug>0) {
      Rprintf("%d: length(%s) is %d\n", i+1, cgeneric_data->ints[i]->name, ilen[i]);
    }
    iaux = INTEGER(VECTOR_ELT(ints, i));
    for(j=0; j<ilen[i]; j++) {
      cgeneric_data->ints[i]->ints[j] = iaux[j];
//      Rprintf("%d ", cgeneric_data->ints[i]->ints[j]);
    }
  //  Rprintf("\n");
  }

  // allocate and collect doubles
  cgeneric_data->n_doubles = nd;
  cgeneric_data->doubles = Calloc(nd, inla_cgeneric_vec_tp *);
  pcaux = (char*) CHAR(STRING_ELT(dnames,0));
  for(i=0; i<nd; i++) {
    cgeneric_data->doubles[i] = Calloc(1, inla_cgeneric_vec_tp);
    pcaux = (char*) CHAR(STRING_ELT(dnames,i));
    cgeneric_data->doubles[i]->name = pcaux;
    cgeneric_data->doubles[i]->len = dlen[i];
    cgeneric_data->doubles[i]->doubles = Calloc(dlen[i], double);
    if(debug>0) {
      Rprintf("%d: length(%s) is %d\n", i+1, cgeneric_data->doubles[i]->name, dlen[i]);
    }
    daux = REAL(VECTOR_ELT(doubles, i));
    for(j=0; j<dlen[i]; j++) {
      cgeneric_data->doubles[i]->doubles[j] = daux[j];
    //  Rprintf("%f ", cgeneric_data->doubles[i]->doubles[j]);
    }
    //Rprintf("\n");
  }

  // allocate and collect chars
  cgeneric_data->n_chars = nc;
  cgeneric_data->chars = Calloc(nc, inla_cgeneric_vec_tp *);
  for(i=0; i<nc; i++) {
    cgeneric_data->chars[i] = Calloc(1, inla_cgeneric_vec_tp);
    pcaux = (char*) CHAR(STRING_ELT(cnames,i));
    cgeneric_data->chars[i]->name = pcaux;
    cgeneric_data->chars[i]->len = clen[i];
    cgeneric_data->chars[i]->chars = Calloc(clen[i] + 1L, char);
    if(debug>0) {
      Rprintf("%d: length(%s) is %d\n", i+1, cgeneric_data->chars[i]->name, clen[i]);
    }
    if(debug>0) {
      Rprintf("%s, ", pcaux);
      Rprintf("%s :", CHAR(STRING_ELT(VECTOR_ELT(chars, i), 0)));
    }
    pcaux = (char*)CHAR(STRING_ELT(VECTOR_ELT(chars, i), 0));
    cgeneric_data->chars[i]->chars = pcaux;
    if(debug>0) {
      Rprintf("%s \n", cgeneric_data->chars[i]->chars);
    }
  }
  // check strings
  assert(cgeneric_data->chars[0]->name == "model");
  cgeneric_model = cgeneric_data->chars[0]->chars;
  assert(cgeneric_data->chars[1]->name == "shlib");
  cgeneric_shlib = cgeneric_data->chars[1]->chars;

  // load lib
  static int ltdl_init = 1;
  if (ltdl_init) {
    lt_dlinit();
    lt_dlerror();
  }
  handle = lt_dlopen(cgeneric_shlib);
  model_func = (inla_cgeneric_func_tp *) lt_dlsym(handle, cgeneric_model);

  int nout = 0;
  SEXP Rret = R_NilValue;

  if(strcmp(CMD, "graph") == 0) {
    ret = model_func(INLA_CGENERIC_GRAPH, theta, cgeneric_data);
    nout = (int)ret[1];
    SEXP ii = PROTECT(allocVector(INTSXP, nout));
    iaux = INTEGER(ii);
    for(i = 0; i < nout; i++) {
      iaux[i] = ret[i];
    }
    SEXP jj = PROTECT(allocVector(INTSXP, nout));
    iaux = INTEGER(jj);
    for(i = 0; i < nout; i++) {
      iaux[i] = ret[nout+i];
    }
    Rret = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(Rret, 0, ii);
    SET_VECTOR_ELT(Rret, 1, jj);
    UNPROTECT(3);
  }

  if(strcmp(CMD, "Q") == 0) {
    ret = model_func(INLA_CGENERIC_Q, theta, cgeneric_data);
    nout = (int)ret[1];
    Rret = PROTECT(allocVector(REALSXP, nout));
    daux = REAL(Rret);
    for(i = 0; i < nout; i++) {
      daux[i] = ret[2+i];
    }
    UNPROTECT(1);
  }

  if(strcmp(CMD, "mu") == 0) {
    ret = model_func(INLA_CGENERIC_MU, theta, cgeneric_data);
    Rret = PROTECT(allocVector(REALSXP, 1));
    REAL(Rret)[0] = ret[0];
    UNPROTECT(1);
  }

  if(strcmp(CMD, "initial") == 0) {
    ret = model_func(INLA_CGENERIC_INITIAL, theta, cgeneric_data);
    nout = (int)ret[0];
    Rret = PROTECT(allocVector(REALSXP, nout));
    daux = REAL(Rret);
    for(i = 0; i < nout; i++) {
      daux[i] = ret[1 + i];
    }
    UNPROTECT(1);
  }

  if(strcmp(CMD, "log_prior") == 0) {
    ret = model_func(INLA_CGENERIC_LOG_PRIOR, theta, cgeneric_data);
    Rret = PROTECT(allocVector(REALSXP, 1));
    REAL(Rret)[0] = ret[0];
    UNPROTECT(1);
  }

  Free(cgeneric_data);
  Free(ret);

  return Rret;
}
