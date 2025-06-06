
/* cgeneric_get.c
 *
 * Copyright (C) 2024 Elias T Krainski
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
#include <R_ext/Utils.h>				       // needed to allow user interrupts
#include <Rdefines.h>
#include <Rinternals.h>
#include "graphpcor.h"

SEXP cgeneric_element_get(SEXP Rcmd, SEXP Stheta, SEXP Sntheta, SEXP ints, SEXP doubles, SEXP chars, SEXP mats, SEXP smats)
{

	int ni = 0, nd = 0, nc = 0, nm = 0, nsm = 0, nout = 0;
	int i, j, debug;
	int *ntheta = INTEGER(Sntheta);
	int nth = length(Stheta);
	if (ntheta[0] > 0) {
		nth /= ntheta[0];
	}

	// initial check: mandatory arguments
	if (isNewList(ints)) {
		ni = length(ints);
		assert(ni > 1);
		if (ni < 2) {
			error("at least length 2 'ints' must be provided!");
		}
	} else {
		error("'ints' must be a list!");
	}

	if (isNewList(chars)) {
		nc = length(chars);
		assert(nc > 1);
		if (nc < 2) {
			error("at least length 2 'chars' must be provided!");
		}
	} else {
		error("'chars' must be a list!");
	}

	// check the optional arguments
	if (isNull(doubles)) {
		nd = 0;
	} else {
		if (isNewList(doubles)) {
			nd = length(doubles);
			assert(nd > 0);
		} else {
			error("'doubles' must be a list!");
		}
	}
	if (isNull(mats)) {
		nm = 0;
	} else {
		if (isNewList(mats)) {
			nm = length(mats);
			assert(nm > 0);
		} else {
			error("'mats' must be a list");
		}
	}
	if (isNull(smats)) {
		nsm = 0;
	} else {
		if (isNewList(smats)) {
			nsm = length(smats);
			assert(nsm > 0);
		} else {
			error("'smats' must be a list");
		}
	}

	// get initial info
	char *CMD = (char *) CHAR(STRING_ELT(Rcmd, 0));
	// n = asInteger(VECTOR_ELT(ints, 0));
	debug = asInteger(VECTOR_ELT(ints, 1));
	if (debug > 0) {
		Rprintf("Rcmd is %s, debug = %d\n", CMD, debug);
		// Rprintf("n = %d, debug = %d\n", n, debug);
		Rprintf("ni = %d, nd = %d, nc = %d, nm = %d, nsm = %d, ntheta = %d, nth = %d\n", ni, nd, nc, nm, nsm, ntheta[0], nth);
	}

	double *theta = NULL;
	if (!isNull(Stheta)) {
		theta = REAL(Stheta);
		if (debug) {
			Rprintf("theta: ");
			if ((ntheta[0]) < 2) {
				for (i = 0; i < nth; i++) {
					Rprintf("%f ", theta[i]);
				}
				Rprintf("\n");
			} else {
				for (j = 0; j < (ntheta[0]); j++) {
					for (i = 0; i < nth; i++) {
						Rprintf("%f ", theta[i]);
					}
					Rprintf("\n");
				}
			}
		}
	}

	// define objects
	int ilen[ni], clen[nc];
	inla_cgeneric_data_tp *cgeneric_data = Calloc(1, inla_cgeneric_data_tp);
	char *cgeneric_shlib;
	char *cgeneric_model;
	int *iaux;
	double *daux, *ret = NULL;
	char *pcaux;
	const char *caux;

	// collect data from ints
	SEXP inames = getAttrib(ints, R_NamesSymbol);
	cgeneric_data->n_ints = ni;
	cgeneric_data->ints = Calloc(ni, inla_cgeneric_vec_tp *);
	for (i = 0; i < ni; i++) {
		ilen[i] = length(VECTOR_ELT(ints, i));
		cgeneric_data->ints[i] = Calloc(1, inla_cgeneric_vec_tp);
		pcaux = (char *) CHAR(STRING_ELT(inames, i));
		cgeneric_data->ints[i]->name = pcaux;
		cgeneric_data->ints[i]->len = ilen[i];
		cgeneric_data->ints[i]->ints = Calloc(ilen[i], int);
		if (debug > 0) {
			caux = CHAR(STRING_ELT(inames, i));
			Rprintf("length(ints[[%d]]), %s, is %d\n", i + 1, caux, ilen[i]);
		}
		iaux = INTEGER(VECTOR_ELT(ints, i));
		for (j = 0; j < ilen[i]; j++) {
			cgeneric_data->ints[i]->ints[j] = iaux[j];
		}
	}

	if (nd > 0) {
		// collect lengths and names from doubles
		int dlen[nd];
		SEXP dnames = getAttrib(doubles, R_NamesSymbol);
		for (i = 0; i < nd; i++) {
			dlen[i] = length(VECTOR_ELT(doubles, i));
			if (debug > 0) {
				caux = CHAR(STRING_ELT(dnames, i));
				Rprintf("length(doubles[[%d]]), %s, is %d\n", i + 1, caux, dlen[i]);
			}
		}
		// allocate and collect doubles
		cgeneric_data->n_doubles = nd;
		cgeneric_data->doubles = Calloc(nd, inla_cgeneric_vec_tp *);
		pcaux = (char *) CHAR(STRING_ELT(dnames, 0));
		for (i = 0; i < nd; i++) {
			cgeneric_data->doubles[i] = Calloc(1, inla_cgeneric_vec_tp);
			pcaux = (char *) CHAR(STRING_ELT(dnames, i));
			cgeneric_data->doubles[i]->name = pcaux;
			cgeneric_data->doubles[i]->len = dlen[i];
			cgeneric_data->doubles[i]->doubles = Calloc(dlen[i], double);
			daux = REAL(VECTOR_ELT(doubles, i));
			for (j = 0; j < dlen[i]; j++) {
				cgeneric_data->doubles[i]->doubles[j] = daux[j];
			}
		}
	}

	// collect data from chars
	SEXP cnames = getAttrib(chars, R_NamesSymbol);
	cgeneric_data->n_chars = nc;
	cgeneric_data->chars = Calloc(nc, inla_cgeneric_vec_tp *);
	for (i = 0; i < nc; i++) {
		clen[i] = length(VECTOR_ELT(chars, i));
		caux = CHAR(STRING_ELT(cnames, i));
		if (debug > 0) {
			Rprintf("length(chars[[%d]]), %s, is %d\n", i + 1, caux, clen[i]);
		}
		assert(clen[i] == 1);
		cgeneric_data->chars[i] = Calloc(1, inla_cgeneric_vec_tp);
		pcaux = (char *) CHAR(STRING_ELT(cnames, i));
		cgeneric_data->chars[i]->name = pcaux;
		cgeneric_data->chars[i]->len = clen[i];
		cgeneric_data->chars[i]->chars = Calloc(clen[i] + 1L, char);
		if (debug > 0) {
			Rprintf("%d: length(%s) is %d\n", i + 1, cgeneric_data->chars[i]->name, clen[i]);
			Rprintf("%s, ", pcaux);
			Rprintf("%s :", CHAR(STRING_ELT(VECTOR_ELT(chars, i), 0)));
		}
		pcaux = (char *) CHAR(STRING_ELT(VECTOR_ELT(chars, i), 0));
		cgeneric_data->chars[i]->chars = pcaux;
		if (debug > 0) {
			Rprintf("%s \n", cgeneric_data->chars[i]->chars);
		}
	}

	// check the mandatory strings
	//  assert(cgeneric_data->chars[0]->name == "model");
	if (strcmp(cgeneric_data->chars[0]->name, "model") != 0)
		error("'chars[[1]]' name is not equal 'model'");
	cgeneric_model = cgeneric_data->chars[0]->chars;
	// assert(cgeneric_data->chars[1]->name == "shlib");
	if (strcmp(cgeneric_data->chars[1]->name, "shlib") != 0)
		error("'chars[[2]]' name is not equal 'shlib'");
	cgeneric_shlib = cgeneric_data->chars[1]->chars;

	if (nm > 0) {
		// collect data from mats
		int mnr[nm], mnc[nm], mlen[nm];
		SEXP mnames = getAttrib(mats, R_NamesSymbol);
		cgeneric_data->n_mats = nm;
		cgeneric_data->mats = Calloc(nm, inla_cgeneric_mat_tp *);
		pcaux = (char *) CHAR(STRING_ELT(mnames, 0));
		for (i = 0; i < nm; i++) {
			daux = REAL(VECTOR_ELT(mats, i));
			mnr[i] = daux[0];
			mnc[i] = daux[1];
			mlen[i] = mnr[i] * mnc[i];
			if (debug > 0) {
				caux = CHAR(STRING_ELT(mnames, i));
				Rprintf("dim(mats[[%d]]), %s, is %d %d\n", i + 1, caux, mnr[i], mnc[i]);
			}
			cgeneric_data->mats[i] = Calloc(1, inla_cgeneric_mat_tp);
			pcaux = (char *) CHAR(STRING_ELT(mnames, i));
			cgeneric_data->mats[i]->name = pcaux;
			cgeneric_data->mats[i]->nrow = mnr[i];
			cgeneric_data->mats[i]->ncol = mnc[i];
			cgeneric_data->mats[i]->x = Calloc(mlen[i], double);
			for (j = 0; j < mlen[i]; j++) {
				cgeneric_data->mats[i]->x[j] = daux[2 + j];
			}
		}

	}

	if (nsm > 0) {
		// collect lengths and names from smats
		int smlen[nsm], smnr[nsm], smnc[nsm], smn[nsm];
		SEXP smnames = getAttrib(smats, R_NamesSymbol);
		cgeneric_data->n_smats = nsm;
		cgeneric_data->smats = Calloc(nsm, inla_cgeneric_smat_tp *);
		pcaux = (char *) CHAR(STRING_ELT(smnames, 0));
		for (i = 0; i < nsm; i++) {
			smlen[i] = length(VECTOR_ELT(smats, i));
			if (debug > 0) {
				caux = CHAR(STRING_ELT(smnames, i));
				Rprintf("length(smats[[%d]]), %s, is %d\n", i + 1, caux, smlen[i]);
			}
			daux = REAL(VECTOR_ELT(smats, i));
			int offset = 0;
			smnr[i] = daux[offset++];
			smnc[i] = daux[offset++];
			smn[i] = daux[offset++];
			cgeneric_data->smats[i] = Calloc(1, inla_cgeneric_smat_tp);
			pcaux = (char *) CHAR(STRING_ELT(smnames, i));
			cgeneric_data->smats[i]->name = pcaux;
			cgeneric_data->smats[i]->nrow = smnr[i];
			cgeneric_data->smats[i]->ncol = smnc[i];
			cgeneric_data->smats[i]->n = smn[i];
			cgeneric_data->smats[i]->i = Calloc(smn[i], int);
			cgeneric_data->smats[i]->j = Calloc(smn[i], int);
			cgeneric_data->smats[i]->x = Calloc(smn[i], double);
			for (j = 0; j < smn[i]; j++) {
				cgeneric_data->smats[i]->i[j] = (int) daux[offset++];
			}
			for (j = 0; j < smn[i]; j++) {
				cgeneric_data->smats[i]->j[j] = (int) daux[offset++];
			}
			for (j = 0; j < smn[i]; j++) {
				cgeneric_data->smats[i]->x[j] = daux[offset++];
			}
		}

	}

	// load lib
	static int ltdl_init = 1;

	if (ltdl_init) {
#pragma omp critical (Name_97a987bafe2f18c620f3cefa664a9f775cb7559b)
		if (ltdl_init) {
			lt_dlinit();
			lt_dlerror();
			ltdl_init = 0;
		}
	}

	lt_dlhandle handle;
	inla_cgeneric_func_tp *model_func = NULL;

	handle = lt_dlopen(cgeneric_shlib);
	model_func = (inla_cgeneric_func_tp *) lt_dlsym(handle, cgeneric_model);
	assert(model_func);

	SEXP Rret = R_NilValue;

	if (strcmp(CMD, "graph") == 0) {
		ret = model_func(INLA_CGENERIC_GRAPH, theta, cgeneric_data);
		nout = (int) ret[1];
		SEXP ii = PROTECT(allocVector(INTSXP, nout));
		iaux = INTEGER(ii);
		for (i = 0; i < nout; i++) {
			iaux[i] = ret[2 + i];
		}
		SEXP jj = PROTECT(allocVector(INTSXP, nout));
		iaux = INTEGER(jj);
		for (i = 0; i < nout; i++) {
			iaux[i] = ret[2 + nout + i];
		}
		Rret = PROTECT(allocVector(VECSXP, 2));
		SET_VECTOR_ELT(Rret, 0, ii);
		SET_VECTOR_ELT(Rret, 1, jj);
		UNPROTECT(3);
		if (debug > 0) {
			Rprintf("graph with n = %d and %d nz\n", (int) ret[0], nout);
		}
	}

	if (strcmp(CMD, "Q") == 0) {
		ret = model_func(INLA_CGENERIC_Q, theta, cgeneric_data);
		nout = (int) ret[1];
		Rret = PROTECT(allocVector(REALSXP, nout));
		daux = REAL(Rret);
		for (i = 0; i < nout; i++) {
			daux[i] = ret[2 + i];
		}
		UNPROTECT(1);
		if (debug > 0) {
			Rprintf("Q with %d nz\n", nout);
		}
	}

	if (strcmp(CMD, "mu") == 0) {
		ret = model_func(INLA_CGENERIC_MU, theta, cgeneric_data);
		Rret = PROTECT(allocVector(REALSXP, 1));
		REAL(Rret)[0] = ret[0];
		UNPROTECT(1);
	}

	if (strcmp(CMD, "initial") == 0) {
		ret = model_func(INLA_CGENERIC_INITIAL, theta, cgeneric_data);
		nout = (int) ret[0];
		if (debug > 0) {
			Rprintf("intial: %d elements (%f)\n", nout, ret[0]);
			if (nout > 0) {
				for (i = 0; i < nout; i++)
					Rprintf("theta[%d] = %f\n", i, ret[1 + i]);
			}
		}
		if (nout > 0) {
			Rret = PROTECT(allocVector(REALSXP, nout));
			daux = REAL(Rret);
			for (i = 0; i < nout; i++) {
				daux[i] = ret[1 + i];
			}
			UNPROTECT(1);
		}
	}

	if (strcmp(CMD, "log_prior") == 0) {
		Rret = PROTECT(allocVector(REALSXP, ntheta[0]));
		for (j = 0; j < (ntheta[0]); j++) {
			ret = model_func(INLA_CGENERIC_LOG_PRIOR, &theta[j * nth], cgeneric_data);
			REAL(Rret)[j] = ret[0];
		}
		UNPROTECT(1);
	}

	free(cgeneric_data);
	free(ret);

	return Rret;
}
