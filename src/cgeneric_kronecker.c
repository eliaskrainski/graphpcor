
/* cgeneric_kronecker.c
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
 *        Elias T Krainski
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 */

#include <ltdl.h>
#include "graphpcor.h"

typedef struct
{
	int nth1;
	lt_dlhandle handle1;
	lt_dlhandle handle2;
	inla_cgeneric_func_tp *model1_func;
	inla_cgeneric_func_tp *model2_func;
}
	cache_tp;

double *inla_cgeneric_kronecker(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{
	// concatenated data approach of the lists

	// Merge cgeneric data from model1 (M1) and model2 (M2)
	// to define a Kronecker product model (KM) with precision as
	// Q = Q1 (x) Q2
	// where
	// Q1 is the precision for M1
	// Q2 is the precision for M2

	// cmd: length 1 string
	// theta: {theta1, theta2}
	// data: {data1, data2}, but data1->ints[0]->ints[0...10] contains
	// n1, ni1, nd1, nc1, nm1, nsm1, m1,
	// ni2, nd2, nc2, nm2, nsm2, m2, n, M}
	// n1 is n for mode1 (M1)
	// ni1 is the number of ints for M1
	// ...
	// m1 is M for M1
	// ni2 is the number of ints for M2
	// ...
	// n is n for KM
	// M is M for KM
	// so that
	// data->ints[<ni1] contain ints for M1
	// data->doubles[<nd1] contain doubles for M1
	// data->chars[<nc1] contain chars for M1
	// data->mats[<nm1] contain mats for M1
	// data->smats[<nsm1] contain smats for M1
	// and
	// data->ints[ni1 .. < ni1+ni2] contain ints for M2
	// data->doubles[nd1 .. < nd1+nd2] contain doubles for M2
	// data->chars[nc1 .. < nc1+nc2] contain chars for M2
	// data->mats[nm1 .. < nm1+nm2] contain mats for M2
	// data->smats[nsm1 .. < nsm1+nsm2] contain smats for M2
	// additionally:
	// data->ints[ni1+ni2] contain nu1 index
	// data->ints[ni1+ni2+1] contain nu2 index
	// data->smatrices[nsm1+nsm2] contains the graph
	// where ->x is the order

	double *ret1 = NULL;				       // to store output from M1.
	double *ret2 = NULL;				       // to store output from M2.
	double *ret = NULL;				       // to return;

	int i, j, k;
	int ni1, nd1, nc1, nm1, nsm1;
	int ni2, nsm2; //, nd2, nc2, nm2, nsm2;

	// int n1 = data->ints[0]->ints[0];
	ni1 = data->ints[0]->ints[1];
	nd1 = data->ints[0]->ints[2];
	nc1 = data->ints[0]->ints[3];
	nm1 = data->ints[0]->ints[4];
	nsm1 = data->ints[0]->ints[5];
	int M1 = data->ints[0]->ints[6];

	// int n2 = data->ints[ni1]->ints[0];
	ni2 = data->ints[0]->ints[7];
//	nd2 = data->ints[0]->ints[8];
//	nc2 = data->ints[0]->ints[9];
//	nm2 = data->ints[0]->ints[10];
	nsm2 = data->ints[0]->ints[11];
	int M2 = data->ints[0]->ints[12];
	int n = data->ints[0]->ints[13];
	int M = data->ints[0]->ints[14];

	assert(ni1 > 1);
	assert(ni2 > 1);
	assert(nc1 > 1);
//	assert(nc2 > 1);
	assert(data->n_chars > 5);

	assert(!strcasecmp(data->ints[1]->name, "debug"));     // this will always be the case
	assert(!strcasecmp(data->ints[ni1 + 1]->name, "debug"));	// this will always be the case
//	int debug = (data->ints[1]->ints[0] | data->ints[ni1 + 1]->ints[0]);

	assert(!strcasecmp(data->ints[ni1 + ni2]->name, "idx1u"));	// this will always be the case
	assert(!strcasecmp(data->ints[ni1 + ni2 + 1]->name, "idx2u"));	// this will always be the case

	assert(!strcasecmp(data->ints[nsm1 + nsm2]->name, "Kgraph"));	// this will always be the case

	inla_cgeneric_data_tp *dataM1 = Calloc(1, inla_cgeneric_data_tp);
	inla_cgeneric_data_tp *dataM2 = Calloc(1, inla_cgeneric_data_tp);

	dataM1->ints = &data->ints[0];
	dataM1->doubles = &data->doubles[0];
	dataM1->chars = &data->chars[2];		       // first two is for KM!
	dataM1->mats = &data->mats[0];
	dataM1->smats = &data->smats[0];

	dataM2->ints = &data->ints[ni1];
	dataM2->doubles = &data->doubles[nd1];
	dataM2->chars = &data->chars[2 + nc1];		       // first two is for KM!
	dataM2->mats = &data->mats[nm1];
	dataM2->smats = &data->smats[nsm1];

	if (!(data->cache)) {
#pragma omp critical (Name_5bd4b7198feb5550e84446518f90d47072338c18)
		if (!(data->cache)) {

			lt_dlinit();
			lt_dlerror();
			cache_tp *c = Calloc(1, cache_tp);

			c->handle1 = lt_dlopen(&dataM1->chars[1]->chars[0]);

			if (strcmp(&dataM1->chars[1]->chars[0], &dataM2->chars[1]->chars[0]) != 0) {
				c->handle2 = lt_dlopen(&dataM2->chars[1]->chars[0]);
			} else {
				c->handle2 = c->handle1;
			}

			c->model1_func = (inla_cgeneric_func_tp *) lt_dlsym(c->handle1, &dataM1->chars[0]->chars[0]);
			c->model2_func = (inla_cgeneric_func_tp *) lt_dlsym(c->handle2, &dataM2->chars[0]->chars[0]);
			assert(c->model1_func);
			assert(c->model2_func);
			c->nth1 = (int) c->model1_func(INLA_CGENERIC_INITIAL, NULL, dataM1)[0];

			data->cache = (void *) c;
		}
	}

	assert(data->cache);
	cache_tp *c = (cache_tp *) data->cache;

	double *theta1 = &theta[0];
	double *theta2 = &theta[c->nth1];

	switch (cmd) {
	case INLA_CGENERIC_VOID:
	{
		assert(!(cmd == INLA_CGENERIC_VOID));
		break;
	}

	case INLA_CGENERIC_GRAPH:
	{

		assert(M == data->smats[nsm1 + nsm2]->n);

		ret = Calloc(2 + 2 * M, double);
		ret[0] = n;
		ret[1] = M;

		// collect i
		for (i = 0; i < M; i++) {
			ret[2 + i] = data->smats[nsm1 + nsm2]->i[i];
		}
		// collect j
		for (i = 0; i < M; i++) {
			ret[2 + M + i] = data->smats[nsm1 + nsm2]->j[i];
		}

		break;
	}

	case INLA_CGENERIC_Q:
	{
		ret = Calloc(2 + M, double);
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;

		ret1 = c->model1_func(INLA_CGENERIC_Q, theta1, dataM1);
		ret2 = c->model2_func(INLA_CGENERIC_Q, theta2, dataM2);

		int nu1 = data->ints[ni1 + ni2]->len;
		int nu2 = data->ints[ni1 + ni2 + 1]->len;

		double retE[M];
		double daux;
		int ox;

		k = 0;
		for (i = 0; i < M1; i++) {
			daux = ret1[2 + i];
			for (j = 0; j < M2; j++) {
				retE[k++] = daux * ret2[2 + j];
			}
		}

		if ((nu1 > 0) & (nu2 > 0)) {
			for (i = 0; i < nu1; i++) {
				daux = ret1[2 + data->ints[ni1 + ni2]->ints[i]];
				for (j = 0; j < nu2; j++) {
					retE[k++] = daux * ret2[2 + data->ints[ni1 + ni2 + 1]->ints[j]];
				}
			}
		}

		assert(k == data->smats[nsm1 + nsm2]->n);
		for (k = 0; k < data->smats[nsm1 + nsm2]->n; k++) {
			ox = (int) data->smats[nsm1 + nsm2]->x[k];
			ret[2 + k] = retE[ox];
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

		ret1 = c->model1_func(INLA_CGENERIC_INITIAL, NULL, dataM1);
		ret2 = c->model2_func(INLA_CGENERIC_INITIAL, NULL, dataM2);

		int nth1 = (int) ret1[0], nth2 = (int) ret2[0];

		ret = Calloc(1 + nth1 + nth2, double);
		ret[0] = nth1 + nth2;

		for (i = 0; i < nth1; i++) {
			ret[1 + i] = ret1[1 + i];
		}

		for (i = 0; i < nth2; i++) {
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
		ret1 = c->model1_func(INLA_CGENERIC_LOG_PRIOR, theta1, dataM1);
		ret2 = c->model2_func(INLA_CGENERIC_LOG_PRIOR, theta2, dataM2);

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
