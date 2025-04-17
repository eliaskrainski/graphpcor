
/* cgeneric_generic0.c
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

#include "graphpcor.h"

double *inla_cgeneric_generic0(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{

	// Q = \tau*R

	double *ret = NULL, prec;

	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	int N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

	assert(!strcasecmp(data->ints[1]->name, "debug"));     // this will always be the case
//	int debug = data->ints[1]->ints[0];

	int fixed = 1;
	if ((data->doubles[0]->doubles[1] > 0) & (data->doubles[0]->doubles[1] < 1)) {
		fixed = 0;
	}

	if (fixed) {
		prec = 1 / pow2(data->doubles[0]->doubles[0]);
	} else {
		if (theta) {
			prec = exp(theta[0]);
		} else {
			prec = NAN;
		}
	}

	// the number of neighbors at the diagonal
	// and the (upper) graph (-1)
	// NOTE: scale this above def one!!!
	assert(data->n_smats > 0);
	assert(!strcasecmp(data->smats[0]->name, "Rgraph"));
	inla_cgeneric_smat_tp *Rgraph = data->smats[0];
	assert(N == Rgraph->nrow);
	assert(N == Rgraph->ncol);
	int M = Rgraph->n;
	/*
	if (debug > 99) {
		printf("M = %d\n", M);
	}
	*/

	switch (cmd) {

	case INLA_CGENERIC_VOID:
	{
		assert(!(cmd == INLA_CGENERIC_VOID));
		break;
	}

	case INLA_CGENERIC_GRAPH:
	{
		int offset = 2;
		ret = Calloc(2 + 2 * M, double);
		ret[0] = N;				       /* dimension */
		ret[1] = M;				       /* number of (i <= j) */
		for (int k = 0; k < M; k++) {
			ret[offset + k] = Rgraph->i[k];
			ret[offset + M + k] = Rgraph->j[k];
		}
		break;
	}

	case INLA_CGENERIC_Q:
	{
		int offset = 2;
		ret = Calloc(2 + M, double);
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;
		for (int k = 0; k < M; k++) {
			ret[offset + k] = prec * Rgraph->x[k];
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
		if (fixed) {
			ret = Calloc(1, double);
			ret[0] = 0;
		} else {
			ret = Calloc(2, double);
			ret[0] = 1;
			ret[1] = 0.0;
		}
		break;
	}

	case INLA_CGENERIC_LOG_NORM_CONST:
	{
		break;
	}

	case INLA_CGENERIC_LOG_PRIOR:
	{
		// return c(LOG_PRIOR), PC-PREC
		ret = Calloc(1, double);
		if (fixed) {
			ret[0] = 0.0;
		} else {
			double val = 0.5 * theta[0];
			assert(!strcasecmp(data->doubles[0]->name, "param"));
			double u = data->doubles[0]->doubles[0];
			double a = data->doubles[0]->doubles[1];
			double l = -log(a) / u;
/*
 			if (debug > 99) {
				printf("th %f, u %f, a %f and l %f\n", theta[0], u, a, l);
			}
*/
			if (u <= 0) {
				ret[0] = 0.0;
			} else {
				if (a <= 0.0) {
					ret[0] = 0.0;
				} else {
					if (a >= 1.0) {
						ret[0] = 0.0;
					} else {
						ret[0] = log(0.5 * l) - l * exp(-val) - val;
					}
				}
			}
		}
		break;
	}

	case INLA_CGENERIC_QUIT:
	default:
		break;
	}

	return (ret);
}
