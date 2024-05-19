
/* cgeneric_defs.h
 *
 * Copyright (C) 2022-2023 Elias Krainski
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

#ifndef __CGENERIC_DEFS_H__
#define __CGENERIC_DEFS_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS					       /* empty */
#define __END_DECLS					       /* empty */
#endif

__BEGIN_DECLS
#include <stddef.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "cgeneric.h"
#include "corgraphs_utils.h"

#if !defined(Calloc)
#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#endif
#define pow2(x) ((x)*(x))
#define pow3(x) (pow2(x)*(x))
#define pow4(x) (pow2(x)*pow2(x))
#ifdef __SUPPORT_SNAN__
#define iszero(x) (fpclassify(x) == FP_ZERO)
#else
#define iszero(x) (((__typeof(x))(x)) == 0)
#endif

inla_cgeneric_func_tp inla_cgeneric_besag;
inla_cgeneric_func_tp inla_cgeneric_corgraphs0;
inla_cgeneric_func_tp inla_cgeneric_corgraphs_sfixed;
inla_cgeneric_func_tp inla_cgeneric_corgraphs_sch;
inla_cgeneric_func_tp inla_cgeneric_corgraphs_cholQ;

__END_DECLS
#endif
