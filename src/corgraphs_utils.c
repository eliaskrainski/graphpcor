
/* corgraphs_utils.c
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

#include "corgraphs_utils.h"

double variance_parent_children_kld(int verbose, int np, int N, int niiv, int *iiv, int *jjv, int *ipar, int *itop, double *sch, double *v2) {

  if(verbose>999)
    printf("working inside variance_parent_children_kld now\n");

  char uplo = 'U';
  int info=0, N2 = N*N, l, i, j, k;
  double hldet0, hldet1, trc, kld = 0.0;
  double v2a[np], v2b[np], vp0[np], vp1[np], s0[N], s1[N];
  double C0[N2], C1[N2], cc0[N2], cc1[N2];

  for(i=0; i<np; i++) {
    v2a[i] = v2[i]; // copy
    v2b[i] = v2[i]; // copy
  }

  for(l=np; l>0; l--) {

    if(verbose>999) {
      printf("l is now %d\n", l);
    }

    v2a[l-1] = 0.0;
    if(l<np) {
      v2b[l]= 0.0;
    }

    if(verbose>999){
      printf("v2a[i]:\n");
      for(i=0; i<np; i++) {
        printf("%2.1f ", v2a[i]);
      }
      printf("\n");
      printf("v2b[i]:\n");
      for(i=0; i<np; i++) {
        printf("%2.1f ", v2b[i]);
      }
      printf("\n");
    }
    for(i=0; i<np; i++) {
      vp0[i] = 0.0;
      vp1[i] = 0.0;
    }
    for(i=0; i<niiv; i++) {
      vp0[iiv[i]] += v2a[jjv[i]];
      vp1[iiv[i]] += v2b[jjv[i]];
    }

    k=0;
    for(i=0; i<N; i++) {
      for(j=0; j<N; j++) {
        if(j==i) {
          C0[k] = 1.0 + vp0[ipar[i]];
          C1[k] = 1.0 + vp1[ipar[i]];
          s0[i] = sqrt(C0[k]);
          s1[i] = sqrt(C1[k]);
        } else {
          C0[k] = sch[i] * sch[j] * vp0[itop[k]];
          C1[k] = sch[i] * sch[j] * vp1[itop[k]];
        }
        k++;
      }
    }

    if(verbose>999){
      printf("vp0[i]:\n");
      for(i=0; i<np; i++) {
        printf("%2.1f ", vp0[i]);
      }
      printf("\nvp1[i]:\n");
      for(i=0; i<np; i++) {
        printf("%2.1f ", vp1[i]);
      }
      printf("\ns0[i]:\n");
      for(i=0; i<N; i++) {
        printf("%2.1f ", s0[i]);
      }
      printf("\ns1[i]:\n");
      for(i=0; i<N; i++) {
        printf("%2.1f ", s1[i]);
      }
      printf("\nC0[i,j]:\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=0; j<N; j++) {
          printf("%2.3f ", C0[k]);
          k++;
        }
        printf("\n");
      }
      printf("C1[i,j]:\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=0; j<N; j++) {
          printf("%2.3f ", C1[k]);
          k++;
        }
        printf("\n");
      }
    }
    k=0;
    for(i=0; i<N; i++) {
      for(j=0; j<N; j++) {
        C0[k] /= (s0[i] * s0[j]);
        C1[k] /= (s1[i] * s1[j]);
        cc0[k] = C0[k]; // copy
        cc1[k] = C1[k]; // copy
        k++;
      }
    }

    if(verbose>999){
      printf("C0[i,j]:\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=0; j<N; j++) {
          printf("%2.3f ", C0[k]);
          k++;
        }
        printf("\n");
      }
      printf("C1[i,j]:\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=0; j<N; j++) {
          printf("%2.3f ", C1[k]);
          k++;
        }
        printf("\n");
      }
    }

    dpotrf_(&uplo, &N, &cc0[0], &N, &info, F_ONE);
    if(verbose>99) {
      printf("INFO for L0 with l = %d is %d\n", l, info);
    }
    dpotrf_(&uplo, &N, &cc1[0], &N, &info, F_ONE);
    if(verbose>99) {
      printf("INFO for L1 with l = %d is %d\n", l, info);
    }

    if(verbose>999){
      printf("L0[i,j]:\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=0; j<N; j++) {
          printf("%2.3f ", C0[k]);
          k++;
        }
        if(verbose>1){
          printf("\n");
        }
      }
      printf("L1[i,j]:\n");
      k=0;
      for(i=0; i<N; i++) {
        for(j=0; j<N; j++) {
          printf("%2.3f ", C1[k]);
          k++;
        }
        if(verbose>1){
          printf("\n");
        }
      }
    }

    hldet0 = 0.0;
    hldet1 = 0.0;
    k=0;
    for(i=0; i<N; i++) {
      hldet0 += log(cc0[k]);
      hldet1 += log(cc1[k]);
      if(verbose>999){
        printf("ld0 = %2.5f, ld1 = %2.5f\n", cc0[k], cc1[k]);
      }
      k += N+1;
    }

    dposv_(&uplo, &N, &N, &C1[0], &N, &C0[0], &N, &info, F_ONE);
    if(verbose>99) {
      printf("INFO for dposv with l = %d is %d\n", l, info);
    }

    // add trace of C1/C0
    k=0;
    for(i=0; i<N; i++) {
      trc += C0[k];
      k += N+1;
    }

    kld += 0.5 * (trc - np) - hldet1 + hldet0;

    if(verbose>99) {
      printf("ld0: %2.4f, ld1: %2.4f, tr(C1/C0): %2.4f, kld:= %2.4f\n",
             hldet0, hldet1, trc, kld);
    }

  }

  return kld;

}
