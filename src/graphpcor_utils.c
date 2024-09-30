
/* graphpcor_utils.c
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
 *        Elias T Krainski
 *        CEMSE Division
 *        King Abdullah University of Science and Technology
 *        Thuwal 23955-6900, Saudi Arabia
 */

#include "graphpcor.h"

void exchangeableU(int n, double r, double *cc) {
  if(n==1) {

    cc[0] = 1.0;

  } else {

    // upper (C), lower (Fortran)
    int i, j, k=0;
    for(i=0; i<n; i++) {
      for(j=i; j<n; j++) {
        if(i==j) {
          cc[k++] = 1.0;
        } else {
          cc[k++] = r;
        }
      }
    }

/*
    k=0;
    for(i=0; i<n; i++) {
      for(j=i; j<n; j++) {
        printf("%f ", dlQ[k++]);
      }
      printf("\n");
    }
*/

  }

}

void dl2Qu(int n, double *d, double *l, double *qu) {
  assert(n>0);
  if(n==1) {
    qu[0] = (d[0])*(d[0]);
  }
  if(n==2) {
    qu[0] = d[0]*d[0];
    qu[1] = d[0]*l[0];
    qu[2] = d[1]*d[1] + l[0]*l[0];
  }
  if(n>2) {
    double qq[n*n];
    dl2fullQ(n, &d[0], &l[0], &qq[0]);
    // copy (as the upper part) to the vector 'qu'
    int i,j,k = 0, k2 = 0;
    for(i=0; i<n; i++) {
      k2 = (n+1)*i;
      for(j=i; j<n; j++) {
        qu[k++] = qq[k2++];
      }
    }
  }

}

void dl2fullQ(int n, double *d, double *l, double *qq) {

    assert(n>0);
    if(n==1) {
      qq[0] = (d[0])*(d[0]);
    }
    if(n==2) {
      qq[0] = d[0]*d[0];
      qq[1] = d[0]*l[0];
      qq[2] = qq[1];
      qq[3] = d[1]*d[1] + l[0]*l[0];
    }
    if(n>2) {

      double aa[n*n], bb[n*n];
      int i,j,k=0, k2=0;
      /*
      for(i=0; i<(n-1); i++) {
        k = (n+1)*i;
        aa[k++] = d[i];
        for(j=(i+1); j<n; j++) {
          aa[k++] = l[k2++];
        }
      }
      aa[n*n-1] = d[n-1];
       */
      k=0; k2=0;
      for(i=0; i<n; i++) {
        for(j=0; j<n; j++) {
          if(j>i) {
            aa[k] = l[k2];
            bb[k] = l[k2];
            k2++;
          } else {
            if(i==j) {
              aa[k] = d[i];
              bb[k] = d[i];
            } else {
              aa[k] = 0.0;
              bb[k] = 0.0;
            }
          }
          k++;
        }
      }

/*
      printf("\nPrinting aa:\n");
      k=0;
      for(i=0; i<n; i++) {
        for(j=0; j<n; j++) {
          printf("%f ", aa[k++]);
        }
        printf("\n");
      }
      printf("\nPrinting bb:\n");
      k=0;
      for(i=0; i<n; i++) {
        for(j=0; j<n; j++) {
          printf("%f ", bb[k++]);
        }
        printf("\n");
      }
*/

      char tra = 'N'; // upper in C is lower in Fortran
      char trb = 'T';
      double alpha = 1, beta = 0;
      dgemm_(&tra, &trb, &n, &n, &n,
             &alpha, &aa[0], &n, &bb[0], &n,
             &beta, &qq[0], &n, F_ONE);

/*
      printf("\nPrinting qq:\n");
      k=0;
      for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
          printf("%f ", qq[k++]);
        }
        printf("\n");
      }
*/

    }

}

void cov2cor(int verbose, int n, double *cc) {
  double s[n];
  int i, j, k=0;
  for(i=0; i<n; i++) {
    s[i] = sqrt(cc[k]);
    k += n+1;
    //    printf("i = %d, k = %d, %2.1f ", i, k, s[i]);
  }

  if(verbose>999) {
    printf("\ns[i]:\n");
    for(i=0; i<n; i++) {
      printf("%2.1f ", s[i]);
    }
  }

  k=0;
  for(i=0; i<n; i++) {
    for(j=0; j<n; j++) {
      cc[k++] /= (s[i] * s[j]);
    }
  }

  if(verbose>9999) {
    printf("cc[i,j]:\n");
    k=0;
    for(i=0; i<n; i++) {
      for(j=0; j<n; j++) {
        printf("%2.3f ", cc[k++]);
      }
      printf("\n");
    }
  }

}

double covariance_kld(int verbose, int n, double *C0, double *C1) {
  char uplo = 'U';
  int info=0, i, k, n2 = n*n;
  double hldet0=0, hldet1=0, trc = 0, kld = 0;
  double cc0[n2], cc1[n2];

  for(k=0; k<n2; k++) {
    cc0[k] = C0[k]; // copy
    cc1[k] = C1[k]; // copy
  }

  dpotrf_(&uplo, &n, &cc0[0], &n, &info, F_ONE);
  dpotrf_(&uplo, &n, &cc1[0], &n, &info, F_ONE);

  k=0;
  for(i=0; i<n; i++) {
    hldet0 += log(cc0[k]);
    hldet1 += log(cc1[k]);
    k += n+1;
  }

  dposv_(&uplo, &n, &n, &C0[0], &n, &C1[0], &n, &info, F_ONE);

  // trace of C1/C0
  k=0;
  for(i=0; i<n; i++) {
    trc += C1[k];
    k += n+1;
  }

  kld += 0.5 * (trc - n) - hldet1 + hldet0;

  if(verbose>99) {
    printf("ld0: %2.4f, ld1: %2.4f, tr(C1/C0): %2.4f, kld:= %2.4f\n",
           hldet0, hldet1, trc, kld);
  }
  return kld;
}

void covariance_parent_children(int verbose, int np, int N, int niiv, int *iiv, int *jjv, int *ipar, int *itop, double *sch, double *v2, double *CC) {

  int i, j, k;
  double vp[np];

  for(i=0; i<np; i++) {
    vp[i] = 0.0;
  }
  for(i=0; i<niiv; i++) {
    vp[iiv[i]] += v2[jjv[i]];
  }

  k=0;
  for(i=0; i<N; i++) {
    for(j=0; j<N; j++) {
      if(j==i) {
        CC[k] = 1.0 + vp[ipar[i]];
      } else {
        CC[k] = sch[i] * sch[j] * vp[itop[k]];
      }
      k++;
    }
  }

  if(verbose>999){
    printf("vp[i]:\n");
    for(i=0; i<np; i++) {
      printf("%2.1f ", vp[i]);
    }
    printf("\nCC[i,j]:\n");
    k=0;
    for(i=0; i<N; i++) {
      for(j=0; j<N; j++) {
        printf("%2.3f ", CC[k]);
        k++;
      }
      printf("\n");
    }
  }

}

void correlation_parent_children(int verbose, int np, int N, int niiv, int *iiv, int *jjv, int *ipar, int *itop, double *sch, double *v2, double *CC) {
  covariance_parent_children(verbose, np, N, niiv, &iiv[0], &jjv[0], &ipar[0], &itop[0], &sch[0], &v2[0], &CC[0]) ;
  cov2cor(verbose, N, &CC[0]) ;
}

void theta_parent_children_kldh(int verbose, int np, int N, int niiv, int *iiv, int *jjv, int *ipar, int *itop, double *sch, double *theta, double hs, double *kld, double *kldh) {

  int i, l, k, N2 = N*N;
  double v2a[np], v2b[np];
  double C0[N2], C1[N2], C0c[N2], C1h[N2];

  for(i=0; i<np; i++) {
    v2a[i] = exp(2.0 * theta[i]);
    v2b[i] = v2a[i];
  }

  correlation_parent_children(verbose, np, N, niiv, &iiv[0], &jjv[0], &ipar[0], &itop[0], &sch[0], &v2a[0], &C1[0]);

  for(l=np; l>0; l--) {

    if(verbose>99) {
      printf("l is now %d\n", l);
    }

    if(l<np) {
      v2a[l] = 0.0; // zero all above last
    }
    v2a[l-1] = exp(2.0 * (theta[l-1] + hs));
    if(verbose>99) {
      printf("C1h:\n");
    }
    correlation_parent_children(verbose, np, N, niiv, &iiv[0], &jjv[0], &ipar[0], &itop[0], &sch[0], &v2a[0], &C1h[0]);

    if(verbose>99) {
      printf("C0:\n");
    }
    v2b[l-1] = 0.0; // from last one is zero
    correlation_parent_children(verbose, np, N, niiv, &iiv[0], &jjv[0], &ipar[0], &itop[0], &sch[0], &v2b[0], &C0[0]);
    for(k=0; k<N2; k++) {
      C0c[k] = C0[k]; // copy
    }

    kld[l-1] = covariance_kld(verbose, N, &C0[0], &C1[0]);

    if(l>1) {      // prepare for next step
      for(k=0; k<N2; k++) {
        C1[k] = C0c[k]; // copy one-param reduced matrix
      }
    }

    kldh[l-1] = covariance_kld(verbose, N, &C0c[0], &C1h[0]);

    if(verbose>999) {
      printf("k,d: %2.5f %2.5f\n", kld[l-1], kldh[l-1]);
    }

  } // end l

}

