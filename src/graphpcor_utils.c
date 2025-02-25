
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

void theta2Qcorrel(int n, int std, double *theta, double *Qcorr) {

  int i, j, k, k0=0;

  // build the lower triangle L_0, Q_0 = L_0(L_0')
  for(i=0; i<n; i++) {
    for(j=0; j<=i; j++) {
      nij2Lk(n,i,j);
      if(i==j) {
        Qcorr[k] = 1.0;
      } else {
        Qcorr[k] = theta[k0++];
      }
    }
  }

  if(std) {
    // chol2inv: to compute V = Q_0^{-1}
    int info;
    char uplo = 'L';
    dpptri_(&uplo, &n, &Qcorr[0], &info, F_ONE);

    // diag(V)^{1/2}
    double si[n];
    k = 0;
    for(i=0; i<n; i++) {
      si[i] = sqrt(Qcorr[k]);
      k += n-i;
    }

    // Q = SVS :  the precision of a correlation matrix
    for(i=0; i<n; i++) {
      k = i;
      for(j=0; j<=i; j++) {
        Qcorr[k] *= 1/(si[i]*si[j]);
        k += (n-j-1);
      }
    }
  }

}

void theta2precision(int n, double *l, double *Q) {
  // return Q = LL' where
  // diag(L) = exp(l_1 ... l_n)
  //   l_n+1 ... l_m, m = n*(n+1)/2
  // goes in the lower triangle of L
  // Example for n = 3
  //      | exp(l_1)      0         0      |
  //  L = |   l_4     exp(l_2)      0      |
  //      |   l_5        l_6     exp(l_3)  |
  int i;
  double d[n];
  for(i=0; i<n; i++) {
    d[i] = exp(l[i]);
  }
  dl2fullQ(n, &d[0], &l[n], &Q[0]);
}

void theta2gamma2Lcorr(int n, double *hldet, double *theta, double *L) {
  // hypershere decomposition, Rapisarda, Brigo and Mercurio (2007).
  // from m=n(n-1)/2 length theta to n x n lower triangle L
  //  for a correlation matrix R = LL'
  // ONLY the lower triangle (with diagonal) is returned.
  // 1. \theta[k] \in (-\infty, \infty) to x[k] \in (0, \pi])
  //    x[k] = pi/(1+exp(-theta[k]))
  // 2. compute cos(x[i,j]) and sin(x[i,j])
  //  co[i,j] = cos(x[i,j])
  //  si[i,j] = sin(x[i,j])
  // 3. for j=1 (j=0 in C)
  //   L[1,1] = 1;
  //  and for i>1
  //   L[i,j] = co[i,j]
  // 4. compute (j>1 | j>0 in C)
  //  si[i,j] = \prod_{k=1}^{j-1} sin(x[i,k])
  // 5. for j>1
  //            { co[i,j]s[i,j-1],  2 <= j <= i-1
  //   L[i,j] = { si[i,j-1]      ,      j=i
  //            { 0              , j+1 <= j <= n
  assert(n>1);
  L[0] = 1.0;
  int i, j, k, k0, k1, n1=n-1;
  int m = (int)(((double)n)*((double)(n1))*0.5);
  double daux, co[m], si[m];

  // x[k] = pi/(1+exp(-theta[k]));
  // compute cos(x[k]) and sin(x[k]);
  for(i=0; i<m; i++) {
    daux = M_PI / ( 1 + exp(-theta[i]) );
    co[i] = cos(daux);
    si[i] = sin(daux);
  }

  if(n==2) {
    L[1] = co[0];
    L[2] = si[0];
    hldet[0] = L[2];
  }
  if(n==3) {
    L[1] = co[0];
    L[2] = co[1];
    L[3] = si[0];
    L[4] = co[2]*si[1];
    L[5] = si[1]*si[2];
    hldet[0] = L[3] + L[5];

  }
  if(n>3) {
    //    printL(L, n, n, "L[i,j]\n");
    //    printL(co, n1, n1, "cos[i,j]\n");
    //    printL(si, n1, n1, "sin[i,j]\n");
    // s[i,j] = \prod_{k=0}^{j-1} sin(x[i,k])

    for(i=1; i<n1; i++) {
      for(j=0; j<i; j++) {
        nij2Lk(n1,i,j);
        k0 = k;
        nij2Lk(n1,i,j+1);
        k1 = k;
        si[k1] *= si[k0]; // \prod{k=0}^{j-1} sin(x[i,k])
      }
    }

    // build L[,1]
    k0 = 0;
    k1 = n;
    for(i=1; i<n; i++) {
      L[i] = co[i-1];       // L[, 1]
      L[k1] = si[k0];       // diag(L)
      k0 += n-i;
      k1 += n-i;
    }

//    printL(si, n-1, n-1, "cum prod of sin[i,j]\n");
// L[lower,2:n] is now just si[2:n,2:j] * co[2:n,1:j-1]
    for(i=1; i<n1; i++) {
      for(j=0; j<i; j++) {
        nij2Lk(n1,i,j); // location of cum prod sin(x[i,j])
        k0 = k;
        nij2Lk(n1,i,j+1);  // loc of cos(x[i,j])
        k1 = k;
        nij2Lk(n,i+1,j+1); // loc of L[i,j]
        L[k] = co[k1] * si[k0];
      }
    }
//    printL(L, n, n, "L[i,j]\n");
  }
}

void theta2gamma2Ucorrel(int n, double *hldet, double *theta, double *cc) {
  double l[n*(n-1)/2];
  theta2gamma2Lcorr(n, &hldet[0], &theta[0], &l[0]);
  int i, j, k=0;
  for(i=0; i<n; i++) {
    cc[k] = l[k];
    k++;
  }
  for(i=1; i<n; i++) {
    for(j=0; j<=i; j++) {
      cc[k++] = 1;
    }
  }
}

void l2L(int n, double *l, double *L) {
  // l: lower triangle vector of length n(n-1)/2 l
  // L: nxn matrix L upper(C) lower (Fortran) triangle contain l
  int i,j,k = 0, k2 = 0;
  for(i=0; i<n; i++) {
    k2 = (n+1)*i;
    for(j=i; j<n; j++) {
      L[k2++] = l[k++];
    }
  }
}

void exchangeableU(int n, double r, double *cc) {
  // build upper (C), lower (Fortran)
  // exchangeable correlation matrix
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
  // build upper (C), lower (Fortran)
  // precision matrix from its (Cholesky?) factorized
  // L[full] = expand( l[lower triangle] )
  // Q = L'L
  // return upper(U)
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
  // build Q = L'L
  // return Q
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

void fillL(int *d, int *m, int *ii, int *jj, double *x) {
  // fill-in lower triangle from graph patterned elements
  int i, j, l, k, p, n = *d, nij = *m;
  for(l=0; l<nij; l++) {
    i = ii[l];
    j = jj[l];
    //  printf("%d: ", l);
    if(j>0) {
      p = j*n + i;
      for(k=0; k<j; k++) {
//        printf("%d: i %d, j %d, k %d ", p, i, j, k);
  //      printf("%2.2f %2.2f %2.2f \n", x[k*n+i], x[k*n+j], x[j*n+j]);
        x[p] -= x[k*n+i] * x[k*n+j] / x[j*n+j];
      }
    }
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
  // return ( tr(C1/C0) - n - |C1| - |C0| ) / 2
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

