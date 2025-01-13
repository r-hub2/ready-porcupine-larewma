#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include "lar.h"

/* for BLAS */
static int IONE = 1 ;
static double ZERO = 0 , ONE = 1 , ONEM = -1;
static char tran = 'N' ;

#define L(i,j) m[act[j]*p+act[i]]
#define B(i) b[act[i]]
#define X(i) x[act[i]]

/* add which to the active set and update the cholesky decomposition 
   variable is not added if dependend by previous variables */
static void ggcholadd(int p, int *info, double *m, int which) {
    int i = info[0], *act = info+1, j, r;
    double tol = sqrt(DOUBLE_EPS); 
    act[i] = which;
    for ( j = 0 ; j < i ; j++ ) {
	for ( r = 0 ; r < j ; r++) L(i,j) -= L(i,r)*L(j,r) ;
	L(i,j) /= L(j,j) ;
    }
    for ( r = 0 ; r < i ; r++) L(i,i) -= L(i,r)*L(i,r) ;
    if ( L(i,i) > tol ) {
	L(i,i) = sqrt(L(i,i)) ;	
	info[0] = i + 1 ;
    } else {
	act[i] = 0 ;
    }
}

/* solve LL'b=x for the active set */
static void ggcholsolve(int p, int *info, double *m, double *b, double *x) {
    int na = info[0], *act = info+1, i , j ;
    for ( i = 0 ; i < p ; i++ ) x[i] = 0 ;
    for ( i = 0 ; i < na ; i++ ) {
	for ( j = 0 , X(i) = B(i) ; j < i ; j++) X(i) -= L(i,j) * X(j) ;
	X(i) /= L(i,i) ;
    }
    for ( i = na-1 ; i >= 0 ; i-- ) {
	for ( j = i+1 ; j < na ; j++) X(i) -= L(j,i) * X(j) ;
	X(i) /= L(i,i) ;
    }
}

/* static void ggcholmult(int p, int *info, double *m, double *b, double *x) { */
/*     int na = info[0], *act = info+1, i , j ; */
/*     for ( i = 0 ; i < p ; i++ ) x[i] = 0 ; */
/*     for ( i = 0 ; i < na ; i++ )  */
/* 	for ( j = i ; j < na ; j++) X(i) += L(j,i) * B(j) ; */
/* } */

#undef L
#undef X
#undef B


/* 
   wlar=(p,est,active,info,sxx,b,cholxx,cov,d,w,a) 
   length=1+1+p+(p+1)+p*p+p+p*p+4*p=3+7*p+2*p*p 
*/
void gglarstart(int p, int est, double *sxx, double *sxy, double *wlar) {
    int i, pp = p*p, *ilar = (int *) wlar , *active = ilar + 2;
    double *psxx=wlar+3+2*p, *b=psxx+p*p, *cholxx=b+p, *cov=cholxx+p*p, *d=cov+p;
    ilar[0] = p ;
    ilar[1] = est ;
    for ( i = 0 ; i < 2*p+1 ; i++) active[i] = 0 ;
    F77_CALL(dcopy)(&pp, sxx , &IONE, psxx, &IONE) ;
    F77_CALL(dcopy)(&pp, sxx , &IONE, cholxx, &IONE) ;
    F77_CALL(dcopy)(&p, sxy , &IONE, cov, &IONE) ;
    for ( i = 0 ; i < p ; i++ ) {
	b[i] = ZERO ;
	d[i] = sqrt(sxx[i*p+i]) ;
	cov[i] /= d[i] ;
    }
}

void gglarnext(double *wlar, double *beta) {
    int i, *ilar = (int *) wlar, p=ilar[0], est=ilar[1], *active=ilar+2, *info=active+p;
    double cmax, g, v, *sxx = wlar+3+2*p, *b = sxx + p*p, *cholxx = b + p ,
	*cov = cholxx + p*p, *d = cov + p, *w = d + p, *a = w + p , tol =sqrt(DOUBLE_EPS) ;
    for ( i = 0 , cmax = 0 ; i < p ; i++) cmax = fmax2( cmax , fabs(cov[i]) ) ;
    if ( cmax < tol ) return ;
    for ( i = 0 ; i < p ; i++ ) {
	if ( fabs(cov[i]) > cmax - tol ) {
	    a[i] = d[i]*sign(cov[i]) ;
	    if ( !active[i] ) {
		active[i] = IONE ;
		ggcholadd(p, info, cholxx, i ) ;
	    }
	}
    }
    ggcholsolve(p, info, cholxx, a, w ) ;
    F77_CALL(dgemv)(&tran, &p, &p, &ONE, sxx, &p, w, &IONE, &ZERO, a, &IONE);
    g = cmax ;
    for ( i = 0 ; i < p ; i++) {
	a[i] /= d[i] ;
	if ( !active[i] ) {
	    v = (cmax - cov[i]) / (1-a[i]) ;
	    if ( v > ZERO ) g = fmin2( g, v ) ;
	    v = (cmax + cov[i]) / (1+a[i]) ;
	    if ( v > ZERO ) g = fmin2( g, v ) ;
	}
    }
    v = -g ;
    F77_CALL(daxpy)(&p, &v, a, &IONE, cov, &IONE);
    F77_CALL(daxpy)(&p, &g , w , &IONE, b, &IONE) ;
    F77_CALL(dcopy)(&p, b, &IONE, beta, &IONE);
    if ( est == ESTLS ) {
	v = cmax - g ;
	F77_CALL(daxpy)(&p, &v , w , &IONE, beta, &IONE) ;	
    }
}

void gglar(int *ipar, double *sxx, double *sxy, double *beta) {
    int i, p = ipar[0], est = ipar[1]; 
    double *wlar = (double *) R_alloc(3+7*p+2*p*p, sizeof(double));
    gglarstart(p, est, sxx, sxy, wlar) ;
    for ( i = 0 ; i < p ; i++, beta += p ) gglarnext(wlar, beta) ;
}
