#include <R.h>
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include "lar.h"

/* for BLAS */
static int IONE = 1 ;
static double ZERO = 0 , ONE = 1 , ONEM = -1, HALF = 0.5, TWO = 2;
static char tran = 'N' ;

/* calling functions are responsable for the [put/get]Rng thing */
static void ggsimy(int t, double *chic, int tau, double *delta, double *choc,
		 int n, double *y, double *e) {
    int i ;
    R_CheckUserInterrupt() ;
    for ( i = 0 ; i < n ; i++) e[i] = norm_rand() ;
    if ( t < tau ) {
	F77_CALL(dgemv)(&tran, &n, &n, &ONE, chic , &n, e, &IONE, &ZERO, y, &IONE);
    } else {
	F77_CALL(dgemv)(&tran, &n, &n, &ONE, choc , &n, e, &IONE, &ZERO, y, &IONE);
	F77_CALL(daxpy)(&n, &ONE, delta , &IONE, y, &IONE );
    }
}

/* control statistics update functions */
typedef int (*chartfn)(int, int , double *,int *, double *, double *);

/* 
LAR 
- ipar=(p,only.mean)
- par=(h,l,S^{-1},F'S^{-1}F,F'S^{-1},a,b)
- stat=(w,s,v,beta,sxy,space,wlar) length(stat)=1+(p+2)+p*p+p+max(n,p)+3+7*p+2*p*p 
*/
static void gglarupd(int t, int n, double *y, int *ipar, double *par, double *stat) {
    int i, p = ipar[0] ;
    double h = par[0] , l = par[1] , l1 = 1 - l, vl = (2-l)/l ,
	*si=par+2, *sxx=si+n*n, *fsi=sxx+p*p, *tmean=fsi+n*p, *tsd=tmean+p+1,
	*s=stat+1, *beta=s+p+2, *sxy=beta+p*p, *space=sxy+p, *wlar=space+n;  
    F77_CALL(dgemv)(&tran, &p, &n, &l, fsi, &p, y, &IONE, &l1, sxy, &IONE);
    gglarstart(p, ESTLS, sxx, sxy, wlar);
    for ( i = 0 ; i < p ; i++, beta += p) {
	gglarnext( wlar, beta ) ;
	F77_CALL(dgemv)(&tran, &p, &p, &ONE, sxx, &p, beta, &IONE, &ZERO, space, &IONE);
	s[i] = vl * F77_CALL(ddot)(&p , beta, &IONE, space, &IONE) ;
    }
    F77_CALL(dgemv)(&tran, &n , &n , &ONE, si , &n, y, &IONE, &ZERO, space, &IONE);
    s[p] = s[p+1] = fmax2(ONE, l1*s[p+1] + l*F77_CALL(ddot)(&n, y, &IONE, space , &IONE)/n);
}

static int ggspclar(int t, int n, double *y, int *ipar, double *par, double *stat) {
    int i , p = ipar[0] , nt = (ipar[1]) ? p : p+1 ;
    double h=par[0], *s=stat+1, 
	*si=par+2, *sxx=si+n*n, *fsi=sxx+p*p, *tmean=fsi+n*p, *tsd=tmean+p+1;
    gglarupd(t, n , y, ipar, par, stat);
    for ( i = 0 ; i <= p ; i++) s[i] = (s[i] - tmean[i]) / tsd[i] ;
    for (i=1, stat[0]=s[0]; i < nt ; i++) stat[0] = fmax2( stat[0], s[i] ) ;
    return( (stat[0] <= h) ? 0 : 1 ) ;
}

void gglarcompmom(int *horiz, double *chic, int *ipar, double *par, double *stat0) {
    int i, j, h = horiz[0], tau = h+1, n = ipar[1], allstat = ipar[2], p = ipar[7] ;
    double *y = (double *) R_alloc( 2*n+allstat , sizeof(double) ),
	*e = y + n, *stat = e + n, *s = stat + 1, 
	*si=par+2, *sxx=si+n*n, *fsi=sxx+p*p, *tmean=fsi+n*p, *tsd=tmean+p+1;
    ipar = ipar + 7 ;
    memcpy( stat , stat0 , allstat * sizeof(double) ) ;
    for ( j = 0 ; j <= p ; j++) tmean[j] = tsd[j] = ZERO ;
    GetRNGstate();
    for ( i = 1 ; i <= h ; i++) {
	ggsimy( i , chic, tau, NULL, NULL , n , y , e ) ;
	gglarupd(i, n , y, ipar, par, stat);
	for ( j = 0 ; j <= p ; j++ ) {
	    tmean[j] += ( s[j] - tmean[j] ) / i ;
	    tsd[j] += ( s[j]*s[j] - tsd[j] ) / i ;
	}
    }
    PutRNGstate();
    for ( j = 0 ; j <= p ; j++ ) tsd[j] = sqrt(tsd[j]-tmean[j]*tmean[j]) ;
    tmean[p-1] = p ;
    tsd[p-1] = sqrt(2*p) ;
}

static chartfn updchart[] = {ggspclar};

static void ggupdatechart( int t , int *ipar, double *par,
			   double *y , double *stat, int *signal) {
    int i, nc = ipar[0], n = ipar[1] ;
    ipar = ipar + 3 ;
    for ( i = 0 ; i < nc ; i++ ) {
	signal[i] = updchart[ipar[0]](t,n,y,ipar+4,par,stat);
	par += ipar[2] ;
	stat += ipar[3] ;
	ipar += ipar[1] ;
    }
}

void ggcomputecharts(int *ipar, double *par, int *ny , double *y, double *stat0, double *ans) {
    int t , nc = ipar[0] , n = ipar[1] , allstat = ipar[2] ,
	*signal = (int *) R_alloc( nc , sizeof(int) ) ;
    double *stat = (double *) R_alloc( allstat, sizeof(double)) ;
    memcpy( stat , stat0 , allstat * sizeof(double) ) ;
    for ( t = 1 ; t <= (*ny) ; t++ , y += n , ans += allstat ) {
	ggupdatechart( t , ipar, par , y , stat, signal ) ;
	memcpy( ans , stat , allstat * sizeof(double) ) ;
    }
}


/* length(work) = length(stat)+2*n ; length(iwork) = nc */
static void ggsimrl1(int *ipar, double *par, double *stat0,
		     double *chic, int tau , double *delta, double *choc,
		     int maxrl, int *rl , int *iwork, double *work) {
    int j, nc = ipar[0], n = ipar[1] , allstat = ipar[2] ,
	done = 0 , t = 1;
    double *e = work , *y = e + n , *stat = y + n ;
    for ( j = 0 ; j < nc ; j++) rl[j] = maxrl ;
    memcpy( stat , stat0 , allstat * sizeof(double) ) ;
    while ( ( done < nc ) && ( t < maxrl ) ) {
	ggsimy(t, chic, tau, delta, choc, n, y , e ) ;
	ggupdatechart( t , ipar, par , y , stat, iwork);
	for ( j = 0 ; j < nc ; j++ ) {
	    if ( rl[j] == maxrl ) {
		if ( iwork[j] ) {
		    done++ ;
		    rl[j] = t ;
		}
	    }
	}
	t++ ;
    }
}

void ggsimrlss(int *ipar, double *par, double *stat0,
	       double *chic, int *tau , double *delta, double *choc,
	       int *irl, int *rl) {
    int i , nc = ipar[0] , n = ipar[1] , allstat = ipar[2] ,
	nrl = irl[0] , maxrl =irl[1] , *iwork = (int *) R_alloc( nc, sizeof(int) );
    double *work = (double *) R_alloc( 2*n + allstat, sizeof(double)) ;
    GetRNGstate();
    for ( i = 0 ; i < nrl ; i++ , rl += nc )
	ggsimrl1(ipar, par, stat0, chic, *tau, delta, choc, maxrl, rl, iwork, work) ;
    PutRNGstate();
}



/* control=(arl0,csi,burnin,iter,A,B) */
void ggmatcharl0one(int *ipar, double *par , double *stat0, double *chic, double *control) {
    int i , rl, signal, iwork, n = ipar[1] , allstat = ipar[2] , maxrl,
	burnin = control[2] + 0.5 , iter = control[3] + 0.5 ;
    double h = 0, gain,  arl0 = control[0] , csi = control[1],  A = control[4] , B = control[5] ,
	*work = (double *) R_alloc( allstat + 2*n , sizeof(double)) ;
    GetRNGstate();
    for ( i = 1 ; i <= iter+burnin ; i++) {
	gain = A / pow( i , B ) ;
	maxrl = arl0 * (1 + csi/gain) ;
	ggsimrl1(ipar, par, stat0, chic, maxrl+1, NULL, NULL, maxrl, &rl, &iwork, work);
	par[0] = fmax2( ZERO, par[0] - gain *( rl - arl0 ) / arl0 ) ;
	if ( i > burnin ) h += ( par[0] - h ) / ( i - burnin ) ;
    }
    PutRNGstate();
    par[0] = h ;
}

