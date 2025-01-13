#define ESTLAR 0
#define ESTLS 1

/* length(wlar)=3+7*p+2*p*p */
void gglarstart(int p, int est, double *sxx, double *sxy, double *wlar);
void gglarnext(double *wlar, double *beta) ;
void gglar(int *ipar, double *sxx, double *sxy, double *beta) ;
