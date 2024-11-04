#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

E_f slamch_(char *cmach);
int slamc1_(integer *beta, integer *t, logical *rnd, logical *ieee1);
int slamc2_(integer *beta, integer *t, logical *rnd, realRopt *eps, integer *emin, realRopt *rmin, integer *emax, realRopt *rmax);
E_f slamc3_(realRopt *a, realRopt *b);
int slamc4_(integer *emin, realRopt *start, integer *base);
int slamc5_(integer *beta, integer *p, integer *emin, logical *ieee, integer *emax, realRopt *rmax);

#ifdef __cplusplus
}
#endif