#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

doublerealRopt dlamch_(char *cmach);
int dlamc1_(integer *beta, integer *t, logical *rnd, logical *ieee1);
int dlamc2_(integer *beta, integer *t, logical *rnd, doublerealRopt *eps, integer *emin, doublerealRopt *rmin, integer *emax, doublerealRopt *rmax);
doublerealRopt dlamc3_(doublerealRopt *a, doublerealRopt *b);
int dlamc4_(integer *emin, doublerealRopt *start, integer *base);
int dlamc5_(integer *beta, integer *p, integer *emin, logical *ieee, integer *emax, doublerealRopt *rmax);

#ifdef __cplusplus
}
#endif