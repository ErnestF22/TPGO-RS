#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zgbsvx_(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublecomplexRopt *ab, integer *ldab, doublecomplexRopt *afb, integer *ldafb, integer *ipiv, char *equed, doublerealRopt *r__, doublerealRopt *c__, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *x, integer *ldx, doublerealRopt *rcond, doublerealRopt *ferr, doublerealRopt *berr, doublecomplexRopt *work, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif