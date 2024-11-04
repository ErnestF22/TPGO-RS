#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zgbrfs_(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublecomplexRopt *ab, integer *ldab, doublecomplexRopt *afb, integer *ldafb, integer *ipiv, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *x, integer *ldx, doublerealRopt *ferr, doublerealRopt *berr, doublecomplexRopt *work, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif