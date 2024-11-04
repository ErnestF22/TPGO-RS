#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sgbrfs_(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, realRopt *ab, integer *ldab, realRopt *afb, integer *ldafb, integer *ipiv, realRopt *b, integer *ldb, realRopt *x, integer *ldx, realRopt *ferr, realRopt *berr, realRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif