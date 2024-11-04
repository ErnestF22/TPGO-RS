#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sgerfs_(char *trans, integer *n, integer *nrhs, realRopt *a, integer *lda, realRopt *af, integer *ldaf, integer *ipiv, realRopt *b, integer *ldb, realRopt *x, integer *ldx, realRopt *ferr, realRopt *berr, realRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif