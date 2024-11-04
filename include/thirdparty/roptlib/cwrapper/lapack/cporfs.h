#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cporfs_(char *uplo, integer *n, integer *nrhs, complexRopt *a, integer *lda, complexRopt *af, integer *ldaf, complexRopt *b, integer *ldb, complexRopt *x, integer *ldx, realRopt *ferr, realRopt *berr, complexRopt *work, realRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif