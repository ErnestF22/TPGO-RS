#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int chprfs_(char *uplo, integer *n, integer *nrhs, complexRopt *ap, complexRopt *afp, integer *ipiv, complexRopt *b, integer *ldb, complexRopt *x, integer *ldx, realRopt *ferr, realRopt *berr, complexRopt *work, realRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif