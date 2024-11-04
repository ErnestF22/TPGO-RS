#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sspsvx_(char *fact, char *uplo, integer *n, integer *nrhs, realRopt *ap, realRopt *afp, integer *ipiv, realRopt *b, integer *ldb, realRopt *x, integer *ldx, realRopt *rcond, realRopt *ferr, realRopt *berr, realRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif