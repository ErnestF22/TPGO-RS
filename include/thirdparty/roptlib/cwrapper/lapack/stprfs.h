#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int stprfs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, realRopt *ap, realRopt *b, integer *ldb, realRopt *x, integer *ldx, realRopt *ferr, realRopt *berr, realRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif