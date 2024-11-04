#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sgtsvx_(char *fact, char *trans, integer *n, integer *nrhs, realRopt *dl, realRopt *d__, realRopt *du, realRopt *dlf, realRopt *df, realRopt *duf, realRopt *du2, integer *ipiv, realRopt *b, integer *ldb, realRopt *x, integer *ldx, realRopt *rcond, realRopt *ferr, realRopt *berr, realRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif