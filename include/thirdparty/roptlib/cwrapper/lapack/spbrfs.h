#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int spbrfs_(char *uplo, integer *n, integer *kd, integer *nrhs, realRopt *ab, integer *ldab, realRopt *afb, integer *ldafb, realRopt *b, integer *ldb, realRopt *x, integer *ldx, realRopt *ferr, realRopt *berr, realRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif