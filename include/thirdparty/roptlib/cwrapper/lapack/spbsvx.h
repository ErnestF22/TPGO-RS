#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int spbsvx_(char *fact, char *uplo, integer *n, integer *kd, integer *nrhs, realRopt *ab, integer *ldab, realRopt *afb, integer *ldafb, char *equed, realRopt *s, realRopt *b, integer *ldb, realRopt *x, integer *ldx, realRopt *rcond, realRopt *ferr, realRopt *berr, realRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif