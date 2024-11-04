#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cpbrfs_(char *uplo, integer *n, integer *kd, integer *nrhs, complexRopt *ab, integer *ldab, complexRopt *afb, integer *ldafb, complexRopt *b, integer *ldb, complexRopt *x, integer *ldx, realRopt *ferr, realRopt *berr, complexRopt *work, realRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif