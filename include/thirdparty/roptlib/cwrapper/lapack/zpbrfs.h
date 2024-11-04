#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zpbrfs_(char *uplo, integer *n, integer *kd, integer *nrhs, doublecomplexRopt *ab, integer *ldab, doublecomplexRopt *afb, integer *ldafb, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *x, integer *ldx, doublerealRopt *ferr, doublerealRopt *berr, doublecomplexRopt *work, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif