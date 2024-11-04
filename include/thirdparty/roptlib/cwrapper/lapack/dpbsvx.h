#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dpbsvx_(char *fact, char *uplo, integer *n, integer *kd, integer *nrhs, doublerealRopt *ab, integer *ldab, doublerealRopt *afb, integer *ldafb, char *equed, doublerealRopt *s, doublerealRopt *b, integer *ldb, doublerealRopt *x, integer *ldx, doublerealRopt *rcond, doublerealRopt *ferr, doublerealRopt *berr, doublerealRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif