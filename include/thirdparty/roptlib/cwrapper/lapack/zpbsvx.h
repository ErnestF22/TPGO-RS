#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zpbsvx_(char *fact, char *uplo, integer *n, integer *kd, integer *nrhs, doublecomplexRopt *ab, integer *ldab, doublecomplexRopt *afb, integer *ldafb, char *equed, doublerealRopt *s, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *x, integer *ldx, doublerealRopt *rcond, doublerealRopt *ferr, doublerealRopt *berr, doublecomplexRopt *work, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif