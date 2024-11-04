#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zppsvx_(char *fact, char *uplo, integer *n, integer *nrhs, doublecomplexRopt *ap, doublecomplexRopt *afp, char *equed, doublerealRopt *s, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *x, integer *ldx, doublerealRopt *rcond, doublerealRopt *ferr, doublerealRopt *berr, doublecomplexRopt *work, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif