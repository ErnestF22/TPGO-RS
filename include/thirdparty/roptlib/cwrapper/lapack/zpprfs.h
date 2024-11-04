#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zpprfs_(char *uplo, integer *n, integer *nrhs, doublecomplexRopt *ap, doublecomplexRopt *afp, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *x, integer *ldx, doublerealRopt *ferr, doublerealRopt *berr, doublecomplexRopt *work, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif