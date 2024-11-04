#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int ztprfs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublecomplexRopt *ap, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *x, integer *ldx, doublerealRopt *ferr, doublerealRopt *berr, doublecomplexRopt *work, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif