#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dporfs_(char *uplo, integer *n, integer *nrhs, doublerealRopt *a, integer *lda, doublerealRopt *af, integer *ldaf, doublerealRopt *b, integer *ldb, doublerealRopt *x, integer *ldx, doublerealRopt *ferr, doublerealRopt *berr, doublerealRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif