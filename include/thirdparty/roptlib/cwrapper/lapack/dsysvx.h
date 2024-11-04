#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dsysvx_(char *fact, char *uplo, integer *n, integer *nrhs, doublerealRopt *a, integer *lda, doublerealRopt *af, integer *ldaf, integer *ipiv, doublerealRopt *b, integer *ldb, doublerealRopt *x, integer *ldx, doublerealRopt *rcond, doublerealRopt *ferr, doublerealRopt *berr, doublerealRopt *work, integer *lwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif