#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int zsysvx_(char *fact, char *uplo, integer *n, integer *nrhs, doublecomplexRopt *a, integer *lda, doublecomplexRopt *af, integer *ldaf, integer *ipiv, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *x, integer *ldx, doublerealRopt *rcond, doublerealRopt *ferr, doublerealRopt *berr, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif