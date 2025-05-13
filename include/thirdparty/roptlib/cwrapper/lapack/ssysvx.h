#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int ssysvx_(char *fact, char *uplo, integer *n, integer *nrhs, realRopt *a, integer *lda, realRopt *af, integer *ldaf, integer *ipiv, realRopt *b, integer *ldb, realRopt *x, integer *ldx, realRopt *rcond, realRopt *ferr, realRopt *berr, realRopt *work, integer *lwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif