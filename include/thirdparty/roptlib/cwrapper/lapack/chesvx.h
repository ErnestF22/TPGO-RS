#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int chesvx_(char *fact, char *uplo, integer *n, integer *nrhs, complexRopt *a, integer *lda, complexRopt *af, integer *ldaf, integer *ipiv, complexRopt *b, integer *ldb, complexRopt *x, integer *ldx, realRopt *rcond, realRopt *ferr, realRopt *berr, complexRopt *work, integer *lwork, realRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif