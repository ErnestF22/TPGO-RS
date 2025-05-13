#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int sgesvx_(char *fact, char *trans, integer *n, integer *nrhs, realRopt *a, integer *lda, realRopt *af, integer *ldaf, integer *ipiv, char *equed, realRopt *r__, realRopt *c__, realRopt *b, integer *ldb, realRopt *x, integer *ldx, realRopt *rcond, realRopt *ferr, realRopt *berr, realRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif