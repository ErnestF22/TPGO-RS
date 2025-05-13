#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dgesvx_(char *fact, char *trans, integer *n, integer *nrhs, doublerealRopt *a, integer *lda, doublerealRopt *af, integer *ldaf, integer *ipiv, char *equed, doublerealRopt *r__, doublerealRopt *c__, doublerealRopt *b, integer *ldb, doublerealRopt *x, integer *ldx, doublerealRopt *rcond, doublerealRopt *ferr, doublerealRopt *berr, doublerealRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif