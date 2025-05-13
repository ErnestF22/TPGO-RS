#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dgerfs_(char *trans, integer *n, integer *nrhs, doublerealRopt *a, integer *lda, doublerealRopt *af, integer *ldaf, integer *ipiv, doublerealRopt *b, integer *ldb, doublerealRopt *x, integer *ldx, doublerealRopt *ferr, doublerealRopt *berr, doublerealRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif