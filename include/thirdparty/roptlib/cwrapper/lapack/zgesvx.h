#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int zgesvx_(char *fact, char *trans, integer *n, integer *nrhs, doublecomplexRopt *a, integer *lda, doublecomplexRopt *af, integer *ldaf, integer *ipiv, char *equed, doublerealRopt *r__, doublerealRopt *c__, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *x, integer *ldx, doublerealRopt *rcond, doublerealRopt *ferr, doublerealRopt *berr, doublecomplexRopt *work, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif