#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dgbsvx_(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublerealRopt *ab, integer *ldab, doublerealRopt *afb, integer *ldafb, integer *ipiv, char *equed, doublerealRopt *r__, doublerealRopt *c__, doublerealRopt *b, integer *ldb, doublerealRopt *x, integer *ldx, doublerealRopt *rcond, doublerealRopt *ferr, doublerealRopt *berr, doublerealRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif