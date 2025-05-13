#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int cgbsvx_(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, complexRopt *ab, integer *ldab, complexRopt *afb, integer *ldafb, integer *ipiv, char *equed, realRopt *r__, realRopt *c__, complexRopt *b, integer *ldb, complexRopt *x, integer *ldx, realRopt *rcond, realRopt *ferr, realRopt *berr, complexRopt *work, realRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif