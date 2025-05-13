#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int cgbrfs_(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, complexRopt *ab, integer *ldab, complexRopt *afb, integer *ldafb, integer *ipiv, complexRopt *b, integer *ldb, complexRopt *x, integer *ldx, realRopt *ferr, realRopt *berr, complexRopt *work, realRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif