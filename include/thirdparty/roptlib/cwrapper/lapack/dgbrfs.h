#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dgbrfs_(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublerealRopt *ab, integer *ldab, doublerealRopt *afb, integer *ldafb, integer *ipiv, doublerealRopt *b, integer *ldb, doublerealRopt *x, integer *ldx, doublerealRopt *ferr, doublerealRopt *berr, doublerealRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif