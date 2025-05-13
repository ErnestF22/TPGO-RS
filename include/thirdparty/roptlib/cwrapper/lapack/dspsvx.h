#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dspsvx_(char *fact, char *uplo, integer *n, integer *nrhs, doublerealRopt *ap, doublerealRopt *afp, integer *ipiv, doublerealRopt *b, integer *ldb, doublerealRopt *x, integer *ldx, doublerealRopt *rcond, doublerealRopt *ferr, doublerealRopt *berr, doublerealRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif