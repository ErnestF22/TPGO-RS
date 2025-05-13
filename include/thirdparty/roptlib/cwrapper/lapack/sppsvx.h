#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int sppsvx_(char *fact, char *uplo, integer *n, integer *nrhs, realRopt *ap, realRopt *afp, char *equed, realRopt *s, realRopt *b, integer *ldb, realRopt *x, integer *ldx, realRopt *rcond, realRopt *ferr, realRopt *berr, realRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif