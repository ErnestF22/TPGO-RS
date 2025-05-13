#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int ctprfs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, complexRopt *ap, complexRopt *b, integer *ldb, complexRopt *x, integer *ldx, realRopt *ferr, realRopt *berr, complexRopt *work, realRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif