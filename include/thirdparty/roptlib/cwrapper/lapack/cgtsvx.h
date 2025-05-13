#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int cgtsvx_(char *fact, char *trans, integer *n, integer *nrhs, complexRopt *dl, complexRopt *d__, complexRopt *du, complexRopt *dlf, complexRopt *df, complexRopt *duf, complexRopt *du2, integer *ipiv, complexRopt *b, integer *ldb, complexRopt *x, integer *ldx, realRopt *rcond, realRopt *ferr, realRopt *berr, complexRopt *work, realRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif