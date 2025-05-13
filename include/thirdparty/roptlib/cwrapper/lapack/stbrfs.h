#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int stbrfs_(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs, realRopt *ab, integer *ldab, realRopt *b, integer *ldb, realRopt *x, integer *ldx, realRopt *ferr, realRopt *berr, realRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif