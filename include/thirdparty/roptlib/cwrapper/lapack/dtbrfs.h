#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dtbrfs_(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs, doublerealRopt *ab, integer *ldab, doublerealRopt *b, integer *ldb, doublerealRopt *x, integer *ldx, doublerealRopt *ferr, doublerealRopt *berr, doublerealRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif