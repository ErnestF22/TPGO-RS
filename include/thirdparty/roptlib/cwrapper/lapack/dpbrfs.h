#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dpbrfs_(char *uplo, integer *n, integer *kd, integer *nrhs, doublerealRopt *ab, integer *ldab, doublerealRopt *afb, integer *ldafb, doublerealRopt *b, integer *ldb, doublerealRopt *x, integer *ldx, doublerealRopt *ferr, doublerealRopt *berr, doublerealRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif