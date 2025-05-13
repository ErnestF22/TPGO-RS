#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int zsprfs_(char *uplo, integer *n, integer *nrhs, doublecomplexRopt *ap, doublecomplexRopt *afp, integer *ipiv, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *x, integer *ldx, doublerealRopt *ferr, doublerealRopt *berr, doublecomplexRopt *work, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif