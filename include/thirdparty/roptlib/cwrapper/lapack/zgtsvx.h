#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int zgtsvx_(char *fact, char *trans, integer *n, integer *nrhs, doublecomplexRopt *dl, doublecomplexRopt *d__, doublecomplexRopt *du, doublecomplexRopt *dlf, doublecomplexRopt *df, doublecomplexRopt *duf, doublecomplexRopt *du2, integer *ipiv, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *x, integer *ldx, doublerealRopt *rcond, doublerealRopt *ferr, doublerealRopt *berr, doublecomplexRopt *work, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif