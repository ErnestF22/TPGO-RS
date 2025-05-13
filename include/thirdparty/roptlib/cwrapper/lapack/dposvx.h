#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dposvx_(char *fact, char *uplo, integer *n, integer *nrhs, doublerealRopt *a, integer *lda, doublerealRopt *af, integer *ldaf, char *equed, doublerealRopt *s, doublerealRopt *b, integer *ldb, doublerealRopt *x, integer *ldx, doublerealRopt *rcond, doublerealRopt *ferr, doublerealRopt *berr, doublerealRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif