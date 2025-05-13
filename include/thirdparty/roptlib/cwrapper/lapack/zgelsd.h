#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int zgelsd_(integer *m, integer *n, integer *nrhs, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, doublerealRopt *s, doublerealRopt *rcond, integer *rank, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif