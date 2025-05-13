#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int zheevx_(char *jobz, char *range, char *uplo, integer *n, doublecomplexRopt *a, integer *lda, doublerealRopt *vl, doublerealRopt *vu, integer *il, integer *iu, doublerealRopt *abstol, integer *m, doublerealRopt *w, doublecomplexRopt *z__, integer *ldz, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *iwork, integer *ifail, integer *info);

#ifdef __cplusplus
}
#endif