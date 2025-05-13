#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int zggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp delctg, char *sense, integer *n, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, integer *sdim, doublecomplexRopt *alpha, doublecomplexRopt *beta, doublecomplexRopt *vsl, integer *ldvsl, doublecomplexRopt *vsr, integer *ldvsr, doublerealRopt *rconde, doublerealRopt *rcondv, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif