#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int zggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *alpha, doublecomplexRopt *beta, doublecomplexRopt *vl, integer *ldvl, doublecomplexRopt *vr, integer *ldvr, integer *ilo, integer *ihi, doublerealRopt *lscale, doublerealRopt *rscale, doublerealRopt *abnrm, doublerealRopt *bbnrm, doublerealRopt *rconde, doublerealRopt *rcondv, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *iwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif