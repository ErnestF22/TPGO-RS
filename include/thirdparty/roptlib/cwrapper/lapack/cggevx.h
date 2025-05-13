#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int cggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, complexRopt *alpha, complexRopt *beta, complexRopt *vl, integer *ldvl, complexRopt *vr, integer *ldvr, integer *ilo, integer *ihi, realRopt *lscale, realRopt *rscale, realRopt *abnrm, realRopt *bbnrm, realRopt *rconde, realRopt *rcondv, complexRopt *work, integer *lwork, realRopt *rwork, integer *iwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif