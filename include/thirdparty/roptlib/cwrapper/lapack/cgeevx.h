#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int cgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, complexRopt *a, integer *lda, complexRopt *w, complexRopt *vl, integer *ldvl, complexRopt *vr, integer *ldvr, integer *ilo, integer *ihi, realRopt *scale, realRopt *abnrm, realRopt *rconde, realRopt *rcondv, complexRopt *work, integer *lwork, realRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif