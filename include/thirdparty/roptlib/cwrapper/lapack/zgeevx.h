#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int zgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, doublecomplexRopt *a, integer *lda, doublecomplexRopt *w, doublecomplexRopt *vl, integer *ldvl, doublecomplexRopt *vr, integer *ldvr, integer *ilo, integer *ihi, doublerealRopt *scale, doublerealRopt *abnrm, doublerealRopt *rconde, doublerealRopt *rcondv, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif