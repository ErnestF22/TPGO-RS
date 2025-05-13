#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *wr, doublerealRopt *wi, doublerealRopt *vl, integer *ldvl, doublerealRopt *vr, integer *ldvr, integer *ilo, integer *ihi, doublerealRopt *scale, doublerealRopt *abnrm, doublerealRopt *rconde, doublerealRopt *rcondv, doublerealRopt *work, integer *lwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif