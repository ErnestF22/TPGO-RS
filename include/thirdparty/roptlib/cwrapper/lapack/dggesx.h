#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp delctg, char *sense, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, integer *sdim, doublerealRopt *alphar, doublerealRopt *alphai, doublerealRopt *beta, doublerealRopt *vsl, integer *ldvsl, doublerealRopt *vsr, integer *ldvsr, doublerealRopt *rconde, doublerealRopt *rcondv, doublerealRopt *work, integer *lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif