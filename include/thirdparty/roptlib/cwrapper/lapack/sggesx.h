#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int sggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, integer *n, realRopt *a, integer *lda, realRopt *b, integer *ldb, integer *sdim, realRopt *alphar, realRopt *alphai, realRopt *beta, realRopt *vsl, integer *ldvsl, realRopt *vsr, integer *ldvsr, realRopt *rconde, realRopt *rcondv, realRopt *work, integer *lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif