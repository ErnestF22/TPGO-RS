#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int cgees_(char *jobvs, char *sort, L_fp select, integer *n, complexRopt *a, integer *lda, integer *sdim, complexRopt *w, complexRopt *vs, integer *ldvs, complexRopt *work, integer *lwork, realRopt *rwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif