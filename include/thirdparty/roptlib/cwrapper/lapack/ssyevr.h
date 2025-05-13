#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int ssyevr_(char *jobz, char *range, char *uplo, integer *n, realRopt *a, integer *lda, realRopt *vl, realRopt *vu, integer *il, integer *iu, realRopt *abstol, integer *m, realRopt *w, realRopt *z__, integer *ldz, integer *isuppz, realRopt *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

#ifdef __cplusplus
}
#endif