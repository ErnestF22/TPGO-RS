#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dsyevr_(char *jobz, char *range, char *uplo, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *vl, doublerealRopt *vu, integer *il, integer *iu, doublerealRopt *abstol, integer *m, doublerealRopt *w, doublerealRopt *z__, integer *ldz, integer *isuppz, doublerealRopt *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

#ifdef __cplusplus
}
#endif