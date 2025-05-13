#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int cgesvd_(char *jobu, char *jobvt, integer *m, integer *n, complexRopt *a, integer *lda, realRopt *s, complexRopt *u, integer *ldu, complexRopt *vt, integer *ldvt, complexRopt *work, integer *lwork, realRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif