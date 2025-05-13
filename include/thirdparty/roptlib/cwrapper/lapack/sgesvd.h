#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int sgesvd_(char *jobu, char *jobvt, integer *m, integer *n, realRopt *a, integer *lda, realRopt *s, realRopt *u, integer *ldu, realRopt *vt, integer *ldvt, realRopt *work, integer *lwork, integer *info);

#ifdef __cplusplus
}
#endif