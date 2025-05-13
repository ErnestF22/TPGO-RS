#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int zgesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublecomplexRopt *a, integer *lda, doublerealRopt *s, doublecomplexRopt *u, integer *ldu, doublecomplexRopt *vt, integer *ldvt, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif