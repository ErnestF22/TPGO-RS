#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dggsvd_(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k, integer *l, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, doublerealRopt *alpha, doublerealRopt *beta, doublerealRopt *u, integer *ldu, doublerealRopt *v, integer *ldv, doublerealRopt *q, integer *ldq, doublerealRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif