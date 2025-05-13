#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int cggsvd_(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k, integer *l, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, realRopt *alpha, realRopt *beta, complexRopt *u, integer *ldu, complexRopt *v, integer *ldv, complexRopt *q, integer *ldq, complexRopt *work, realRopt *rwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif