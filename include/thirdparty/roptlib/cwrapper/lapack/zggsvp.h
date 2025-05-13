#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int zggsvp_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, doublerealRopt *tola, doublerealRopt *tolb, integer *k, integer *l, doublecomplexRopt *u, integer *ldu, doublecomplexRopt *v, integer *ldv, doublecomplexRopt *q, integer *ldq, integer *iwork, doublerealRopt *rwork, doublecomplexRopt *tau, doublecomplexRopt *work, integer *info);

#ifdef __cplusplus
}
#endif