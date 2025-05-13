#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int dggsvp_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, doublerealRopt *tola, doublerealRopt *tolb, integer *k, integer *l, doublerealRopt *u, integer *ldu, doublerealRopt *v, integer *ldv, doublerealRopt *q, integer *ldq, integer *iwork, doublerealRopt *tau, doublerealRopt *work, integer *info);

#ifdef __cplusplus
}
#endif