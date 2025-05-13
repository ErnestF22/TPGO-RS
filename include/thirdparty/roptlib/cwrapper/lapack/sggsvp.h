#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int sggsvp_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, realRopt *a, integer *lda, realRopt *b, integer *ldb, realRopt *tola, realRopt *tolb, integer *k, integer *l, realRopt *u, integer *ldu, realRopt *v, integer *ldv, realRopt *q, integer *ldq, integer *iwork, realRopt *tau, realRopt *work, integer *info);

#ifdef __cplusplus
}
#endif