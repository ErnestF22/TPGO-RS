#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cggsvp_(char *jobu, char *jobv, char *jobq, integer *m, integer *p, integer *n, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, realRopt *tola, realRopt *tolb, integer *k, integer *l, complexRopt *u, integer *ldu, complexRopt *v, integer *ldv, complexRopt *q, integer *ldq, integer *iwork, realRopt *rwork, complexRopt *tau, complexRopt *work, integer *info);

#ifdef __cplusplus
}
#endif