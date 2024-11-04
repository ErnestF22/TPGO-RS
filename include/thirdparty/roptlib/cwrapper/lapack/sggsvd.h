#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sggsvd_(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k, integer *l, realRopt *a, integer *lda, realRopt *b, integer *ldb, realRopt *alpha, realRopt *beta, realRopt *u, integer *ldu, realRopt *v, integer *ldv, realRopt *q, integer *ldq, realRopt *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif