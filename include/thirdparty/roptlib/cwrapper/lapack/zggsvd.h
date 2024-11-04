#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zggsvd_(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k, integer *l, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, doublerealRopt *alpha, doublerealRopt *beta, doublecomplexRopt *u, integer *ldu, doublecomplexRopt *v, integer *ldv, doublecomplexRopt *q, integer *ldq, doublecomplexRopt *work, doublerealRopt *rwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif