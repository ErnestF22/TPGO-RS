#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zhgeqz_(char *job, char *compq, char *compz, integer *n, integer *ilo, integer *ihi, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *alpha, doublecomplexRopt *beta, doublecomplexRopt *q, integer *ldq, doublecomplexRopt *z__, integer *ldz, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif