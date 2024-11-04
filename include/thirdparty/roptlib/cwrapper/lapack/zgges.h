#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zgges_(char *jobvsl, char *jobvsr, char *sort, L_fp delctg, integer *n, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, integer *sdim, doublecomplexRopt *alpha, doublecomplexRopt *beta, doublecomplexRopt *vsl, integer *ldvsl, doublecomplexRopt *vsr, integer *ldvsr, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif