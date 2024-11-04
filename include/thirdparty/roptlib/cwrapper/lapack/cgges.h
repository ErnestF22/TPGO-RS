#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgges_(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, integer *n, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, integer *sdim, complexRopt *alpha, complexRopt *beta, complexRopt *vsl, integer *ldvsl, complexRopt *vsr, integer *ldvsr, complexRopt *work, integer *lwork, realRopt *rwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif