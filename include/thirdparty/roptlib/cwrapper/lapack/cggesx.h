#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, integer *n, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, integer *sdim, complexRopt *alpha, complexRopt *beta, complexRopt *vsl, integer *ldvsl, complexRopt *vsr, integer *ldvsr, realRopt *rconde, realRopt *rcondv, complexRopt *work, integer *lwork, realRopt *rwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif