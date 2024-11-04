#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sgges_(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, integer *n, realRopt *a, integer *lda, realRopt *b, integer *ldb, integer *sdim, realRopt *alphar, realRopt *alphai, realRopt *beta, realRopt *vsl, integer *ldvsl, realRopt *vsr, integer *ldvsr, realRopt *work, integer *lwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif