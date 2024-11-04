#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgeesx_(char *jobvs, char *sort, L_fp select, char *sense, integer *n, complexRopt *a, integer *lda, integer *sdim, complexRopt *w, complexRopt *vs, integer *ldvs, realRopt *rconde, realRopt *rcondv, complexRopt *work, integer *lwork, realRopt *rwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif