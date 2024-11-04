#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sgeesx_(char *jobvs, char *sort, L_fp select, char *sense, integer *n, realRopt *a, integer *lda, integer *sdim, realRopt *wr, realRopt *wi, realRopt *vs, integer *ldvs, realRopt *rconde, realRopt *rcondv, realRopt *work, integer *lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif