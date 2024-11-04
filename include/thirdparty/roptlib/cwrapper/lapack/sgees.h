#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sgees_(char *jobvs, char *sort, L_fp select, integer *n, realRopt *a, integer *lda, integer *sdim, realRopt *wr, realRopt *wi, realRopt *vs, integer *ldvs, realRopt *work, integer *lwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif