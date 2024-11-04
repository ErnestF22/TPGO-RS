#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dgees_(char *jobvs, char *sort, L_fp select, integer *n, doublerealRopt *a, integer *lda, integer *sdim, doublerealRopt *wr, doublerealRopt *wi, doublerealRopt *vs, integer *ldvs, doublerealRopt *work, integer *lwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif