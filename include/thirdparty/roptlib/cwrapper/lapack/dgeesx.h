#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dgeesx_(char *jobvs, char *sort, L_fp select, char *sense, integer *n, doublerealRopt *a, integer *lda, integer *sdim, doublerealRopt *wr, doublerealRopt *wi, doublerealRopt *vs, integer *ldvs, doublerealRopt *rconde, doublerealRopt *rcondv, doublerealRopt *work, integer *lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif