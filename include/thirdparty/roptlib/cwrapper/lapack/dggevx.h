#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, doublerealRopt *alphar, doublerealRopt *alphai, doublerealRopt *beta, doublerealRopt *vl, integer *ldvl, doublerealRopt *vr, integer *ldvr, integer *ilo, integer *ihi, doublerealRopt *lscale, doublerealRopt *rscale, doublerealRopt *abnrm, doublerealRopt *bbnrm, doublerealRopt *rconde, doublerealRopt *rcondv, doublerealRopt *work, integer *lwork, integer *iwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif