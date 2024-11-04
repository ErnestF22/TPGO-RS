#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, realRopt *a, integer *lda, realRopt *b, integer *ldb, realRopt *alphar, realRopt *alphai, realRopt *beta, realRopt *vl, integer *ldvl, realRopt *vr, integer *ldvr, integer *ilo, integer *ihi, realRopt *lscale, realRopt *rscale, realRopt *abnrm, realRopt *bbnrm, realRopt *rconde, realRopt *rcondv, realRopt *work, integer *lwork, integer *iwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif