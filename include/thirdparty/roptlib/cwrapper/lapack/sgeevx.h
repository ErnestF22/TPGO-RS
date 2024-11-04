#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, realRopt *a, integer *lda, realRopt *wr, realRopt *wi, realRopt *vl, integer *ldvl, realRopt *vr, integer *ldvr, integer *ilo, integer *ihi, realRopt *scale, realRopt *abnrm, realRopt *rconde, realRopt *rcondv, realRopt *work, integer *lwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif