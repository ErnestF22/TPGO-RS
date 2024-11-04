#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int chpgvx_(integer *itype, char *jobz, char *range, char *uplo, integer *n, complexRopt *ap, complexRopt *bp, realRopt *vl, realRopt *vu, integer *il, integer *iu, realRopt *abstol, integer *m, realRopt *w, complexRopt *z__, integer *ldz, complexRopt *work, realRopt *rwork, integer *iwork, integer *ifail, integer *info);

#ifdef __cplusplus
}
#endif