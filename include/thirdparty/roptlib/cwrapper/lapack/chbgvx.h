#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int chbgvx_(char *jobz, char *range, char *uplo, integer *n, integer *ka, integer *kb, complexRopt *ab, integer *ldab, complexRopt *bb, integer *ldbb, complexRopt *q, integer *ldq, realRopt *vl, realRopt *vu, integer *il, integer *iu, realRopt *abstol, integer *m, realRopt *w, complexRopt *z__, integer *ldz, complexRopt *work, realRopt *rwork, integer *iwork, integer *ifail, integer *info);

#ifdef __cplusplus
}
#endif