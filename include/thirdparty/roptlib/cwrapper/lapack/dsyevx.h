#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dsyevx_(char *jobz, char *range, char *uplo, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *vl, doublerealRopt *vu, integer *il, integer *iu, doublerealRopt *abstol, integer *m, doublerealRopt *w, doublerealRopt *z__, integer *ldz, doublerealRopt *work, integer *lwork, integer *iwork, integer *ifail, integer *info);

#ifdef __cplusplus
}
#endif