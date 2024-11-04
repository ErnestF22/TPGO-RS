#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int ssbevx_(char *jobz, char *range, char *uplo, integer *n, integer *kd, realRopt *ab, integer *ldab, realRopt *q, integer *ldq, realRopt *vl, realRopt *vu, integer *il, integer *iu, realRopt *abstol, integer *m, realRopt *w, realRopt *z__, integer *ldz, realRopt *work, integer *iwork, integer *ifail, integer *info);

#ifdef __cplusplus
}
#endif