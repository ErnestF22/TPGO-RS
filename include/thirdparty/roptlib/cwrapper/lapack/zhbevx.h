#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zhbevx_(char *jobz, char *range, char *uplo, integer *n, integer *kd, doublecomplexRopt *ab, integer *ldab, doublecomplexRopt *q, integer *ldq, doublerealRopt *vl, doublerealRopt *vu, integer *il, integer *iu, doublerealRopt *abstol, integer *m, doublerealRopt *w, doublecomplexRopt *z__, integer *ldz, doublecomplexRopt *work, doublerealRopt *rwork, integer *iwork, integer *ifail, integer *info);

#ifdef __cplusplus
}
#endif