#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int chbevd_(char *jobz, char *uplo, integer *n, integer *kd, complexRopt *ab, integer *ldab, realRopt *w, complexRopt *z__, integer *ldz, complexRopt *work, integer *lwork, realRopt *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);

#ifdef __cplusplus
}
#endif