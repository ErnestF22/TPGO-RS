#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgelsd_(integer *m, integer *n, integer *nrhs, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, realRopt *s, realRopt *rcond, integer *rank, complexRopt *work, integer *lwork, realRopt *rwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif