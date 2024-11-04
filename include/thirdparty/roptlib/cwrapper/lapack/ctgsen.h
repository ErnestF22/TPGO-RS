#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int ctgsen_(integer *ijob, logical *wantq, logical *wantz, logical *select, integer *n, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, complexRopt *alpha, complexRopt *beta, complexRopt *q, integer *ldq, complexRopt *z__, integer *ldz, integer *m, realRopt *pl, realRopt *pr, realRopt *dif, complexRopt *work, integer *lwork, integer *iwork, integer *liwork, integer *info);

#ifdef __cplusplus
}
#endif