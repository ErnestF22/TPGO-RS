#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgbmv_(char *trans, integer *m, integer *n, integer *kl, integer *ku, complexRopt *alpha, complexRopt *a, integer *lda, complexRopt *x, integer *incx, complexRopt *beta, complexRopt *y, integer *incy);

#ifdef __cplusplus
}
#endif