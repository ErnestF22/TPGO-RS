#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgbcon_(char *norm, integer *n, integer *kl, integer *ku, complexRopt *ab, integer *ldab, integer *ipiv, realRopt *anorm, realRopt *rcond, complexRopt *work, realRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif