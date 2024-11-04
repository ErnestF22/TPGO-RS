#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgbtrs_(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, complexRopt *ab, integer *ldab, integer *ipiv, complexRopt *b, integer *ldb, integer *info);

#ifdef __cplusplus
}
#endif