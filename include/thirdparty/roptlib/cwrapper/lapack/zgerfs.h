#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zgerfs_(char *trans, integer *n, integer *nrhs, doublecomplexRopt *a, integer *lda, doublecomplexRopt *af, integer *ldaf, integer *ipiv, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *x, integer *ldx, doublerealRopt *ferr, doublerealRopt *berr, doublecomplexRopt *work, doublerealRopt *rwork, integer *info);

#ifdef __cplusplus
}
#endif