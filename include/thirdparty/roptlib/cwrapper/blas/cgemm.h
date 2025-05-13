#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int cgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, complexRopt *alpha, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, complexRopt *beta, complexRopt *c__, integer *ldc);

#ifdef __cplusplus
}
#endif