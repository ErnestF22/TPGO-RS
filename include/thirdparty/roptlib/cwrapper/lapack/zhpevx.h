#ifdef __cplusplus
extern "C"
{
#endif

#include "f2c.h"

    int zhpevx_(char *jobz, char *range, char *uplo, integer *n, doublecomplexRopt *ap, doublerealRopt *vl, doublerealRopt *vu, integer *il, integer *iu, doublerealRopt *abstol, integer *m, doublerealRopt *w, doublecomplexRopt *z__, integer *ldz, doublecomplexRopt *work, doublerealRopt *rwork, integer *iwork, integer *ifail, integer *info);

#ifdef __cplusplus
}
#endif