#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dgelsd_(integer *m, integer *n, integer *nrhs, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, doublerealRopt *s, doublerealRopt *rcond, integer *rank, doublerealRopt *work, integer *lwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif