/*
Some blas and lapack functions are not easy to use. This class defines wrapper functions
which can be used to solve some problems or decomposition, e.g., SVD decomposition, Sylevster equation.
More functions will be added in this class.

-----WH
*/

#include "others_BlasLapackCppWrapper.h"

/*Define the namespace */
namespace ROPTLIB {
	//#include <dasum.h>
	//#include <daxpy.h>
	void axpy_(integer *n, complexRopt *ca, complexRopt *cx, integer *incx, complexRopt *cy, integer *incy)
	{
#ifndef MATLAB_MEX_FILE
		caxpy_(n, ca, cx, incx, cy, incy);
#else
		caxpy_(n, (float *)ca, (float *)cx, incx, (float *)cy, incy);
#endif
	};
	void axpy_(integer *n, doublerealRopt *da, doublerealRopt *dx, integer *incx, doublerealRopt *dy, integer *incy)
	{
		daxpy_(n, da, dx, incx, dy, incy);
	};
	void axpy_(integer *n, realRopt *sa, realRopt *sx, integer *incx, realRopt *sy, integer *incy)
	{
		saxpy_(n, sa, sx, incx, sy, incy);
	};
	void axpy_(integer *n, doublecomplexRopt *za, doublecomplexRopt *zx, integer *incx, doublecomplexRopt *zy, integer *incy)
	{
#ifndef MATLAB_MEX_FILE
		zaxpy_(n, za, zx, incx, zy, incy);
#else
		zaxpy_(n, (double *)za, (double *)zx, incx, (double *)zy, incy);
#endif
	};

	//#include <dcabs1.h>
	//#include <dcopy.h>
	void copy_(integer *n, complexRopt *cx, integer *incx, complexRopt *cy, integer *incy)
	{
#ifndef MATLAB_MEX_FILE
		ccopy_(n, cx, incx, cy, incy);
#else
		ccopy_(n, (float *)cx, incx, (float *)cy, incy);
#endif
	};
	void copy_(integer *n, doublerealRopt *cx, integer *incx, doublerealRopt *cy, integer *incy)
	{
		dcopy_(n, cx, incx, cy, incy);
	};
	void copy_(integer *n, realRopt *cx, integer *incx, realRopt *cy, integer *incy)
	{
		scopy_(n, cx, incx, cy, incy);
	};
	void copy_(integer *n, doublecomplexRopt *cx, integer *incx, doublecomplexRopt *cy, integer *incy)
	{
#ifndef MATLAB_MEX_FILE
		zcopy_(n, cx, incx, cy, incy);
#else
		zcopy_(n, (double *)cx, incx, (double *)cy, incy);
#endif
	};

	//#include <ddot.h>
	doublerealRopt dot_(integer *n, doublerealRopt *dx, integer *incx, doublerealRopt *dy, integer *incy)
	{
		return ddot_(n, dx, incx, dy, incy);
	}
	realRopt dot_(integer *n, realRopt *dx, integer *incx, realRopt *dy, integer *incy)
	{
		return sdot_(n, dx, incx, dy, incy);
	}
	void dotu_(complexRopt * ret_val, integer *n, complexRopt *cx, integer *incx, complexRopt *cy, integer *incy)
	{
#ifndef MATLAB_MEX_FILE
		return cdotu_(ret_val, n, cx, incx, cy, incy);
#else
		ret_val[0] = cdotu_(n, (float *)cx, incx, (float *)cy, incy);
#endif
	};
	void dotu_(doublecomplexRopt * ret_val, integer *n, doublecomplexRopt *cx, integer *incx, doublecomplexRopt *cy, integer *incy)
	{
#ifndef MATLAB_MEX_FILE
		return zdotu_(ret_val, n, cx, incx, cy, incy);
#else
		ret_val[0] = zdotu_(n, (double *)cx, incx, (double *)cy, incy);
#endif
	};
	void dotc_(complexRopt * ret_val, integer *n, complexRopt *cx, integer *incx, complexRopt *cy, integer *incy)
	{
#ifndef MATLAB_MEX_FILE
		return cdotc_(ret_val, n, cx, incx, cy, incy);
#else
		ret_val[0] = cdotc_(n, (float *)cx, incx, (float *)cy, incy);
#endif
	};
	void dotc_(doublecomplexRopt * ret_val, integer *n, doublecomplexRopt *cx, integer *incx, doublecomplexRopt *cy, integer *incy)
	{
#ifndef MATLAB_MEX_FILE
		return zdotc_(ret_val, n, cx, incx, cy, incy);
#else
		ret_val[0] = zdotc_(n, (double *)cx, incx, (double *)cy, incy);
#endif
	};

	//#include <dgbmv.h>
	//#include <dgemm.h>
	void gemm_(char *transa, char *transb, integer *m, integer *n, integer *k, complexRopt *alpha, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, complexRopt *beta, complexRopt *c__, integer *ldc)
	{
#ifndef MATLAB_MEX_FILE
		cgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
#else
		cgemm_(transa, transb, m, n, k, (float *)alpha, (float *)a, lda, (float *)b, ldb, (float *)beta, (float *)c__, ldc);
#endif
	};
	void gemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublerealRopt *alpha, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, doublerealRopt *beta, doublerealRopt *c__, integer *ldc)
	{
		dgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
	};
	void gemm_(char *transa, char *transb, integer *m, integer *n, integer *k, realRopt *alpha, realRopt *a, integer *lda, realRopt *b, integer *ldb, realRopt *beta, realRopt *c__, integer *ldc)
	{
		sgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
	};
	void gemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublecomplexRopt *alpha, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *beta, doublecomplexRopt *c__, integer *ldc)
	{
#ifndef MATLAB_MEX_FILE
		zgemm_(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
#else
		zgemm_(transa, transb, m, n, k, (double *)alpha, (double *)a, lda, (double *)b, ldb, (double *)beta, (double *)c__, ldc);
#endif
	};

	//#include <dgemv.h>
	void gemv_(char *trans, integer *m, integer *n, complexRopt *alpha, complexRopt *a, integer *lda, complexRopt *x, integer *incx, complexRopt *beta, complexRopt *y, integer *incy)
	{
#ifndef MATLAB_MEX_FILE
		cgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
		cgemv_(trans, m, n, (float *)alpha, (float *)a, lda, (float *)x, incx, (float *)beta, (float *)y, incy);
#endif
	};
	void gemv_(char *trans, integer *m, integer *n, doublerealRopt *alpha, doublerealRopt *a, integer *lda, doublerealRopt *x, integer *incx, doublerealRopt *beta, doublerealRopt *y, integer *incy)
	{
		dgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
	};
	void gemv_(char *trans, integer *m, integer *n, realRopt *alpha, realRopt *a, integer *lda, realRopt *x, integer *incx, realRopt *beta, realRopt *y, integer *incy)
	{
		sgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
	};
	void gemv_(char *trans, integer *m, integer *n, doublecomplexRopt *alpha, doublecomplexRopt *a, integer *lda, doublecomplexRopt *x, integer *incx, doublecomplexRopt *beta, doublecomplexRopt *y, integer *incy)
	{
#ifndef MATLAB_MEX_FILE
		zgemv_(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
		zgemv_(trans, m, n, (double *)alpha, (double *)a, lda, (double *)x, incx, (double *)beta, (double *)y, incy);
#endif
	};


	//#include <dger.h>
	void ger_(integer *m, integer *n, doublerealRopt *alpha, doublerealRopt *x, integer *incx, doublerealRopt *y, integer *incy, doublerealRopt *a, integer *lda)
	{
		dger_(m, n, alpha, x, incx, y, incy, a, lda);
	};
	void ger_(integer *m, integer *n, realRopt *alpha, realRopt *x, integer *incx, realRopt *y, integer *incy, realRopt *a, integer *lda)
	{
		sger_(m, n, alpha, x, incx, y, incy, a, lda);
	};

	void gerc_(integer *m, integer *n, complexRopt *alpha, complexRopt *x, integer *incx, complexRopt *y, integer *incy, complexRopt *a, integer *lda)
	{
#ifndef MATLAB_MEX_FILE
		cgerc_(m, n, alpha, x, incx, y, incy, a, lda);
#else
		cgerc_(m, n, (float *)alpha, (float *)x, incx, (float *)y, incy, (float *)a, lda);
#endif
	};
	void gerc_(integer *m, integer *n, doublecomplexRopt *alpha, doublecomplexRopt *x, integer *incx, doublecomplexRopt *y, integer *incy, doublecomplexRopt *a, integer *lda)
	{
#ifndef MATLAB_MEX_FILE
		zgerc_(m, n, alpha, x, incx, y, incy, a, lda);
#else
		zgerc_(m, n, (double *)alpha, (double *)x, incx, (double *)y, incy, (double *)a, lda);
#endif
	};

	void ger_(integer *m, integer *n, complexRopt *alpha, complexRopt *x, integer *incx, complexRopt *y, integer *incy, complexRopt *a, integer *lda)
	{
#ifndef MATLAB_MEX_FILE
		cgeru_(m, n, alpha, x, incx, y, incy, a, lda);
#else
		cgeru_(m, n, (float *)alpha, (float *)x, incx, (float *)y, incy, (float *)a, lda);
#endif
	};
	void ger_(integer *m, integer *n, doublecomplexRopt *alpha, doublecomplexRopt *x, integer *incx, doublecomplexRopt *y, integer *incy, doublecomplexRopt *a, integer *lda)
	{
#ifndef MATLAB_MEX_FILE
		zgeru_(m, n, alpha, x, incx, y, incy, a, lda);
#else
		zgeru_(m, n, (double *)alpha, (double *)x, incx, (double *)y, incy, (double *)a, lda);
#endif
	};

	//#include <dnrm2.h>
	doublerealRopt nrm2_(integer *n, doublerealRopt *x, integer *incx)
	{
		return dnrm2_(n, x, incx);
	};
	realRopt nrm2_(integer *n, realRopt *x, integer *incx)
	{
		return snrm2_(n, x, incx);
	};

	//#include <drot.h>
	//#include <drotg.h>
	//#include <dsbmv.h>
	//#include <dscal.h>
	void scal_(integer *n, complexRopt *ca, complexRopt *cx, integer *incx)
	{
#ifndef MATLAB_MEX_FILE
		cscal_(n, ca, cx, incx);
#else
		cscal_(n, (float *)ca, (float *)cx, incx);
#endif
	};
	void scal_(integer *n, doublerealRopt *ca, doublerealRopt *cx, integer *incx)
	{
		dscal_(n, ca, cx, incx);
	};
	void scal_(integer *n, realRopt *ca, realRopt *cx, integer *incx)
	{
		sscal_(n, ca, cx, incx);
	};
	void scal_(integer *n, doublecomplexRopt *ca, doublecomplexRopt *cx, integer *incx)
	{
#ifndef MATLAB_MEX_FILE
		zscal_(n, ca, cx, incx);
#else
		zscal_(n, (double *)ca, (double *)cx, incx);
#endif
	};

	//#include <dspmv.h>
	//#include <dspr.h>
	//#include <dspr2.h>
	//#include <dswap.h>
	//#include <dsymm.h>
	//#include <dsymv.h>
	void symv_(char *uplo, integer *n, complexRopt *alpha, complexRopt *a, integer *lda, complexRopt *x, integer *incx, complexRopt *beta, complexRopt *y, integer *incy)
	{
#ifndef MATLAB_MEX_FILE
		csymv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
#else
		csymv_(uplo, n, (float *)alpha, (float *)a, lda, (float *)x, incx, (float *)beta, (float *)y, incy);
#endif
	};
	void symv_(char *uplo, integer *n, doublerealRopt *alpha, doublerealRopt *a, integer *lda, doublerealRopt *x, integer *incx, doublerealRopt *beta, doublerealRopt *y, integer *incy)
	{
		dsymv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
	};
	void symv_(char *uplo, integer *n, realRopt *alpha, realRopt *a, integer *lda, realRopt *x, integer *incx, realRopt *beta, realRopt *y, integer *incy)
	{
		ssymv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
	};
	void symv_(char *uplo, integer *n, doublecomplexRopt *alpha, doublecomplexRopt *a, integer *lda, doublecomplexRopt *x, integer *incx, doublecomplexRopt *beta, doublecomplexRopt *y, integer *incy)
	{
#ifndef MATLAB_MEX_FILE
		zsymv_(uplo, n, alpha, a, lda, x, incx, beta, y, incy);
#else
		zsymv_(uplo, n, (double *)alpha, (double *)a, lda, (double *)x, incx, (double *)beta, (double *)y, incy);
#endif
	};

	//#include <dsyr.h>
	//#include <dsyr2.h>
	//#include <dsyr2k.h>
	//#include <dsyrk.h>
	//#include <dtbmv.h>
	//#include <dtbsv.h>
	//#include <dtpmv.h>
	//#include <dtpsv.h>
	//#include <dtrmm.h>
	//#include <dtrmv.h>
	//#include <dtrsm.h>
	void trsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublerealRopt *alpha, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb)
	{
		dtrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
	};
	void trsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, realRopt *alpha, realRopt *a, integer *lda, realRopt *b, integer *ldb)
	{
		strsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
	};
    void trsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublecomplexRopt *alpha, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb)
    {
#ifndef MATLAB_MEX_FILE
        ztrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
#else
        ztrsm_(side, uplo, transa, diag, m, n, (double *) alpha, (double *) a, lda, (double *) b, ldb);
#endif
    };
    void trsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, complexRopt *alpha, complexRopt *a, integer *lda, complexRopt *b, integer *ldb)
    {
#ifndef MATLAB_MEX_FILE
        ctrsm_(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
#else
        ctrsm_(side, uplo, transa, diag, m, n, (float *) alpha, (float *) a, lda, (float *) b, ldb);
#endif
    };
	//#include <dtrsv.h>
	//#include <dzasum.h>
	//#include <dznrm2.h>



	//#include <dbdsdc.h>
	//#include <dbdsqr.h>
	//#include <ddisna.h>
	//#include <dgbbrd.h>
	//#include <dgbcon.h>
	//#include <dgbequ.h>
	//#include <dgbrfs.h>
	//#include <dgbsv.h>
	//#include <dgbsvx.h>
	//#include <dgbtf2.h>
	//#include <dgbtrf.h>
	//#include <dgbtrs.h>
	//#include <dgebak.h>
	//#include <dgebal.h>
	//#include <dgebd2.h>
	//#include <dgebrd.h>
	//#include <dgecon.h>
	//#include <dgeequ.h>
	//#include <dgees.h>
	void gees_(char *jobvs, char *sort, L_fp select, integer *n, complexRopt *a, integer *lda, integer *sdim, complexRopt *w, complexRopt *vs, integer *ldvs, complexRopt *work, integer *lwork, realRopt *rwork, logical *bwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cgees_(jobvs, sort, select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, info);
#else
		cgees_(jobvs, sort, select, n, (float *)a, lda, sdim, (float *)w, (float *)vs, ldvs, (float *)work, lwork, rwork, (integer *)bwork, info);
#endif
	};
	void gees_(char *jobvs, char *sort, L_fp select, integer *n, doublerealRopt *a, integer *lda, integer *sdim, doublerealRopt *wr, doublerealRopt *wi, doublerealRopt *vs, integer *ldvs, doublerealRopt *work, integer *lwork, logical *bwork, integer *info)
	{
		dgees_(jobvs, sort, select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, info);
	};
	void gees_(char *jobvs, char *sort, L_fp select, integer *n, realRopt *a, integer *lda, integer *sdim, realRopt *wr, realRopt *wi, realRopt *vs, integer *ldvs, realRopt *work, integer *lwork, logical *bwork, integer *info)
	{
		sgees_(jobvs, sort, select, n, a, lda, sdim, wr, wi, vs, ldvs, work, lwork, bwork, info);
	};
	void gees_(char *jobvs, char *sort, L_fp select, integer *n, doublecomplexRopt *a, integer *lda, integer *sdim, doublecomplexRopt *w, doublecomplexRopt *vs, integer *ldvs, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, logical *bwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zgees_(jobvs, sort, select, n, a, lda, sdim, w, vs, ldvs, work, lwork, rwork, bwork, info);
#else
		zgees_(jobvs, sort, select, n, (double *)a, lda, sdim, (double *)w, (double *)vs, ldvs, (double *)work, lwork, rwork, (integer *)bwork, info);
#endif
	};

	//#include <dgeesx.h>
	//#include <dgeev.h>
	//#include <dgeevx.h>
	//#include <dgegs.h>
	//#include <dgegv.h>
	//#include <dgehd2.h>
	//#include <dgehrd.h>
	//#include <dgelq2.h>
	//#include <dgelqf.h>
	//#include <dgels.h>
	//#include <dgelsd.h>
	//#include <dgelss.h>
	//#include <dgelsx.h>
	//#include <dgelsy.h>
	//#include <dgeql2.h>
	//#include <dgeqlf.h>
	//#include <dgeqp3.h>
	void geqp3_(integer *m, integer *n, complexRopt *a, integer *lda, integer *jpvt, complexRopt *tau, complexRopt *work, integer *lwork, realRopt *rwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
#else
		cgeqp3_(m, n, (float *)a, lda, jpvt, (float *)tau, (float *)work, lwork, rwork, info);
#endif
	};
	void geqp3_(integer *m, integer *n, doublerealRopt *a, integer *lda, integer *jpvt, doublerealRopt *tau, doublerealRopt *work, integer *lwork, integer *info)
	{
		dgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, info);
	};
	void geqp3_(integer *m, integer *n, realRopt *a, integer *lda, integer *jpvt, realRopt *tau, realRopt *work, integer *lwork, integer *info)
	{
		sgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, info);
	};
	void geqp3_(integer *m, integer *n, doublecomplexRopt *a, integer *lda, integer *jpvt, doublecomplexRopt *tau, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
#else
		zgeqp3_(m, n, (double *)a, lda, jpvt, (double *)tau, (double *)work, lwork, rwork, info);
#endif
	};
	//#include <dgeqpf.h>
	//#include <dgeqr2.h>
	//#include <dgeqrf.h>
	//#include <dgerfs.h>
	//#include <dgerq2.h>
	//#include <dgerqf.h>
	//#include <dgesc2.h>
	//#include <dgesdd.h>
	void gesdd_(char *jobz, integer *m, integer *n, complexRopt *a, integer *lda, realRopt *s, complexRopt *u, integer *ldu, complexRopt *vt, integer *ldvt, complexRopt *work, integer *lwork, realRopt *rwork, integer *iwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
#else
		cgesdd_(jobz, m, n, (realRopt *)a, lda, s, (realRopt *)u, ldu, (realRopt *)vt, ldvt, (realRopt *)work, lwork, rwork, iwork, info);
#endif
	}
	void gesdd_(char *jobz, integer *m, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *s, doublerealRopt *u, integer *ldu, doublerealRopt *vt, integer *ldvt, doublerealRopt *work, integer *lwork, integer *iwork, integer *info)
	{
		dgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
	}
	void gesdd_(char *jobz, integer *m, integer *n, realRopt *a, integer *lda, realRopt *s, realRopt *u, integer *ldu, realRopt *vt, integer *ldvt, realRopt *work, integer *lwork, integer *iwork, integer *info)
	{
		sgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
	}
	void gesdd_(char *jobz, integer *m, integer *n, doublecomplexRopt *a, integer *lda, doublerealRopt *s, doublecomplexRopt *u, integer *ldu, doublecomplexRopt *vt, integer *ldvt, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *iwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
#else
		zgesdd_(jobz, m, n, (doublerealRopt *)a, lda, s, (doublerealRopt *)u, ldu, (doublerealRopt *)vt, ldvt, (doublerealRopt *)work, lwork, rwork, iwork, info);
#endif
	}
	//#include <dgesv.h>
	void gesv_(integer *n, integer *nrhs, complexRopt *a, integer *lda, integer *ipiv, complexRopt *b, integer *ldb, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
#else
		cgesv_(n, nrhs, (realRopt *)a, lda, ipiv, (realRopt *)b, ldb, info);
#endif
	};
	void gesv_(integer *n, integer *nrhs, doublerealRopt *a, integer *lda, integer *ipiv, doublerealRopt *b, integer *ldb, integer *info)
	{
		dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
	};
	void gesv_(integer *n, integer *nrhs, realRopt *a, integer *lda, integer *ipiv, realRopt *b, integer *ldb, integer *info)
	{
		sgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
	};
	void gesv_(integer *n, integer *nrhs, doublecomplexRopt *a, integer *lda, integer *ipiv, doublecomplexRopt *b, integer *ldb, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
#else
		zgesv_(n, nrhs, (doublerealRopt *)a, lda, ipiv, (doublerealRopt *)b, ldb, info);
#endif
	};

	//#include <dgesvd.h>
	void gesvd_(char *jobu, char *jobvt, integer *m, integer *n, complexRopt *a, integer *lda, realRopt *s, complexRopt *u, integer *ldu, complexRopt *vt, integer *ldvt, complexRopt *work, integer *lwork, realRopt *rwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
#else
		cgesvd_(jobu, jobvt, m, n, (float *)a, lda, s, (float *)u, ldu, (float *)vt, ldvt, (float *)work, lwork, rwork, info);
#endif
	};
	void gesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *s, doublerealRopt *u, integer *ldu, doublerealRopt *vt, integer *ldvt, doublerealRopt *work, integer *lwork, integer *info)
	{
		dgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
	};
	void gesvd_(char *jobu, char *jobvt, integer *m, integer *n, realRopt *a, integer *lda, realRopt *s, realRopt *u, integer *ldu, realRopt *vt, integer *ldvt, realRopt *work, integer *lwork, integer *info)
	{
		sgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
	};
	void gesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublecomplexRopt *a, integer *lda, doublerealRopt *s, doublecomplexRopt *u, integer *ldu, doublecomplexRopt *vt, integer *ldvt, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zgesvd_(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
#else
		zgesvd_(jobu, jobvt, m, n, (double *)a, lda, s, (double *)u, ldu, (double *)vt, ldvt, (double *)work, lwork, rwork, info);
#endif
	};

	//#include <dgesvx.h>
	//#include <dgetc2.h>
	//#include <dgetf2.h>
	//#include <dgetrf.h>
	void getrf_(integer *m, integer *n, complexRopt *a, integer *lda, integer *ipiv, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cgetrf_(m, n, a, lda, ipiv, info);
#else
		cgetrf_(m, n, (realRopt *)a, lda, ipiv, info);
#endif
	}
	void getrf_(integer *m, integer *n, doublerealRopt *a, integer *lda, integer *ipiv, integer *info)
	{
		dgetrf_(m, n, a, lda, ipiv, info);
	}
	void getrf_(integer *m, integer *n, realRopt *a, integer *lda, integer *ipiv, integer *info)
	{
		sgetrf_(m, n, a, lda, ipiv, info);
	}
	void getrf_(integer *m, integer *n, doublecomplexRopt *a, integer *lda, integer *ipiv, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zgetrf_(m, n, a, lda, ipiv, info);
#else
		zgetrf_(m, n, (doublerealRopt *)a, lda, ipiv, info);
#endif
	}

	//#include <dgetri.h>
	void getri_(integer *n, complexRopt *a, integer *lda, integer *ipiv, complexRopt *work, integer *lwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cgetri_(n, a, lda, ipiv, work, lwork, info);
#else
		cgetri_(n, (realRopt *)a, lda, ipiv, (realRopt *)work, lwork, info);
#endif
	};
	void getri_(integer *n, doublerealRopt *a, integer *lda, integer *ipiv, doublerealRopt *work, integer *lwork, integer *info)
	{
		dgetri_(n, a, lda, ipiv, work, lwork, info);
	};
	void getri_(integer *n, realRopt *a, integer *lda, integer *ipiv, realRopt *work, integer *lwork, integer *info)
	{
		sgetri_(n, a, lda, ipiv, work, lwork, info);
	};
	void getri_(integer *n, doublecomplexRopt *a, integer *lda, integer *ipiv, doublecomplexRopt *work, integer *lwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zgetri_(n, a, lda, ipiv, work, lwork, info);
#else
		zgetri_(n, (doublerealRopt *)a, lda, ipiv, (doublerealRopt *)work, lwork, info);
#endif
	};

	//#include <dgetrs.h>
	void getrs_(char *trans, integer *n, integer *nrhs, complexRopt *a, integer *lda, integer *ipiv, complexRopt *b, integer *ldb, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
#else
		cgetrs_(trans, n, nrhs, (realRopt *)a, lda, ipiv, (realRopt *)b, ldb, info);
#endif
	}
	void getrs_(char *trans, integer *n, integer *nrhs, doublerealRopt *a, integer *lda, integer *ipiv, doublerealRopt *b, integer *ldb, integer *info)
	{
		dgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
	}
	void getrs_(char *trans, integer *n, integer *nrhs, realRopt *a, integer *lda, integer *ipiv, realRopt *b, integer *ldb, integer *info)
	{
		sgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
	}
	void getrs_(char *trans, integer *n, integer *nrhs, doublecomplexRopt *a, integer *lda, integer *ipiv, doublecomplexRopt *b, integer *ldb, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
#else
		zgetrs_(trans, n, nrhs, (doublerealRopt *)a, lda, ipiv, (doublerealRopt *)b, ldb, info);
#endif
	}
	//#include <dggbak.h>
	//#include <dggbal.h>
	//#include <dgges.h>
	//#include <dggesx.h>
	//#include <dggev.h>
	//#include <dggevx.h>
	//#include <dggglm.h>
	//#include <dgghrd.h>
	//#include <dgglse.h>
	//#include <dggqrf.h>
	//#include <dggrqf.h>
	//#include <dggsvd.h>
	//#include <dggsvp.h>
	//#include <dgtcon.h>
	//#include <dgtrfs.h>
	//#include <dgtsv.h>
	//#include <dgtsvx.h>
	//#include <dgttrf.h>
	//#include <dgttrs.h>
	//#include <dgtts2.h>
	//#include <dhgeqz.h>
	//#include <dhsein.h>
	//#include <dhseqr.h>
	//#include <dlabad.h>
	//#include <dlabrd.h>
	//#include <dlacon.h>
	//#include <dlacpy.h>
	//#include <dladiv.h>
	//#include <dlae2.h>
	//#include <dlaebz.h>
	//#include <dlaed0.h>
	//#include <dlaed1.h>
	//#include <dlaed2.h>
	//#include <dlaed3.h>
	//#include <dlaed4.h>
	//#include <dlaed5.h>
	//#include <dlaed6.h>
	//#include <dlaed7.h>
	//#include <dlaed8.h>
	//#include <dlaed9.h>
	//#include <dlaeda.h>
	//#include <dlaein.h>
	//#include <dlaev2.h>
	//#include <dlaexc.h>
	//#include <dlag2.h>
	//#include <dlags2.h>
	//#include <dlagtf.h>
	//#include <dlagtm.h>
	//#include <dlagts.h>
	//#include <dlagv2.h>
	//#include <dlahqr.h>
	//#include <dlahrd.h>
	//#include <dlaic1.h>
	//#include <dlaln2.h>
	//#include <dlals0.h>
	//#include <dlalsa.h>
	//#include <dlalsd.h>
	//#include <dlamch.h>
	//#include <dlamrg.h>
	//#include <dlangb.h>
	//#include <dlange.h>
	//#include <dlangt.h>
	//#include <dlanhs.h>
	//#include <dlansb.h>
	//#include <dlansp.h>
	//#include <dlanst.h>
	//#include <dlansy.h>
	//#include <dlantb.h>
	//#include <dlantp.h>
	//#include <dlantr.h>
	//#include <dlanv2.h>
	//#include <dlapll.h>
	//#include <dlapmt.h>
	void lapmt_(logical *forwrd, integer *m, integer *n, complexRopt *x, integer *ldx, integer *k)
	{
#ifndef MATLAB_MEX_FILE
		clapmt_(forwrd, m, n, x, ldx, k);
#else
		clapmt_((integer *)forwrd, m, n, (realRopt *)x, ldx, k);
#endif
	};
	void lapmt_(logical *forwrd, integer *m, integer *n, doublerealRopt *x, integer *ldx, integer *k)
	{
#ifndef MATLAB_MEX_FILE
		dlapmt_(forwrd, m, n, x, ldx, k);
#else
		dlapmt_((integer *)forwrd, m, n, x, ldx, k);
#endif
	};
	void lapmt_(logical *forwrd, integer *m, integer *n, realRopt *x, integer *ldx, integer *k)
	{
#ifndef MATLAB_MEX_FILE
		slapmt_(forwrd, m, n, x, ldx, k);
#else
		slapmt_((integer *)forwrd, m, n, x, ldx, k);
#endif
	};
	void lapmt_(logical *forwrd, integer *m, integer *n, doublecomplexRopt *x, integer *ldx, integer *k)
	{
#ifndef MATLAB_MEX_FILE
		zlapmt_(forwrd, m, n, x, ldx, k);
#else
		zlapmt_((integer *)forwrd, m, n, (doublerealRopt *)x, ldx, k);
#endif
	};

	//#include <dlapy2.h>
	//#include <dlapy3.h>
	//#include <dlaqgb.h>
	//#include <dlaqge.h>
	//#include <dlaqp2.h>
	//#include <dlaqps.h>
	//#include <dlaqsb.h>
	//#include <dlaqsp.h>
	//#include <dlaqsy.h>
	//#include <dlaqtr.h>
	//#include <dlar1v.h>
	//#include <dlar2v.h>
	//#include <dlarf.h>
	//#include <dlarfb.h>
	//#include <dlarfg.h>
	//#include <dlarft.h>
	//#include <dlarfx.h>
	void larfx_(char *side, integer *m, integer *n, complexRopt *v, complexRopt *tau, complexRopt *c__, integer *ldc, complexRopt *work)
	{
#ifndef MATLAB_MEX_FILE
		clarfx_(side, m, n, v, tau, c__, ldc, work);
#else
		clarfx_(side, m, n, (realRopt *)v, (realRopt *)tau, (realRopt *)c__, ldc, (realRopt *)work);
#endif
	};
	void larfx_(char *side, integer *m, integer *n, doublerealRopt *v, doublerealRopt *tau, doublerealRopt *c__, integer *ldc, doublerealRopt *work)
	{
		dlarfx_(side, m, n, v, tau, c__, ldc, work);
	};
	void larfx_(char *side, integer *m, integer *n, realRopt *v, realRopt *tau, realRopt *c__, integer *ldc, realRopt *work)
	{
		slarfx_(side, m, n, v, tau, c__, ldc, work);
	};
	void larfx_(char *side, integer *m, integer *n, doublecomplexRopt *v, doublecomplexRopt *tau, doublecomplexRopt *c__, integer *ldc, doublecomplexRopt *work)
	{
#ifndef MATLAB_MEX_FILE
		zlarfx_(side, m, n, v, tau, c__, ldc, work);
#else
		zlarfx_(side, m, n, (doublerealRopt *)v, (doublerealRopt *)tau, (doublerealRopt *)c__, ldc, (doublerealRopt *)work);
#endif
	};
	//#include <dlargv.h>
	//#include <dlarnv.h>
	//#include <dlarrb.h>
	//#include <dlarre.h>
	//#include <dlarrf.h>
	//#include <dlarrv.h>
	//#include <dlartg.h>
	//#include <dlartv.h>
	//#include <dlaruv.h>
	//#include <dlarz.h>
	//#include <dlarzb.h>
	//#include <dlarzt.h>
	//#include <dlas2.h>
	//#include <dlascl.h>
	//#include <dlasd0.h>
	//#include <dlasd1.h>
	//#include <dlasd2.h>
	//#include <dlasd3.h>
	//#include <dlasd4.h>
	//#include <dlasd5.h>
	//#include <dlasd6.h>
	//#include <dlasd7.h>
	//#include <dlasd8.h>
	//#include <dlasd9.h>
	//#include <dlasda.h>
	//#include <dlasdq.h>
	//#include <dlasdt.h>
	//#include <dlaset.h>
	//#include <dlasq1.h>
	//#include <dlasq2.h>
	//#include <dlasq3.h>
	//#include <dlasq4.h>
	//#include <dlasq5.h>
	//#include <dlasq6.h>
	//#include <dlasr.h>
	//#include <dlasrt.h>
	//#include <dlassq.h>
	//#include <dlasv2.h>
	//#include <dlaswp.h>
	//#include <dlasy2.h>
	//#include <dlasyf.h>
	//#include <dlatbs.h>
	//#include <dlatdf.h>
	//#include <dlatps.h>
	//#include <dlatrd.h>
	//#include <dlatrs.h>
	//#include <dlatrz.h>
	//#include <dlatzm.h>
	//#include <dlauu2.h>
	//#include <dlauum.h>
	//#include <dopgtr.h>
	//#include <dopmtr.h>
	//#include <dorg2l.h>
	//#include <dorg2r.h>
	//#include <dorgbr.h>
	//#include <dorghr.h>
	//#include <dorgl2.h>
	//#include <dorglq.h>
	//#include <dorgql.h>
	//#include <dorgqr.h>
	void orgqr_(integer *m, integer *n, integer *k, doublerealRopt *a, integer *lda, doublerealRopt *tau, doublerealRopt *work, integer *lwork, integer *info)
	{
		dorgqr_(m, n, k, a, lda, tau, work, lwork, info);
	};
	void orgqr_(integer *m, integer *n, integer *k, realRopt *a, integer *lda, realRopt *tau, realRopt *work, integer *lwork, integer *info)
	{
		sorgqr_(m, n, k, a, lda, tau, work, lwork, info);
	};

	//#include <dorgr2.h>
	//#include <dorgrq.h>
	//#include <dorgtr.h>
	//#include <dorm2l.h>
	//#include <dorm2r.h>
	//#include <dormbr.h>
	//#include <dormhr.h>
	//#include <dorml2.h>
	//#include <dormlq.h>
	//#include <dormql.h>
	//#include <dormqr.h>
	void ormqr_(char *side, char *trans, integer *m, integer *n, integer *k, doublerealRopt *a, integer *lda, doublerealRopt *tau, doublerealRopt *c__, integer *ldc, doublerealRopt *work, integer *lwork, integer *info)
	{
		dormqr_(side, trans, m, n, k, a, lda, tau, c__, ldc, work, lwork, info);
	};
	void ormqr_(char *side, char *trans, integer *m, integer *n, integer *k, realRopt *a, integer *lda, realRopt *tau, realRopt *c__, integer *ldc, realRopt *work, integer *lwork, integer *info)
	{
		sormqr_(side, trans, m, n, k, a, lda, tau, c__, ldc, work, lwork, info);
	};

	//#include <dormr2.h>
	//#include <dormr3.h>
	//#include <dormrq.h>
	//#include <dormrz.h>
	//#include <dormtr.h>
	//#include <dpbcon.h>
	//#include <dpbequ.h>
	//#include <dpbrfs.h>
	//#include <dpbstf.h>
	//#include <dpbsv.h>
	//#include <dpbsvx.h>
	//#include <dpbtf2.h>
	//#include <dpbtrf.h>
	//#include <dpbtrs.h>
	//#include <dpocon.h>
	//#include <dpoequ.h>
	//#include <dporfs.h>
	//#include <dposv.h>
	//#include <dposvx.h>
	//#include <dpotf2.h>
	//#include <dpotrf.h>
	void potrf_(char *uplo, integer *n, complexRopt *a, integer *lda, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cpotrf_(uplo, n, a, lda, info);
#else
		cpotrf_(uplo, n, (float *)a, lda, info);
#endif
	};
	void potrf_(char *uplo, integer *n, doublerealRopt *a, integer *lda, integer *info)
	{
		dpotrf_(uplo, n, a, lda, info);
	};
	void potrf_(char *uplo, integer *n, realRopt *a, integer *lda, integer *info)
	{
		spotrf_(uplo, n, a, lda, info);
	};
	void potrf_(char *uplo, integer *n, doublecomplexRopt *a, integer *lda, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zpotrf_(uplo, n, a, lda, info);
#else
		zpotrf_(uplo, n, (double *)a, lda, info);
#endif
	};

	//#include <dpotri.h>
	//#include <dpotrs.h>
	void potrs_(char *uplo, integer *n, integer *nrhs, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cpotrs_(uplo, n, nrhs, a, lda, b, ldb, info);
#else
		cpotrs_(uplo, n, nrhs, (float *)a, lda, (float *)b, ldb, info);
#endif
	};
	void potrs_(char *uplo, integer *n, integer *nrhs, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, integer *info)
	{
		dpotrs_(uplo, n, nrhs, a, lda, b, ldb, info);
	};
	void potrs_(char *uplo, integer *n, integer *nrhs, realRopt *a, integer *lda, realRopt *b, integer *ldb, integer *info)
	{
		spotrs_(uplo, n, nrhs, a, lda, b, ldb, info);
	};
	void potrs_(char *uplo, integer *n, integer *nrhs, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zpotrs_(uplo, n, nrhs, a, lda, b, ldb, info);
#else
		zpotrs_(uplo, n, nrhs, (double *)a, lda, (double *)b, ldb, info);
#endif
	};

	//#include <dppcon.h>
	//#include <dppequ.h>
	//#include <dpprfs.h>
	//#include <dppsv.h>
	//#include <dppsvx.h>
	//#include <dpptrf.h>
	//#include <dpptri.h>
	//#include <dpptrs.h>
	//#include <dptcon.h>
	//#include <dpteqr.h>
	//#include <dptrfs.h>
	//#include <dptsv.h>
	//#include <dptsvx.h>
	//#include <dpttrf.h>
	//#include <dpttrs.h>
	//#include <dptts2.h>
	//#include <drscl.h>
	//#include <dsbev.h>
	//#include <dsbevd.h>
	//#include <dsbevx.h>
	//#include <dsbgst.h>
	//#include <dsbgv.h>
	//#include <dsbgvd.h>
	//#include <dsbgvx.h>
	//#include <dsbtrd.h>
	//#include <dsecnd.h>
	//#include <dspcon.h>
	//#include <dspev.h>
	//#include <dspevd.h>
	//#include <dspevx.h>
	//#include <dspgst.h>
	//#include <dspgv.h>
	//#include <dspgvd.h>
	//#include <dspgvx.h>
	//#include <dsprfs.h>
	//#include <dspsv.h>
	//#include <dspsvx.h>
	//#include <dsptrd.h>
	//#include <dsptrf.h>
	//#include <dsptri.h>
	//#include <dsptrs.h>
	//#include <dstebz.h>
	//#include <dstedc.h>
	//#include <dstegr.h>
	//#include <dstein.h>
	//#include <dsteqr.h>
	//#include <dsterf.h>
	//#include <dstev.h>
	//#include <dstevd.h>
	//#include <dstevr.h>
	//#include <dstevx.h>
	//#include <dsycon.h>
	//#include <dsyev.h>
	void syev_(char *jobz, char *uplo, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *w, doublerealRopt *work, integer *lwork, integer *info)
	{
		dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
	};
	void syev_(char *jobz, char *uplo, integer *n, realRopt *a, integer *lda, realRopt *w, realRopt *work, integer *lwork, integer *info)
	{
		ssyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
	};

	//#include <dsyevd.h>
	void syevd_(char *jobz, char *uplo, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *w, doublerealRopt *work, integer *lwork, integer *iwork, integer *liwork, integer *info)
	{
		dsyevd_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
	};
	void syevd_(char *jobz, char *uplo, integer *n, realRopt *a, integer *lda, realRopt *w, realRopt *work, integer *lwork, integer *iwork, integer *liwork, integer *info)
	{
		ssyevd_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
	};
	//#include <dsyevr.h>
	void syevr_(char *jobz, char *range, char *uplo, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *vl, doublerealRopt *vu, integer *il, integer *iu, doublerealRopt *abstol, integer *m, doublerealRopt *w, doublerealRopt *z__, integer *ldz, integer *isuppz, doublerealRopt *work, integer *lwork, integer *iwork, integer *liwork, integer *info)
	{
		dsyevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z__, ldz, isuppz, work, lwork, iwork, liwork, info);
	};
	void syevr_(char *jobz, char *range, char *uplo, integer *n, realRopt *a, integer *lda, realRopt *vl, realRopt *vu, integer *il, integer *iu, realRopt *abstol, integer *m, realRopt *w, realRopt *z__, integer *ldz, integer *isuppz, realRopt *work, integer *lwork, integer *iwork, integer *liwork, integer *info)
	{
		ssyevr_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z__, ldz, isuppz, work, lwork, iwork, liwork, info);
	};
	//#include <dsyevx.h>
	void syevx_(char *jobz, char *range, char *uplo, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *vl, doublerealRopt *vu, integer *il, integer *iu, doublerealRopt *abstol, integer *m, doublerealRopt *w, doublerealRopt *z__, integer *ldz, doublerealRopt *work, integer *lwork, integer *iwork, integer *ifail, integer *info)
	{
		dsyevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z__, ldz, work, lwork, iwork, ifail, info);
	};
	void syevx_(char *jobz, char *range, char *uplo, integer *n, realRopt *a, integer *lda, realRopt *vl, realRopt *vu, integer *il, integer *iu, realRopt *abstol, integer *m, realRopt *w, realRopt *z__, integer *ldz, realRopt *work, integer *lwork, integer *iwork, integer *ifail, integer *info)
	{
		ssyevx_(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z__, ldz, work, lwork, iwork, ifail, info);
	};
	//#include <dsygs2.h>
	//#include <dsygst.h>
	//#include <dsygv.h>
	//#include <dsygvd.h>
	//#include <dsygvx.h>
	//#include <dsyrfs.h>
	//#include <dsysv.h>
	//#include <dsysvx.h>
	//#include <dsytd2.h>
	//#include <dsytf2.h>
	//#include <dsytrd.h>
	//#include <dsytrf.h>
	//#include <dsytri.h>
	//#include <dsytrs.h>
	//#include <dtbcon.h>
	//#include <dtbrfs.h>
	//#include <dtbtrs.h>
	//#include <dtgevc.h>
	//#include <dtgex2.h>
	//#include <dtgexc.h>
	//#include <dtgsen.h>
	//#include <dtgsja.h>
	//#include <dtgsna.h>
	//#include <dtgsy2.h>
	//#include <dtgsyl.h>
	void tgsyl_(char *trans, integer *ijob, integer *m, integer *n, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, complexRopt *c__, integer *ldc, complexRopt *d__, integer *ldd, complexRopt *e, integer *lde, complexRopt *f, integer *ldf, realRopt *scale, realRopt *dif, complexRopt *work, integer *lwork, integer *iwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		ctgsyl_(trans, ijob, m, n, a, lda, b, ldb, c__, ldc, d__, ldd, e, lde, f, ldf, scale, dif, work, lwork, iwork, info);
#else
		ctgsyl_(trans, ijob, m, n, (float *)a, lda, (float *)b, ldb, (float *)c__, ldc, (float *)d__, ldd, (float *)e, lde, (float *)f, ldf, scale, dif, (float *)work, lwork, iwork, info);
#endif
	};
	void tgsyl_(char *trans, integer *ijob, integer *m, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, doublerealRopt *c__, integer *ldc, doublerealRopt *d__, integer *ldd, doublerealRopt *e, integer *lde, doublerealRopt *f, integer *ldf, doublerealRopt *scale, doublerealRopt *dif, doublerealRopt *work, integer *lwork, integer *iwork, integer *info)
	{
		dtgsyl_(trans, ijob, m, n, a, lda, b, ldb, c__, ldc, d__, ldd, e, lde, f, ldf, scale, dif, work, lwork, iwork, info);
	};
	void tgsyl_(char *trans, integer *ijob, integer *m, integer *n, realRopt *a, integer *lda, realRopt *b, integer *ldb, realRopt *c__, integer *ldc, realRopt *d__, integer *ldd, realRopt *e, integer *lde, realRopt *f, integer *ldf, realRopt *scale, realRopt *dif, realRopt *work, integer *lwork, integer *iwork, integer *info)
	{
		stgsyl_(trans, ijob, m, n, a, lda, b, ldb, c__, ldc, d__, ldd, e, lde, f, ldf, scale, dif, work, lwork, iwork, info);
	};
	void tgsyl_(char *trans, integer *ijob, integer *m, integer *n, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *c__, integer *ldc, doublecomplexRopt *d__, integer *ldd, doublecomplexRopt *e, integer *lde, doublecomplexRopt *f, integer *ldf, doublerealRopt *scale, doublerealRopt *dif, doublecomplexRopt *work, integer *lwork, integer *iwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		ztgsyl_(trans, ijob, m, n, a, lda, b, ldb, c__, ldc, d__, ldd, e, lde, f, ldf, scale, dif, work, lwork, iwork, info);
#else
		ztgsyl_(trans, ijob, m, n, (double *)a, lda, (double *)b, ldb, (double *)c__, ldc, (double *)d__, ldd, (double *)e, lde, (double *)f, ldf, scale, dif, (double *)work, lwork, iwork, info);
#endif
	};
	//#include <dtpcon.h>
	//#include <dtprfs.h>
	//#include <dtptri.h>
	//#include <dtptrs.h>
	//#include <dtrcon.h>
	//#include <dtrevc.h>
	//#include <dtrexc.h>
	//#include <dtrrfs.h>
	//#include <dtrsen.h>
	//#include <dtrsna.h>
	//#include <dtrsyl.h>
	//#include <dtrti2.h>
	//#include <dtrtri.h>
	//#include <dtrtrs.h>
	void trtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		ctrtrs_(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
#else
		ctrtrs_(uplo, trans, diag, n, nrhs, (float *)a, lda, (float *)b, ldb, info);
#endif
	};
	void trtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, integer *info)
	{
		dtrtrs_(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
	};
	void trtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, realRopt *a, integer *lda, realRopt *b, integer *ldb, integer *info)
	{
		strtrs_(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
	};
	void trtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		ztrtrs_(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info);
#else
		ztrtrs_(uplo, trans, diag, n, nrhs, (double *)a, lda, (double *)b, ldb, info);
#endif
	};
	//#include <dtzrqf.h>
	//#include <dtzrzf.h>
	//#include <dzsum1.h>

		//zgegs_
	void gegs_(char *jobvsl, char *jobvsr, integer *n, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, complexRopt *alpha, complexRopt *beta, complexRopt *vsl, integer *ldvsl, complexRopt *vsr, integer *ldvsr, complexRopt *work, integer *lwork, realRopt *rwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cgegs_(jobvsl, jobvsr, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, info);
#else
		cgegs_(jobvsl, jobvsr, n, (realRopt *)a, lda, (realRopt *)b, ldb, (realRopt *)alpha, (realRopt *)beta, (realRopt *)vsl, ldvsl, (realRopt *)vsr, ldvsr, (realRopt *)work, lwork, rwork, info);
#endif
	};
	void gegs_(char *jobvsl, char *jobvsr, integer *n, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *alpha, doublecomplexRopt *beta, doublecomplexRopt *vsl, integer *ldvsl, doublecomplexRopt *vsr, integer *ldvsr, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zgegs_(jobvsl, jobvsr, n, a, lda, b, ldb, alpha, beta, vsl, ldvsl, vsr, ldvsr, work, lwork, rwork, info);
#else
		zgegs_(jobvsl, jobvsr, n, (doublerealRopt *)a, lda, (doublerealRopt *)b, ldb, (doublerealRopt *)alpha, (doublerealRopt *)beta, (doublerealRopt *)vsl, ldvsl, (doublerealRopt *)vsr, ldvsr, (doublerealRopt *)work, lwork, rwork, info);
#endif
	};

	//zunmqr_
	void unmqr_(char *side, char *trans, integer *m, integer *n, integer *k, complexRopt *a, integer *lda, complexRopt *tau, complexRopt *c__, integer *ldc, complexRopt *work, integer *lwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cunmqr_(side, trans, m, n, k, a, lda, tau, c__, ldc, work, lwork, info);
#else
		cunmqr_(side, trans, m, n, k, (float *)a, lda, (float *)tau, (float *)c__, ldc, (float *)work, lwork, info);
#endif
	};
	void unmqr_(char *side, char *trans, integer *m, integer *n, integer *k, doublecomplexRopt *a, integer *lda, doublecomplexRopt *tau, doublecomplexRopt *c__, integer *ldc, doublecomplexRopt *work, integer *lwork, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zunmqr_(side, trans, m, n, k, a, lda, tau, c__, ldc, work, lwork, info);
#else
		zunmqr_(side, trans, m, n, k, (double *)a, lda, (double *)tau, (double *)c__, ldc, (double *)work, lwork, info);
#endif
	};

	//zpotri_
	void potri_(char *uplo, integer *n, complexRopt *a, integer *lda, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		cpotri_(uplo, n, a, lda, info);
#else
		cpotri_(uplo, n, (realRopt *)a, lda, info);
#endif
	};
	void potri_(char *uplo, integer *n, doublecomplexRopt *a, integer *lda, integer *info)
	{
#ifndef MATLAB_MEX_FILE
		zpotri_(uplo, n, a, lda, info);
#else
		zpotri_(uplo, n, (doublerealRopt *)a, lda, info);
#endif
	};


    //zheev_
    void heev_(char *jobz, char *uplo, integer *n, complexRopt *A, integer *LDA, realRopt *W, complexRopt *work, integer *lwork, realRopt *rwork, integer *info)
    {
    #ifndef MATLAB_MEX_FILE
        cheev_(jobz, uplo, n, A, LDA, W, work, lwork, rwork, info);
    #else
        cheev_(jobz, uplo, n, (realRopt *) A, LDA, W, (realRopt *) work, lwork, rwork, info);
    #endif
    };

    void heev_(char *jobz, char *uplo, integer *n, doublecomplexRopt *A, integer *LDA, doublerealRopt *W, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *info)
    {
    #ifndef MATLAB_MEX_FILE
        zheev_(jobz, uplo, n, A, LDA, W, work, lwork, rwork, info);
    #else
        zheev_(jobz, uplo, n, (doublerealRopt *) A, LDA, W, (doublerealRopt *) work, lwork, rwork, info);
    #endif
    };

    // zungqr_
    void ungqr_(integer *m, integer *n, integer *k, complexRopt *a, integer *lda, complexRopt *tau, complexRopt *work, integer *lwork, integer *info)
    {
#ifndef MATLAB_MEX_FILE
        cungqr_(m, n, k, a, lda, tau, work, lwork, info);
#else
        cungqr_(m, n, k, (realRopt *) a, lda, (realRopt *) tau, (realRopt *) work, lwork, info);
#endif
    };

    void ungqr_(integer *m, integer *n, integer *k, doublecomplexRopt *a, integer *lda, doublecomplexRopt *tau, doublecomplexRopt *work, integer *lwork, integer *info)
    {
#ifndef MATLAB_MEX_FILE
        zungqr_(m, n, k, a, lda, tau, work, lwork, info);
#else
        zungqr_(m, n, k, (doublerealRopt *) a, lda, (doublerealRopt *) tau, (doublerealRopt *) work, lwork, info);
#endif
    };

	//BLAS_dusmm
	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa,
		int nrhs, float alpha, blas_sparse_matrix A, const float *b, int ldb,
		float *c, int ldc)
	{
		BLAS_susmm(order, transa, nrhs, alpha, A, b, ldb, c, ldc);
	};
	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa,
		int nrhs, double alpha, blas_sparse_matrix A, const double *b,
		int ldb, double *c, int ldc)
	{
		BLAS_dusmm(order, transa, nrhs, alpha, A, b, ldb, c, ldc);
	};
	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa,
		int nrhs, const complexRopt *alpha, blas_sparse_matrix A, const complexRopt *b,
		int ldb, complexRopt *c, int ldc)
	{
		BLAS_cusmm(order, transa, nrhs, alpha, A, b, ldb, c, ldc);
	};
	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa,
		int nrhs, const doublecomplexRopt *alpha, blas_sparse_matrix A, const doublecomplexRopt *b,
		int ldb, doublecomplexRopt *c, int ldc)
	{
		BLAS_zusmm(order, transa, nrhs, alpha, A, b, ldb, c, ldc);
	};
    
//#if ! defined(__APPLE__) || defined(MATLAB_MEX_FILE)
//	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, integer nrhs, float alpha, blas_sparse_matrix A, const float *b, integer ldb, float *c, integer ldc)
//	{
//		BLAS_susmm(order, transa, nrhs, alpha, A, b, ldb, c, ldc);
//	};
//	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, integer nrhs, double alpha, blas_sparse_matrix A, const double *b, integer ldb, double *c, integer ldc)
//	{
//		BLAS_dusmm(order, transa, nrhs, alpha, A, b, ldb, c, ldc);
//	};
//	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, integer nrhs, const complexRopt *alpha, blas_sparse_matrix A, const complexRopt *b, integer ldb, complexRopt *c, integer ldc)
//	{
//		BLAS_cusmm(order, transa, nrhs, alpha, A, b, ldb, c, ldc);
//	};
//	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, integer nrhs, const doublecomplexRopt *alpha, blas_sparse_matrix A, const doublecomplexRopt *b, integer ldb, doublecomplexRopt *c, integer ldc)
//	{
//		BLAS_zusmm(order, transa, nrhs, alpha, A, b, ldb, c, ldc);
//	};
//#endif

	void BLAS_uscr_insert_entries(blas_sparse_matrix A, integer nz, const float *val, const integer *indx, const integer *jndx)
	{
		int *indxptr = new int[nz * 2];
		int *jndxptr = indxptr + nz;
		for (integer i = 0; i < nz; i++)
		{
			indxptr[i] = indx[i];
			jndxptr[i] = jndx[i];
		}
		BLAS_suscr_insert_entries(A, nz, val, indxptr, jndxptr);
		delete[] indxptr;
	}
	void BLAS_uscr_insert_entries(blas_sparse_matrix A, integer nz, const double *val, const integer *indx, const integer *jndx)
	{
		int *indxptr = new int[nz * 2];
		int *jndxptr = indxptr + nz;
		for (integer i = 0; i < nz; i++)
		{
			indxptr[i] = indx[i];
			jndxptr[i] = jndx[i];
		}
		BLAS_duscr_insert_entries(A, nz, val, indxptr, jndxptr);
		delete[] indxptr;
	}
	void BLAS_uscr_insert_entries(blas_sparse_matrix A, integer nz, const complexRopt *val, const integer *indx, const integer *jndx)
	{
		int *indxptr = new int[nz * 2];
		int *jndxptr = indxptr + nz;
		for (integer i = 0; i < nz; i++)
		{
			indxptr[i] = indx[i];
			jndxptr[i] = jndx[i];
		}
		BLAS_cuscr_insert_entries(A, nz, val, indxptr, jndxptr);
		delete[] indxptr;
	}
	void BLAS_uscr_insert_entries(blas_sparse_matrix A, integer nz, const doublecomplexRopt *val, const integer *indx, const integer *jndx)
	{
		int *indxptr = new int[nz * 2];
		int *jndxptr = indxptr + nz;
		for (integer i = 0; i < nz; i++)
		{
			indxptr[i] = indx[i];
			jndxptr[i] = jndx[i];
		}
		BLAS_zuscr_insert_entries(A, nz, val, indxptr, jndxptr);
		delete[] indxptr;
	}

//#if ! defined(__APPLE__) || defined(MATLAB_MEX_FILE)
//	void BLAS_uscr_insert_entries(blas_sparse_matrix A, int nz, const float *val, const int *indx, const int *jndx)
//	{
//		BLAS_suscr_insert_entries(A, nz, val, indx, jndx);
//	}
//	void BLAS_uscr_insert_entries(blas_sparse_matrix A, int nz, const double *val, const int *indx, const int *jndx)
//	{
//		BLAS_duscr_insert_entries(A, nz, val, indx, jndx);
//	}
//	void BLAS_uscr_insert_entries(blas_sparse_matrix A, int nz, const complexRopt *val, const int *indx, const int *jndx)
//	{
//		BLAS_cuscr_insert_entries(A, nz, val, indx, jndx);
//	}
//	void BLAS_uscr_insert_entries(blas_sparse_matrix A, int nz, const doublecomplexRopt *val, const int *indx, const int *jndx)
//	{
//		BLAS_zuscr_insert_entries(A, nz, val, indx, jndx);
//	}
//#endif
}; /*end of ROPTLIB namespace*/
