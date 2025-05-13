/*
Some blas and lapack functions are not easy to use. This class defines wrapper functions
which can be used to solve some problems or decomposition, e.g., SVD decomposition, Sylevster equation.
More functions will be added in this class.

-----WH
*/
#pragma once
#ifndef BLASLAPACKCPPWRAPPER_H
#define BLASLAPACKCPPWRAPPER_H

#include "others_SparseBLAS_blas_sparse.h"

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "blas.h"
#include "lapack.h"
#include "matrix.h"
#define logical integer
#define integer ptrdiff_t
#define realRopt float
#define doublerealRopt double
typedef integer (*L_fp)();
#endif

#ifndef MATLAB_MEX_FILE
/* blas and lapack related
BLAS */
#include <thirdparty/roptlib/cwrapper/blas/caxpy.h>
#include <thirdparty/roptlib/cwrapper/blas/ccopy.h>
#include <thirdparty/roptlib/cwrapper/blas/cdotc.h>
#include <thirdparty/roptlib/cwrapper/blas/cdotu.h>
#include <thirdparty/roptlib/cwrapper/blas/cgbmv.h>
#include <thirdparty/roptlib/cwrapper/blas/cgemm.h>
#include <thirdparty/roptlib/cwrapper/blas/cgemv.h>
#include <thirdparty/roptlib/cwrapper/blas/cgerc.h>
#include <thirdparty/roptlib/cwrapper/blas/cgeru.h>
#include <thirdparty/roptlib/cwrapper/blas/chbmv.h>
#include <thirdparty/roptlib/cwrapper/blas/chemm.h>
#include <thirdparty/roptlib/cwrapper/blas/chemv.h>
#include <thirdparty/roptlib/cwrapper/blas/cher.h>
#include <thirdparty/roptlib/cwrapper/blas/cher2.h>
#include <thirdparty/roptlib/cwrapper/blas/cher2k.h>
#include <thirdparty/roptlib/cwrapper/blas/cherk.h>
#include <thirdparty/roptlib/cwrapper/blas/chpmv.h>
#include <thirdparty/roptlib/cwrapper/blas/chpr.h>
#include <thirdparty/roptlib/cwrapper/blas/chpr2.h>
#include <thirdparty/roptlib/cwrapper/blas/crotg.h>
#include <thirdparty/roptlib/cwrapper/blas/cscal.h>
#include <thirdparty/roptlib/cwrapper/blas/csscal.h>
#include <thirdparty/roptlib/cwrapper/blas/cswap.h>
#include <thirdparty/roptlib/cwrapper/blas/csymm.h>
#include <thirdparty/roptlib/cwrapper/blas/csyr2k.h>
#include <thirdparty/roptlib/cwrapper/blas/csyrk.h>
#include <thirdparty/roptlib/cwrapper/blas/ctbmv.h>
#include <thirdparty/roptlib/cwrapper/blas/ctbsv.h>
#include <thirdparty/roptlib/cwrapper/blas/ctpmv.h>
#include <thirdparty/roptlib/cwrapper/blas/ctpsv.h>
#include <thirdparty/roptlib/cwrapper/blas/ctrmm.h>
#include <thirdparty/roptlib/cwrapper/blas/ctrmv.h>
#include <thirdparty/roptlib/cwrapper/blas/ctrsm.h>
#include <thirdparty/roptlib/cwrapper/blas/ctrsv.h>
#include <thirdparty/roptlib/cwrapper/blas/dasum.h>
#include <thirdparty/roptlib/cwrapper/blas/daxpy.h>
#include <thirdparty/roptlib/cwrapper/blas/dcabs1.h>
#include <thirdparty/roptlib/cwrapper/blas/dcopy.h>
#include <thirdparty/roptlib/cwrapper/blas/ddot.h>
#include <thirdparty/roptlib/cwrapper/blas/dgbmv.h>
#include <thirdparty/roptlib/cwrapper/blas/dgemm.h>
#include <thirdparty/roptlib/cwrapper/blas/dgemv.h>
#include <thirdparty/roptlib/cwrapper/blas/dger.h>
#include <thirdparty/roptlib/cwrapper/blas/dnrm2.h>
#include <thirdparty/roptlib/cwrapper/blas/drot.h>
#include <thirdparty/roptlib/cwrapper/blas/drotg.h>
#include <thirdparty/roptlib/cwrapper/blas/dsbmv.h>
#include <thirdparty/roptlib/cwrapper/blas/dscal.h>
#include <thirdparty/roptlib/cwrapper/blas/dspmv.h>
#include <thirdparty/roptlib/cwrapper/blas/dspr.h>
#include <thirdparty/roptlib/cwrapper/blas/dspr2.h>
#include <thirdparty/roptlib/cwrapper/blas/dswap.h>
#include <thirdparty/roptlib/cwrapper/blas/dsymm.h>
#include <thirdparty/roptlib/cwrapper/blas/dsymv.h>
#include <thirdparty/roptlib/cwrapper/blas/dsyr.h>
#include <thirdparty/roptlib/cwrapper/blas/dsyr2.h>
#include <thirdparty/roptlib/cwrapper/blas/dsyr2k.h>
#include <thirdparty/roptlib/cwrapper/blas/dsyrk.h>
#include <thirdparty/roptlib/cwrapper/blas/dtbmv.h>
#include <thirdparty/roptlib/cwrapper/blas/dtbsv.h>
#include <thirdparty/roptlib/cwrapper/blas/dtpmv.h>
#include <thirdparty/roptlib/cwrapper/blas/dtpsv.h>
#include <thirdparty/roptlib/cwrapper/blas/dtrmm.h>
#include <thirdparty/roptlib/cwrapper/blas/dtrmv.h>
#include <thirdparty/roptlib/cwrapper/blas/dtrsm.h>
#include <thirdparty/roptlib/cwrapper/blas/dtrsv.h>
#include <thirdparty/roptlib/cwrapper/blas/dzasum.h>
#include <thirdparty/roptlib/cwrapper/blas/dznrm2.h>
#include <thirdparty/roptlib/cwrapper/blas/f2c.h>
#include <thirdparty/roptlib/cwrapper/blas/icamax.h>
#include <thirdparty/roptlib/cwrapper/blas/idamax.h>
#include <thirdparty/roptlib/cwrapper/blas/isamax.h>
#include <thirdparty/roptlib/cwrapper/blas/izamax.h>
#include <thirdparty/roptlib/cwrapper/blas/lsame.h>
#include <thirdparty/roptlib/cwrapper/blas/sasum.h>
#include <thirdparty/roptlib/cwrapper/blas/saxpy.h>
#include <thirdparty/roptlib/cwrapper/blas/scasum.h>
#include <thirdparty/roptlib/cwrapper/blas/scnrm2.h>
#include <thirdparty/roptlib/cwrapper/blas/scopy.h>
#include <thirdparty/roptlib/cwrapper/blas/sdot.h>
#include <thirdparty/roptlib/cwrapper/blas/sgbmv.h>
#include <thirdparty/roptlib/cwrapper/blas/sgemm.h>
#include <thirdparty/roptlib/cwrapper/blas/sgemv.h>
#include <thirdparty/roptlib/cwrapper/blas/sger.h>
#include <thirdparty/roptlib/cwrapper/blas/snrm2.h>
#include <thirdparty/roptlib/cwrapper/blas/srot.h>
#include <thirdparty/roptlib/cwrapper/blas/srotg.h>
#include <thirdparty/roptlib/cwrapper/blas/ssbmv.h>
#include <thirdparty/roptlib/cwrapper/blas/sscal.h>
#include <thirdparty/roptlib/cwrapper/blas/sspmv.h>
#include <thirdparty/roptlib/cwrapper/blas/sspr.h>
#include <thirdparty/roptlib/cwrapper/blas/sspr2.h>
#include <thirdparty/roptlib/cwrapper/blas/sswap.h>
#include <thirdparty/roptlib/cwrapper/blas/ssymm.h>
#include <thirdparty/roptlib/cwrapper/blas/ssymv.h>
#include <thirdparty/roptlib/cwrapper/blas/ssyr.h>
#include <thirdparty/roptlib/cwrapper/blas/ssyr2.h>
#include <thirdparty/roptlib/cwrapper/blas/ssyr2k.h>
#include <thirdparty/roptlib/cwrapper/blas/ssyrk.h>
#include <thirdparty/roptlib/cwrapper/blas/stbmv.h>
#include <thirdparty/roptlib/cwrapper/blas/stbsv.h>
#include <thirdparty/roptlib/cwrapper/blas/stpmv.h>
#include <thirdparty/roptlib/cwrapper/blas/stpsv.h>
#include <thirdparty/roptlib/cwrapper/blas/strmm.h>
#include <thirdparty/roptlib/cwrapper/blas/strmv.h>
#include <thirdparty/roptlib/cwrapper/blas/strsm.h>
#include <thirdparty/roptlib/cwrapper/blas/strsv.h>
#include <thirdparty/roptlib/cwrapper/blas/xerbla.h>
#include <thirdparty/roptlib/cwrapper/blas/zaxpy.h>
#include <thirdparty/roptlib/cwrapper/blas/zcopy.h>
#include <thirdparty/roptlib/cwrapper/blas/zdotc.h>
#include <thirdparty/roptlib/cwrapper/blas/zdotu.h>
#include <thirdparty/roptlib/cwrapper/blas/zdscal.h>
#include <thirdparty/roptlib/cwrapper/blas/zgbmv.h>
#include <thirdparty/roptlib/cwrapper/blas/zgemm.h>
#include <thirdparty/roptlib/cwrapper/blas/zgemv.h>
#include <thirdparty/roptlib/cwrapper/blas/zgerc.h>
#include <thirdparty/roptlib/cwrapper/blas/zgeru.h>
#include <thirdparty/roptlib/cwrapper/blas/zhbmv.h>
#include <thirdparty/roptlib/cwrapper/blas/zhemm.h>
#include <thirdparty/roptlib/cwrapper/blas/zhemv.h>
#include <thirdparty/roptlib/cwrapper/blas/zher.h>
#include <thirdparty/roptlib/cwrapper/blas/zher2.h>
#include <thirdparty/roptlib/cwrapper/blas/zher2k.h>
#include <thirdparty/roptlib/cwrapper/blas/zherk.h>
#include <thirdparty/roptlib/cwrapper/blas/zhpmv.h>
#include <thirdparty/roptlib/cwrapper/blas/zhpr.h>
#include <thirdparty/roptlib/cwrapper/blas/zhpr2.h>
#include <thirdparty/roptlib/cwrapper/blas/zrotg.h>
#include <thirdparty/roptlib/cwrapper/blas/zscal.h>
#include <thirdparty/roptlib/cwrapper/blas/zswap.h>
#include <thirdparty/roptlib/cwrapper/blas/zsymm.h>
#include <thirdparty/roptlib/cwrapper/blas/zsyr2k.h>
#include <thirdparty/roptlib/cwrapper/blas/zsyrk.h>
#include <thirdparty/roptlib/cwrapper/blas/ztbmv.h>
#include <thirdparty/roptlib/cwrapper/blas/ztbsv.h>
#include <thirdparty/roptlib/cwrapper/blas/ztpmv.h>
#include <thirdparty/roptlib/cwrapper/blas/ztpsv.h>
#include <thirdparty/roptlib/cwrapper/blas/ztrmm.h>
#include <thirdparty/roptlib/cwrapper/blas/ztrmv.h>
#include <thirdparty/roptlib/cwrapper/blas/ztrsm.h>
#include <thirdparty/roptlib/cwrapper/blas/ztrsv.h>

/* LAPACK */
#include <thirdparty/roptlib/cwrapper/lapack/cbdsqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgbbrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgbcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgbequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgbrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgbsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgbsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgbtf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgbtrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgbtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgebak.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgebal.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgebd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgebrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgecon.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgeequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgees.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgeesx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgeev.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgeevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgegs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgegv.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgehd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgehrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgelq2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgelqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgels.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgelsd.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgelss.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgelsx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgelsy.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgeql2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgeqlf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgeqp3.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgeqpf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgeqr2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgeqrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgerfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgerq2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgerqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgesc2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgesdd.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgesv.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgesvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgesvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgetc2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgetf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgetrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgetri.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgetrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cggbak.h>
#include <thirdparty/roptlib/cwrapper/lapack/cggbal.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgges.h>
#include <thirdparty/roptlib/cwrapper/lapack/cggesx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cggev.h>
#include <thirdparty/roptlib/cwrapper/lapack/cggevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cggglm.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgghrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgglse.h>
#include <thirdparty/roptlib/cwrapper/lapack/cggqrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cggrqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cggsvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/cggsvp.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgtcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgtrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgtsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgtsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgttrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgttrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cgtts2.h>
#include <thirdparty/roptlib/cwrapper/lapack/chbev.h>
#include <thirdparty/roptlib/cwrapper/lapack/chbevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/chbevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/chbgst.h>
#include <thirdparty/roptlib/cwrapper/lapack/chbgv.h>
#include <thirdparty/roptlib/cwrapper/lapack/chbgvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/chbgvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/chbtrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/checon.h>
#include <thirdparty/roptlib/cwrapper/lapack/cheev.h>
#include <thirdparty/roptlib/cwrapper/lapack/cheevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/cheevr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cheevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/chegs2.h>
#include <thirdparty/roptlib/cwrapper/lapack/chegst.h>
#include <thirdparty/roptlib/cwrapper/lapack/chegv.h>
#include <thirdparty/roptlib/cwrapper/lapack/chegvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/chegvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cherfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/chesv.h>
#include <thirdparty/roptlib/cwrapper/lapack/chesvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/chetd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/chetf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/chetrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/chetrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/chetri.h>
#include <thirdparty/roptlib/cwrapper/lapack/chetrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/chgeqz.h>
#include <thirdparty/roptlib/cwrapper/lapack/chpcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/chpev.h>
#include <thirdparty/roptlib/cwrapper/lapack/chpevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/chpevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/chpgst.h>
#include <thirdparty/roptlib/cwrapper/lapack/chpgv.h>
#include <thirdparty/roptlib/cwrapper/lapack/chpgvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/chpgvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/chprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/chpsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/chpsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/chptrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/chptrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/chptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/chptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/chsein.h>
#include <thirdparty/roptlib/cwrapper/lapack/chseqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/clabrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/clacgv.h>
#include <thirdparty/roptlib/cwrapper/lapack/clacon.h>
#include <thirdparty/roptlib/cwrapper/lapack/clacp2.h>
#include <thirdparty/roptlib/cwrapper/lapack/clacpy.h>
#include <thirdparty/roptlib/cwrapper/lapack/clacrm.h>
#include <thirdparty/roptlib/cwrapper/lapack/clacrt.h>
#include <thirdparty/roptlib/cwrapper/lapack/cladiv.h>
#include <thirdparty/roptlib/cwrapper/lapack/claed0.h>
#include <thirdparty/roptlib/cwrapper/lapack/claed7.h>
#include <thirdparty/roptlib/cwrapper/lapack/claed8.h>
#include <thirdparty/roptlib/cwrapper/lapack/claein.h>
#include <thirdparty/roptlib/cwrapper/lapack/claesy.h>
#include <thirdparty/roptlib/cwrapper/lapack/claev2.h>
#include <thirdparty/roptlib/cwrapper/lapack/clags2.h>
#include <thirdparty/roptlib/cwrapper/lapack/clagtm.h>
#include <thirdparty/roptlib/cwrapper/lapack/clahef.h>
#include <thirdparty/roptlib/cwrapper/lapack/clahqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/clahrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/claic1.h>
#include <thirdparty/roptlib/cwrapper/lapack/clals0.h>
#include <thirdparty/roptlib/cwrapper/lapack/clalsa.h>
#include <thirdparty/roptlib/cwrapper/lapack/clalsd.h>
#include <thirdparty/roptlib/cwrapper/lapack/clangb.h>
#include <thirdparty/roptlib/cwrapper/lapack/clange.h>
#include <thirdparty/roptlib/cwrapper/lapack/clangt.h>
#include <thirdparty/roptlib/cwrapper/lapack/clanhb.h>
#include <thirdparty/roptlib/cwrapper/lapack/clanhe.h>
#include <thirdparty/roptlib/cwrapper/lapack/clanhp.h>
#include <thirdparty/roptlib/cwrapper/lapack/clanhs.h>
#include <thirdparty/roptlib/cwrapper/lapack/clanht.h>
#include <thirdparty/roptlib/cwrapper/lapack/clansb.h>
#include <thirdparty/roptlib/cwrapper/lapack/clansp.h>
#include <thirdparty/roptlib/cwrapper/lapack/clansy.h>
#include <thirdparty/roptlib/cwrapper/lapack/clantb.h>
#include <thirdparty/roptlib/cwrapper/lapack/clantp.h>
#include <thirdparty/roptlib/cwrapper/lapack/clantr.h>
#include <thirdparty/roptlib/cwrapper/lapack/clapll.h>
#include <thirdparty/roptlib/cwrapper/lapack/clapmt.h>
#include <thirdparty/roptlib/cwrapper/lapack/claqgb.h>
#include <thirdparty/roptlib/cwrapper/lapack/claqge.h>
#include <thirdparty/roptlib/cwrapper/lapack/claqhb.h>
#include <thirdparty/roptlib/cwrapper/lapack/claqhe.h>
#include <thirdparty/roptlib/cwrapper/lapack/claqhp.h>
#include <thirdparty/roptlib/cwrapper/lapack/claqp2.h>
#include <thirdparty/roptlib/cwrapper/lapack/claqps.h>
#include <thirdparty/roptlib/cwrapper/lapack/claqsb.h>
#include <thirdparty/roptlib/cwrapper/lapack/claqsp.h>
#include <thirdparty/roptlib/cwrapper/lapack/claqsy.h>
#include <thirdparty/roptlib/cwrapper/lapack/clar1v.h>
#include <thirdparty/roptlib/cwrapper/lapack/clar2v.h>
#include <thirdparty/roptlib/cwrapper/lapack/clarcm.h>
#include <thirdparty/roptlib/cwrapper/lapack/clarf.h>
#include <thirdparty/roptlib/cwrapper/lapack/clarfb.h>
#include <thirdparty/roptlib/cwrapper/lapack/clarfg.h>
#include <thirdparty/roptlib/cwrapper/lapack/clarft.h>
#include <thirdparty/roptlib/cwrapper/lapack/clarfx.h>
#include <thirdparty/roptlib/cwrapper/lapack/clargv.h>
#include <thirdparty/roptlib/cwrapper/lapack/clarnv.h>
#include <thirdparty/roptlib/cwrapper/lapack/clarrv.h>
#include <thirdparty/roptlib/cwrapper/lapack/clartg.h>
#include <thirdparty/roptlib/cwrapper/lapack/clartv.h>
#include <thirdparty/roptlib/cwrapper/lapack/clarz.h>
#include <thirdparty/roptlib/cwrapper/lapack/clarzb.h>
#include <thirdparty/roptlib/cwrapper/lapack/clarzt.h>
#include <thirdparty/roptlib/cwrapper/lapack/clascl.h>
#include <thirdparty/roptlib/cwrapper/lapack/claset.h>
#include <thirdparty/roptlib/cwrapper/lapack/clasr.h>
#include <thirdparty/roptlib/cwrapper/lapack/classq.h>
#include <thirdparty/roptlib/cwrapper/lapack/claswp.h>
#include <thirdparty/roptlib/cwrapper/lapack/clasyf.h>
#include <thirdparty/roptlib/cwrapper/lapack/clatbs.h>
#include <thirdparty/roptlib/cwrapper/lapack/clatdf.h>
#include <thirdparty/roptlib/cwrapper/lapack/clatps.h>
#include <thirdparty/roptlib/cwrapper/lapack/clatrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/clatrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/clatrz.h>
#include <thirdparty/roptlib/cwrapper/lapack/clatzm.h>
#include <thirdparty/roptlib/cwrapper/lapack/clauu2.h>
#include <thirdparty/roptlib/cwrapper/lapack/clauum.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpbcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpbequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpbrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpbstf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpbsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpbsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpbtf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpbtrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpbtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpocon.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpoequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/cporfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cposv.h>
#include <thirdparty/roptlib/cwrapper/lapack/cposvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpotf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpotrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpotri.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpotrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cppcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/cppequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cppsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/cppsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpptrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cptcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpteqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cptrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cptsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/cptsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpttrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cpttrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cptts2.h>
#include <thirdparty/roptlib/cwrapper/lapack/crot.h>
#include <thirdparty/roptlib/cwrapper/lapack/cspcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/cspmv.h>
#include <thirdparty/roptlib/cwrapper/lapack/cspr.h>
#include <thirdparty/roptlib/cwrapper/lapack/csprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/cspsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/cspsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/csptrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/csptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/csptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/csrot.h>
#include <thirdparty/roptlib/cwrapper/lapack/csrscl.h>
#include <thirdparty/roptlib/cwrapper/lapack/cstedc.h>
#include <thirdparty/roptlib/cwrapper/lapack/cstegr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cstein.h>
#include <thirdparty/roptlib/cwrapper/lapack/csteqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/csycon.h>
#include <thirdparty/roptlib/cwrapper/lapack/csymv.h>
#include <thirdparty/roptlib/cwrapper/lapack/csyr.h>
#include <thirdparty/roptlib/cwrapper/lapack/csyrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/csysv.h>
#include <thirdparty/roptlib/cwrapper/lapack/csysvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/csytf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/csytrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/csytri.h>
#include <thirdparty/roptlib/cwrapper/lapack/csytrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctbcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctbrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctbtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctgevc.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctgex2.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctgexc.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctgsen.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctgsja.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctgsna.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctgsy2.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctgsyl.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctpcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctrcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctrevc.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctrexc.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctrrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctrsen.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctrsna.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctrsyl.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctrti2.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctrtri.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctrtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctzrqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/ctzrzf.h>
#include <thirdparty/roptlib/cwrapper/lapack/cung2l.h>
#include <thirdparty/roptlib/cwrapper/lapack/cung2r.h>
#include <thirdparty/roptlib/cwrapper/lapack/cungbr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunghr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cungl2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunglq.h>
#include <thirdparty/roptlib/cwrapper/lapack/cungql.h>
#include <thirdparty/roptlib/cwrapper/lapack/cungqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cungr2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cungrq.h>
#include <thirdparty/roptlib/cwrapper/lapack/cungtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunm2l.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunm2r.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunmbr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunmhr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunml2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunmlq.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunmql.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunmqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunmr2.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunmr3.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunmrq.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunmrz.h>
#include <thirdparty/roptlib/cwrapper/lapack/cunmtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cupgtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/cupmtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dbdsdc.h>
#include <thirdparty/roptlib/cwrapper/lapack/dbdsqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/ddisna.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgbbrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgbcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgbequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgbrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgbsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgbsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgbtf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgbtrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgbtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgebak.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgebal.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgebd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgebrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgecon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgeequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgees.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgeesx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgeev.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgeevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgegs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgegv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgehd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgehrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgelq2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgelqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgels.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgelsd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgelss.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgelsx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgelsy.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgeql2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgeqlf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgeqp3.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgeqpf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgeqr2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgeqrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgerfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgerq2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgerqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgesc2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgesdd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgesv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgesvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgesvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgetc2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgetf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgetrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgetri.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgetrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dggbak.h>
#include <thirdparty/roptlib/cwrapper/lapack/dggbal.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgges.h>
#include <thirdparty/roptlib/cwrapper/lapack/dggesx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dggev.h>
#include <thirdparty/roptlib/cwrapper/lapack/dggevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dggglm.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgghrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgglse.h>
#include <thirdparty/roptlib/cwrapper/lapack/dggqrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dggrqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dggsvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dggsvp.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgtcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgtrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgtsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgtsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgttrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgttrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dgtts2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dhgeqz.h>
#include <thirdparty/roptlib/cwrapper/lapack/dhsein.h>
#include <thirdparty/roptlib/cwrapper/lapack/dhseqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlabad.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlabrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlacon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlacpy.h>
#include <thirdparty/roptlib/cwrapper/lapack/dladiv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlae2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaebz.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaed0.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaed1.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaed2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaed3.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaed4.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaed5.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaed6.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaed7.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaed8.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaed9.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaeda.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaein.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaev2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaexc.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlag2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlags2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlagtf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlagtm.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlagts.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlagv2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlahqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlahrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaic1.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaln2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlals0.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlalsa.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlalsd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlamch.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlamrg.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlangb.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlange.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlangt.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlanhs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlansb.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlansp.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlanst.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlansy.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlantb.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlantp.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlantr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlanv2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlapll.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlapmt.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlapy2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlapy3.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaqgb.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaqge.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaqp2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaqps.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaqsb.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaqsp.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaqsy.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaqtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlar1v.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlar2v.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarfb.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarfg.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarft.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarfx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlargv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarnv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarrb.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarre.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarrv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlartg.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlartv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaruv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarz.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarzb.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlarzt.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlas2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlascl.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasd0.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasd1.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasd3.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasd4.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasd5.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasd6.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasd7.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasd8.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasd9.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasda.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasdq.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasdt.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaset.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasq1.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasq2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasq3.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasq4.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasq5.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasq6.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasrt.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlassq.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasv2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlaswp.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasy2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlasyf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlatbs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlatdf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlatps.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlatrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlatrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlatrz.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlatzm.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlauu2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dlauum.h>
#include <thirdparty/roptlib/cwrapper/lapack/dopgtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dopmtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorg2l.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorg2r.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorgbr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorghr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorgl2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorglq.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorgql.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorgqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorgr2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorgrq.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorgtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorm2l.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorm2r.h>
#include <thirdparty/roptlib/cwrapper/lapack/dormbr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dormhr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dorml2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dormlq.h>
#include <thirdparty/roptlib/cwrapper/lapack/dormql.h>
#include <thirdparty/roptlib/cwrapper/lapack/dormqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dormr2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dormr3.h>
#include <thirdparty/roptlib/cwrapper/lapack/dormrq.h>
#include <thirdparty/roptlib/cwrapper/lapack/dormrz.h>
#include <thirdparty/roptlib/cwrapper/lapack/dormtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpbcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpbequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpbrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpbstf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpbsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpbsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpbtf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpbtrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpbtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpocon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpoequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/dporfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dposv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dposvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpotf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpotrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpotri.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpotrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dppcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dppequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dppsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dppsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpptrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dptcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpteqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dptrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dptsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dptsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpttrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dpttrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dptts2.h>
#include <thirdparty/roptlib/cwrapper/lapack/drscl.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsbev.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsbevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsbevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsbgst.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsbgv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsbgvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsbgvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsbtrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsecnd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dspcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dspev.h>
#include <thirdparty/roptlib/cwrapper/lapack/dspevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dspevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dspgst.h>
#include <thirdparty/roptlib/cwrapper/lapack/dspgv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dspgvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dspgvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dspsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dspsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsptrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsptrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dstebz.h>
#include <thirdparty/roptlib/cwrapper/lapack/dstedc.h>
#include <thirdparty/roptlib/cwrapper/lapack/dstegr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dstein.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsteqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsterf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dstev.h>
#include <thirdparty/roptlib/cwrapper/lapack/dstevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dstevr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dstevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsycon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsyev.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsyevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsyevr.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsyevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsygs2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsygst.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsygv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsygvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsygvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsyrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsysv.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsysvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsytd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsytf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsytrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsytrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsytri.h>
#include <thirdparty/roptlib/cwrapper/lapack/dsytrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtbcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtbrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtbtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtgevc.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtgex2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtgexc.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtgsen.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtgsja.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtgsna.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtgsy2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtgsyl.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtpcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtrcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtrevc.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtrexc.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtrrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtrsen.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtrsna.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtrsyl.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtrti2.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtrtri.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtrtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtzrqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dtzrzf.h>
#include <thirdparty/roptlib/cwrapper/lapack/dzsum1.h>
#include <thirdparty/roptlib/cwrapper/lapack/f2c.h>
#include <thirdparty/roptlib/cwrapper/lapack/icmax1.h>
#include <thirdparty/roptlib/cwrapper/lapack/ieeeck.h>
#include <thirdparty/roptlib/cwrapper/lapack/ilaenv.h>
#include <thirdparty/roptlib/cwrapper/lapack/izmax1.h>
#include <thirdparty/roptlib/cwrapper/lapack/lsame.h>
#include <thirdparty/roptlib/cwrapper/lapack/lsamen.h>
#include <thirdparty/roptlib/cwrapper/lapack/sbdsdc.h>
#include <thirdparty/roptlib/cwrapper/lapack/sbdsqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/scsum1.h>
#include <thirdparty/roptlib/cwrapper/lapack/sdisna.h>
#include <thirdparty/roptlib/cwrapper/lapack/second.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgbbrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgbcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgbequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgbrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgbsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgbsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgbtf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgbtrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgbtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgebak.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgebal.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgebd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgebrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgecon.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgeequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgees.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgeesx.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgeev.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgeevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgegs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgegv.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgehd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgehrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgelq2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgelqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgels.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgelsd.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgelss.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgelsx.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgelsy.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgeql2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgeqlf.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgeqp3.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgeqpf.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgeqr2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgeqrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgerfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgerq2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgerqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgesc2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgesdd.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgesv.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgesvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgesvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgetc2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgetf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgetrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgetri.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgetrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sggbak.h>
#include <thirdparty/roptlib/cwrapper/lapack/sggbal.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgges.h>
#include <thirdparty/roptlib/cwrapper/lapack/sggesx.h>
#include <thirdparty/roptlib/cwrapper/lapack/sggev.h>
#include <thirdparty/roptlib/cwrapper/lapack/sggevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/sggglm.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgghrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgglse.h>
#include <thirdparty/roptlib/cwrapper/lapack/sggqrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/sggrqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/sggsvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/sggsvp.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgtcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgtrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgtsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgtsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgttrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgttrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sgtts2.h>
#include <thirdparty/roptlib/cwrapper/lapack/shgeqz.h>
#include <thirdparty/roptlib/cwrapper/lapack/shsein.h>
#include <thirdparty/roptlib/cwrapper/lapack/shseqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/slabad.h>
#include <thirdparty/roptlib/cwrapper/lapack/slabrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/slacon.h>
#include <thirdparty/roptlib/cwrapper/lapack/slacpy.h>
#include <thirdparty/roptlib/cwrapper/lapack/sladiv.h>
#include <thirdparty/roptlib/cwrapper/lapack/slae2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaebz.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaed0.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaed1.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaed2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaed3.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaed4.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaed5.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaed6.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaed7.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaed8.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaed9.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaeda.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaein.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaev2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaexc.h>
#include <thirdparty/roptlib/cwrapper/lapack/slag2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slags2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slagtf.h>
#include <thirdparty/roptlib/cwrapper/lapack/slagtm.h>
#include <thirdparty/roptlib/cwrapper/lapack/slagts.h>
#include <thirdparty/roptlib/cwrapper/lapack/slagv2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slahqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/slahrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaic1.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaln2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slals0.h>
#include <thirdparty/roptlib/cwrapper/lapack/slalsa.h>
#include <thirdparty/roptlib/cwrapper/lapack/slalsd.h>
#include <thirdparty/roptlib/cwrapper/lapack/slamch.h>
#include <thirdparty/roptlib/cwrapper/lapack/slamrg.h>
#include <thirdparty/roptlib/cwrapper/lapack/slangb.h>
#include <thirdparty/roptlib/cwrapper/lapack/slange.h>
#include <thirdparty/roptlib/cwrapper/lapack/slangt.h>
#include <thirdparty/roptlib/cwrapper/lapack/slanhs.h>
#include <thirdparty/roptlib/cwrapper/lapack/slansb.h>
#include <thirdparty/roptlib/cwrapper/lapack/slansp.h>
#include <thirdparty/roptlib/cwrapper/lapack/slanst.h>
#include <thirdparty/roptlib/cwrapper/lapack/slansy.h>
#include <thirdparty/roptlib/cwrapper/lapack/slantb.h>
#include <thirdparty/roptlib/cwrapper/lapack/slantp.h>
#include <thirdparty/roptlib/cwrapper/lapack/slantr.h>
#include <thirdparty/roptlib/cwrapper/lapack/slanv2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slapll.h>
#include <thirdparty/roptlib/cwrapper/lapack/slapmt.h>
#include <thirdparty/roptlib/cwrapper/lapack/slapy2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slapy3.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaqgb.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaqge.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaqp2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaqps.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaqsb.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaqsp.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaqsy.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaqtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/slar1v.h>
#include <thirdparty/roptlib/cwrapper/lapack/slar2v.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarf.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarfb.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarfg.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarft.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarfx.h>
#include <thirdparty/roptlib/cwrapper/lapack/slargv.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarnv.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarrb.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarre.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarrv.h>
#include <thirdparty/roptlib/cwrapper/lapack/slartg.h>
#include <thirdparty/roptlib/cwrapper/lapack/slartv.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaruv.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarz.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarzb.h>
#include <thirdparty/roptlib/cwrapper/lapack/slarzt.h>
#include <thirdparty/roptlib/cwrapper/lapack/slas2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slascl.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasd0.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasd1.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasd3.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasd4.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasd5.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasd6.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasd7.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasd8.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasd9.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasda.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasdq.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasdt.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaset.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasq1.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasq2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasq3.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasq4.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasq5.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasq6.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasr.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasrt.h>
#include <thirdparty/roptlib/cwrapper/lapack/slassq.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasv2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slaswp.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasy2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slasyf.h>
#include <thirdparty/roptlib/cwrapper/lapack/slatbs.h>
#include <thirdparty/roptlib/cwrapper/lapack/slatdf.h>
#include <thirdparty/roptlib/cwrapper/lapack/slatps.h>
#include <thirdparty/roptlib/cwrapper/lapack/slatrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/slatrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/slatrz.h>
#include <thirdparty/roptlib/cwrapper/lapack/slatzm.h>
#include <thirdparty/roptlib/cwrapper/lapack/slauu2.h>
#include <thirdparty/roptlib/cwrapper/lapack/slauum.h>
#include <thirdparty/roptlib/cwrapper/lapack/sopgtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/sopmtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorg2l.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorg2r.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorgbr.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorghr.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorgl2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorglq.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorgql.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorgqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorgr2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorgrq.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorgtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorm2l.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorm2r.h>
#include <thirdparty/roptlib/cwrapper/lapack/sormbr.h>
#include <thirdparty/roptlib/cwrapper/lapack/sormhr.h>
#include <thirdparty/roptlib/cwrapper/lapack/sorml2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sormlq.h>
#include <thirdparty/roptlib/cwrapper/lapack/sormql.h>
#include <thirdparty/roptlib/cwrapper/lapack/sormqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/sormr2.h>
#include <thirdparty/roptlib/cwrapper/lapack/sormr3.h>
#include <thirdparty/roptlib/cwrapper/lapack/sormrq.h>
#include <thirdparty/roptlib/cwrapper/lapack/sormrz.h>
#include <thirdparty/roptlib/cwrapper/lapack/sormtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/spbcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/spbequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/spbrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/spbstf.h>
#include <thirdparty/roptlib/cwrapper/lapack/spbsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/spbsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/spbtf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/spbtrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/spbtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/spocon.h>
#include <thirdparty/roptlib/cwrapper/lapack/spoequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/sporfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sposv.h>
#include <thirdparty/roptlib/cwrapper/lapack/sposvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/spotf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/spotrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/spotri.h>
#include <thirdparty/roptlib/cwrapper/lapack/spotrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sppcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/sppequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/spprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sppsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/sppsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/spptrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/spptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/spptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sptcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/spteqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/sptrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sptsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/sptsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/spttrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/spttrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sptts2.h>
#include <thirdparty/roptlib/cwrapper/lapack/srscl.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssbev.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssbevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssbevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssbgst.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssbgv.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssbgvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssbgvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssbtrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/sspcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/sspev.h>
#include <thirdparty/roptlib/cwrapper/lapack/sspevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/sspevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/sspgst.h>
#include <thirdparty/roptlib/cwrapper/lapack/sspgv.h>
#include <thirdparty/roptlib/cwrapper/lapack/sspgvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/sspgvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sspsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/sspsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssptrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssptrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/sstebz.h>
#include <thirdparty/roptlib/cwrapper/lapack/sstedc.h>
#include <thirdparty/roptlib/cwrapper/lapack/sstegr.h>
#include <thirdparty/roptlib/cwrapper/lapack/sstein.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssteqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssterf.h>
#include <thirdparty/roptlib/cwrapper/lapack/sstev.h>
#include <thirdparty/roptlib/cwrapper/lapack/sstevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/sstevr.h>
#include <thirdparty/roptlib/cwrapper/lapack/sstevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssycon.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssyev.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssyevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssyevr.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssyevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssygs2.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssygst.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssygv.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssygvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssygvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssyrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssysv.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssysvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssytd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssytf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssytrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssytrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssytri.h>
#include <thirdparty/roptlib/cwrapper/lapack/ssytrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/stbcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/stbrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/stbtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/stgevc.h>
#include <thirdparty/roptlib/cwrapper/lapack/stgex2.h>
#include <thirdparty/roptlib/cwrapper/lapack/stgexc.h>
#include <thirdparty/roptlib/cwrapper/lapack/stgsen.h>
#include <thirdparty/roptlib/cwrapper/lapack/stgsja.h>
#include <thirdparty/roptlib/cwrapper/lapack/stgsna.h>
#include <thirdparty/roptlib/cwrapper/lapack/stgsy2.h>
#include <thirdparty/roptlib/cwrapper/lapack/stgsyl.h>
#include <thirdparty/roptlib/cwrapper/lapack/stpcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/stprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/stptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/stptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/strcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/strevc.h>
#include <thirdparty/roptlib/cwrapper/lapack/strexc.h>
#include <thirdparty/roptlib/cwrapper/lapack/strrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/strsen.h>
#include <thirdparty/roptlib/cwrapper/lapack/strsna.h>
#include <thirdparty/roptlib/cwrapper/lapack/strsyl.h>
#include <thirdparty/roptlib/cwrapper/lapack/strti2.h>
#include <thirdparty/roptlib/cwrapper/lapack/strtri.h>
#include <thirdparty/roptlib/cwrapper/lapack/strtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/stzrqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/stzrzf.h>
#include <thirdparty/roptlib/cwrapper/lapack/xerbla.h>
#include <thirdparty/roptlib/cwrapper/lapack/zbdsqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zdrot.h>
#include <thirdparty/roptlib/cwrapper/lapack/zdrscl.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgbbrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgbcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgbequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgbrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgbsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgbsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgbtf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgbtrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgbtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgebak.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgebal.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgebd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgebrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgecon.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgeequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgees.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgeesx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgeev.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgeevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgegs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgegv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgehd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgehrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgelq2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgelqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgels.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgelsd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgelss.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgelsx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgelsy.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgeql2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgeqlf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgeqp3.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgeqpf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgeqr2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgeqrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgerfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgerq2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgerqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgesc2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgesdd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgesv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgesvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgesvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgetc2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgetf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgetrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgetri.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgetrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zggbak.h>
#include <thirdparty/roptlib/cwrapper/lapack/zggbal.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgges.h>
#include <thirdparty/roptlib/cwrapper/lapack/zggesx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zggev.h>
#include <thirdparty/roptlib/cwrapper/lapack/zggevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zggglm.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgghrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgglse.h>
#include <thirdparty/roptlib/cwrapper/lapack/zggqrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zggrqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zggsvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zggsvp.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgtcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgtrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgtsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgtsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgttrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgttrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zgtts2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhbev.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhbevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhbevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhbgst.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhbgv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhbgvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhbgvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhbtrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhecon.h>
#include <thirdparty/roptlib/cwrapper/lapack/zheev.h>
#include <thirdparty/roptlib/cwrapper/lapack/zheevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zheevr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zheevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhegs2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhegst.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhegv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhegvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhegvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zherfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhesv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhesvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhetd2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhetf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhetrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhetrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhetri.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhetrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhgeqz.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhpcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhpev.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhpevd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhpevx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhpgst.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhpgv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhpgvd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhpgvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhpsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhpsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhptrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhptrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhsein.h>
#include <thirdparty/roptlib/cwrapper/lapack/zhseqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlabrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlacgv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlacon.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlacp2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlacpy.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlacrm.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlacrt.h>
#include <thirdparty/roptlib/cwrapper/lapack/zladiv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaed0.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaed7.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaed8.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaein.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaesy.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaev2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlags2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlagtm.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlahef.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlahqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlahrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaic1.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlals0.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlalsa.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlalsd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlangb.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlange.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlangt.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlanhb.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlanhe.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlanhp.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlanhs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlanht.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlansb.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlansp.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlansy.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlantb.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlantp.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlantr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlapll.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlapmt.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaqgb.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaqge.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaqhb.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaqhe.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaqhp.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaqp2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaqps.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaqsb.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaqsp.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaqsy.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlar1v.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlar2v.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlarcm.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlarf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlarfb.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlarfg.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlarft.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlarfx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlargv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlarnv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlarrv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlartg.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlartv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlarz.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlarzb.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlarzt.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlascl.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaset.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlasr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlassq.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlaswp.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlasyf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlatbs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlatdf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlatps.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlatrd.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlatrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlatrz.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlatzm.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlauu2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zlauum.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpbcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpbequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpbrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpbstf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpbsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpbsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpbtf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpbtrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpbtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpocon.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpoequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/zporfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zposv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zposvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpotf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpotrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpotri.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpotrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zppcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/zppequ.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zppsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zppsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpptrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zptcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpteqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zptrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zptsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zptsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpttrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zpttrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zptts2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zrot.h>
#include <thirdparty/roptlib/cwrapper/lapack/zspcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/zspmv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zspr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zspsv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zspsvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsptrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zstedc.h>
#include <thirdparty/roptlib/cwrapper/lapack/zstegr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zstein.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsteqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsycon.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsymv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsyr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsyrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsysv.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsysvx.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsytf2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsytrf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsytri.h>
#include <thirdparty/roptlib/cwrapper/lapack/zsytrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztbcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztbrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztbtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztgevc.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztgex2.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztgexc.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztgsen.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztgsja.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztgsna.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztgsy2.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztgsyl.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztpcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztprfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztptri.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztptrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztrcon.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztrevc.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztrexc.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztrrfs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztrsen.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztrsna.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztrsyl.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztrti2.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztrtri.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztrtrs.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztzrqf.h>
#include <thirdparty/roptlib/cwrapper/lapack/ztzrzf.h>
#include <thirdparty/roptlib/cwrapper/lapack/zung2l.h>
#include <thirdparty/roptlib/cwrapper/lapack/zung2r.h>
#include <thirdparty/roptlib/cwrapper/lapack/zungbr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunghr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zungl2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunglq.h>
#include <thirdparty/roptlib/cwrapper/lapack/zungql.h>
#include <thirdparty/roptlib/cwrapper/lapack/zungqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zungr2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zungrq.h>
#include <thirdparty/roptlib/cwrapper/lapack/zungtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunm2l.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunm2r.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunmbr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunmhr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunml2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunmlq.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunmql.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunmqr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunmr2.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunmr3.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunmrq.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunmrz.h>
#include <thirdparty/roptlib/cwrapper/lapack/zunmtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zupgtr.h>
#include <thirdparty/roptlib/cwrapper/lapack/zupmtr.h>
#endif /* end of ifndef MATLAB_MEX_FILE */
#ifdef MATLAB_MEX_FILE

#include "mex.h"
#include "blas.h"
#include "lapack.h"

#define caxpy_ caxpy
#define ccopy_ ccopy
#define cdotc_ cdotc
#define cdotu_ cdotu
#define cgbmv_ cgbmv
#define cgemm_ cgemm
#define cgemv_ cgemv
#define cgerc_ cgerc
#define cgeru_ cgeru
#define chbmv_ chbmv
#define chemm_ chemm
#define chemv_ chemv
#define cher_ cher
#define cher2_ cher2
#define cher2k_ cher2k
#define cherk_ cherk
#define chpmv_ chpmv
#define chpr_ chpr
#define chpr2_ chpr2
#define crotg_ crotg
#define cscal_ cscal
#define csscal_ csscal
#define cswap_ cswap
#define csymm_ csymm
#define csyr2k_ csyr2k
#define csyrk_ csyrk
#define ctbmv_ ctbmv
#define ctbsv_ ctbsv
#define ctpmv_ ctpmv
#define ctpsv_ ctpsv
#define ctrmm_ ctrmm
#define ctrmv_ ctrmv
#define ctrsm_ ctrsm
#define ctrsv_ ctrsv
#define dasum_ dasum
#define daxpy_ daxpy
#define dcabs1_ dcabs1
#define dcopy_ dcopy
#define ddot_ ddot
#define dgbmv_ dgbmv
#define dgemm_ dgemm
#define dgemv_ dgemv
#define dger_ dger
#define dnrm2_ dnrm2
#define drot_ drot
#define drotg_ drotg
#define dsbmv_ dsbmv
#define dscal_ dscal
#define dspmv_ dspmv
#define dspr_ dspr
#define dspr2_ dspr2
#define dswap_ dswap
#define dsymm_ dsymm
#define dsymv_ dsymv
#define dsyr_ dsyr
#define dsyr2_ dsyr2
#define dsyr2k_ dsyr2k
#define dsyrk_ dsyrk
#define dtbmv_ dtbmv
#define dtbsv_ dtbsv
#define dtpmv_ dtpmv
#define dtpsv_ dtpsv
#define dtrmm_ dtrmm
#define dtrmv_ dtrmv
#define dtrsm_ dtrsm
#define dtrsv_ dtrsv
#define dzasum_ dzasum
#define dznrm2_ dznrm2
#define f2c_ f2c
#define icamax_ icamax
#define idamax_ idamax
#define isamax_ isamax
#define izamax_ izamax
#define lsame_ lsame
#define sasum_ sasum
#define saxpy_ saxpy
#define scasum_ scasum
#define scnrm2_ scnrm2
#define scopy_ scopy
#define sdot_ sdot
#define sgbmv_ sgbmv
#define sgemm_ sgemm
#define sgemv_ sgemv
#define sger_ sger
#define snrm2_ snrm2
#define srot_ srot
#define srotg_ srotg
#define ssbmv_ ssbmv
#define sscal_ sscal
#define sspmv_ sspmv
#define sspr_ sspr
#define sspr2_ sspr2
#define sswap_ sswap
#define ssymm_ ssymm
#define ssymv_ ssymv
#define ssyr_ ssyr
#define ssyr2_ ssyr2
#define ssyr2k_ ssyr2k
#define ssyrk_ ssyrk
#define stbmv_ stbmv
#define stbsv_ stbsv
#define stpmv_ stpmv
#define stpsv_ stpsv
#define strmm_ strmm
#define strmv_ strmv
#define strsm_ strsm
#define strsv_ strsv
#define xerbla_ xerbla
#define zaxpy_ zaxpy
#define zcopy_ zcopy
#define zdotc_ zdotc
#define zdotu_ zdotu
#define zdscal_ zdscal
#define zgbmv_ zgbmv
#define zgemm_ zgemm
#define zgemv_ zgemv
#define zgerc_ zgerc
#define zgeru_ zgeru
#define zhbmv_ zhbmv
#define zhemm_ zhemm
#define zhemv_ zhemv
#define zher_ zher
#define zher2_ zher2
#define zher2k_ zher2k
#define zherk_ zherk
#define zhpmv_ zhpmv
#define zhpr_ zhpr
#define zhpr2_ zhpr2
#define zrotg_ zrotg
#define zscal_ zscal
#define zswap_ zswap
#define zsymm_ zsymm
#define zsyr2k_ zsyr2k
#define zsyrk_ zsyrk
#define ztbmv_ ztbmv
#define ztbsv_ ztbsv
#define ztpmv_ ztpmv
#define ztpsv_ ztpsv
#define ztrmm_ ztrmm
#define ztrmv_ ztrmv
#define ztrsm_ ztrsm
#define ztrsv_ ztrsv

/* LAPACK */
#define cbdsqr_ cbdsqr
#define cgbbrd_ cgbbrd
#define cgbcon_ cgbcon
#define cgbequ_ cgbequ
#define cgbrfs_ cgbrfs
#define cgbsv_ cgbsv
#define cgbsvx_ cgbsvx
#define cgbtf2_ cgbtf2
#define cgbtrf_ cgbtrf
#define cgbtrs_ cgbtrs
#define cgebak_ cgebak
#define cgebal_ cgebal
#define cgebd2_ cgebd2
#define cgebrd_ cgebrd
#define cgecon_ cgecon
#define cgeequ_ cgeequ
#define cgees_ cgees
#define cgeesx_ cgeesx
#define cgeev_ cgeev
#define cgeevx_ cgeevx
#define cgegs_ cgegs
#define cgegv_ cgegv
#define cgehd2_ cgehd2
#define cgehrd_ cgehrd
#define cgelq2_ cgelq2
#define cgelqf_ cgelqf
#define cgels_ cgels
#define cgelsd_ cgelsd
#define cgelss_ cgelss
#define cgelsx_ cgelsx
#define cgelsy_ cgelsy
#define cgeql2_ cgeql2
#define cgeqlf_ cgeqlf
#define cgeqp3_ cgeqp3
#define cgeqpf_ cgeqpf
#define cgeqr2_ cgeqr2
#define cgeqrf_ cgeqrf
#define cgerfs_ cgerfs
#define cgerq2_ cgerq2
#define cgerqf_ cgerqf
#define cgesc2_ cgesc2
#define cgesdd_ cgesdd
#define cgesv_ cgesv
#define cgesvd_ cgesvd
#define cgesvx_ cgesvx
#define cgetc2_ cgetc2
#define cgetf2_ cgetf2
#define cgetrf_ cgetrf
#define cgetri_ cgetri
#define cgetrs_ cgetrs
#define cggbak_ cggbak
#define cggbal_ cggbal
#define cgges_ cgges
#define cggesx_ cggesx
#define cggev_ cggev
#define cggevx_ cggevx
#define cggglm_ cggglm
#define cgghrd_ cgghrd
#define cgglse_ cgglse
#define cggqrf_ cggqrf
#define cggrqf_ cggrqf
#define cggsvd_ cggsvd
#define cggsvp_ cggsvp
#define cgtcon_ cgtcon
#define cgtrfs_ cgtrfs
#define cgtsv_ cgtsv
#define cgtsvx_ cgtsvx
#define cgttrf_ cgttrf
#define cgttrs_ cgttrs
#define cgtts2_ cgtts2
#define chbev_ chbev
#define chbevd_ chbevd
#define chbevx_ chbevx
#define chbgst_ chbgst
#define chbgv_ chbgv
#define chbgvd_ chbgvd
#define chbgvx_ chbgvx
#define chbtrd_ chbtrd
#define checon_ checon
#define cheev_ cheev
#define cheevd_ cheevd
#define cheevr_ cheevr
#define cheevx_ cheevx
#define chegs2_ chegs2
#define chegst_ chegst
#define chegv_ chegv
#define chegvd_ chegvd
#define chegvx_ chegvx
#define cherfs_ cherfs
#define chesv_ chesv
#define chesvx_ chesvx
#define chetd2_ chetd2
#define chetf2_ chetf2
#define chetrd_ chetrd
#define chetrf_ chetrf
#define chetri_ chetri
#define chetrs_ chetrs
#define chgeqz_ chgeqz
#define chpcon_ chpcon
#define chpev_ chpev
#define chpevd_ chpevd
#define chpevx_ chpevx
#define chpgst_ chpgst
#define chpgv_ chpgv
#define chpgvd_ chpgvd
#define chpgvx_ chpgvx
#define chprfs_ chprfs
#define chpsv_ chpsv
#define chpsvx_ chpsvx
#define chptrd_ chptrd
#define chptrf_ chptrf
#define chptri_ chptri
#define chptrs_ chptrs
#define chsein_ chsein
#define chseqr_ chseqr
#define clabrd_ clabrd
#define clacgv_ clacgv
#define clacon_ clacon
#define clacp2_ clacp2
#define clacpy_ clacpy
#define clacrm_ clacrm
#define clacrt_ clacrt
#define cladiv_ cladiv
#define claed0_ claed0
#define claed7_ claed7
#define claed8_ claed8
#define claein_ claein
#define claesy_ claesy
#define claev2_ claev2
#define clags2_ clags2
#define clagtm_ clagtm
#define clahef_ clahef
#define clahqr_ clahqr
#define clahrd_ clahrd
#define claic1_ claic1
#define clals0_ clals0
#define clalsa_ clalsa
#define clalsd_ clalsd
#define clangb_ clangb
#define clange_ clange
#define clangt_ clangt
#define clanhb_ clanhb
#define clanhe_ clanhe
#define clanhp_ clanhp
#define clanhs_ clanhs
#define clanht_ clanht
#define clansb_ clansb
#define clansp_ clansp
#define clansy_ clansy
#define clantb_ clantb
#define clantp_ clantp
#define clantr_ clantr
#define clapll_ clapll
#define clapmt_ clapmt
#define claqgb_ claqgb
#define claqge_ claqge
#define claqhb_ claqhb
#define claqhe_ claqhe
#define claqhp_ claqhp
#define claqp2_ claqp2
#define claqps_ claqps
#define claqsb_ claqsb
#define claqsp_ claqsp
#define claqsy_ claqsy
#define clar1v_ clar1v
#define clar2v_ clar2v
#define clarcm_ clarcm
#define clarf_ clarf
#define clarfb_ clarfb
#define clarfg_ clarfg
#define clarft_ clarft
#define clarfx_ clarfx
#define clargv_ clargv
#define clarnv_ clarnv
#define clarrv_ clarrv
#define clartg_ clartg
#define clartv_ clartv
#define clarz_ clarz
#define clarzb_ clarzb
#define clarzt_ clarzt
#define clascl_ clascl
#define claset_ claset
#define clasr_ clasr
#define classq_ classq
#define claswp_ claswp
#define clasyf_ clasyf
#define clatbs_ clatbs
#define clatdf_ clatdf
#define clatps_ clatps
#define clatrd_ clatrd
#define clatrs_ clatrs
#define clatrz_ clatrz
#define clatzm_ clatzm
#define clauu2_ clauu2
#define clauum_ clauum
#define cpbcon_ cpbcon
#define cpbequ_ cpbequ
#define cpbrfs_ cpbrfs
#define cpbstf_ cpbstf
#define cpbsv_ cpbsv
#define cpbsvx_ cpbsvx
#define cpbtf2_ cpbtf2
#define cpbtrf_ cpbtrf
#define cpbtrs_ cpbtrs
#define cpocon_ cpocon
#define cpoequ_ cpoequ
#define cporfs_ cporfs
#define cposv_ cposv
#define cposvx_ cposvx
#define cpotf2_ cpotf2
#define cpotrf_ cpotrf
#define cpotri_ cpotri
#define cpotrs_ cpotrs
#define cppcon_ cppcon
#define cppequ_ cppequ
#define cpprfs_ cpprfs
#define cppsv_ cppsv
#define cppsvx_ cppsvx
#define cpptrf_ cpptrf
#define cpptri_ cpptri
#define cpptrs_ cpptrs
#define cptcon_ cptcon
#define cpteqr_ cpteqr
#define cptrfs_ cptrfs
#define cptsv_ cptsv
#define cptsvx_ cptsvx
#define cpttrf_ cpttrf
#define cpttrs_ cpttrs
#define cptts2_ cptts2
#define crot_ crot
#define cspcon_ cspcon
#define cspmv_ cspmv
#define cspr_ cspr
#define csprfs_ csprfs
#define cspsv_ cspsv
#define cspsvx_ cspsvx
#define csptrf_ csptrf
#define csptri_ csptri
#define csptrs_ csptrs
#define csrot_ csrot
#define csrscl_ csrscl
#define cstedc_ cstedc
#define cstegr_ cstegr
#define cstein_ cstein
#define csteqr_ csteqr
#define csycon_ csycon
#define csymv_ csymv
#define csyr_ csyr
#define csyrfs_ csyrfs
#define csysv_ csysv
#define csysvx_ csysvx
#define csytf2_ csytf2
#define csytrf_ csytrf
#define csytri_ csytri
#define csytrs_ csytrs
#define ctbcon_ ctbcon
#define ctbrfs_ ctbrfs
#define ctbtrs_ ctbtrs
#define ctgevc_ ctgevc
#define ctgex2_ ctgex2
#define ctgexc_ ctgexc
#define ctgsen_ ctgsen
#define ctgsja_ ctgsja
#define ctgsna_ ctgsna
#define ctgsy2_ ctgsy2
#define ctgsyl_ ctgsyl
#define ctpcon_ ctpcon
#define ctprfs_ ctprfs
#define ctptri_ ctptri
#define ctptrs_ ctptrs
#define ctrcon_ ctrcon
#define ctrevc_ ctrevc
#define ctrexc_ ctrexc
#define ctrrfs_ ctrrfs
#define ctrsen_ ctrsen
#define ctrsna_ ctrsna
#define ctrsyl_ ctrsyl
#define ctrti2_ ctrti2
#define ctrtri_ ctrtri
#define ctrtrs_ ctrtrs
#define ctzrqf_ ctzrqf
#define ctzrzf_ ctzrzf
#define cung2l_ cung2l
#define cung2r_ cung2r
#define cungbr_ cungbr
#define cunghr_ cunghr
#define cungl2_ cungl2
#define cunglq_ cunglq
#define cungql_ cungql
#define cungqr_ cungqr
#define cungr2_ cungr2
#define cungrq_ cungrq
#define cungtr_ cungtr
#define cunm2l_ cunm2l
#define cunm2r_ cunm2r
#define cunmbr_ cunmbr
#define cunmhr_ cunmhr
#define cunml2_ cunml2
#define cunmlq_ cunmlq
#define cunmql_ cunmql
#define cunmqr_ cunmqr
#define cunmr2_ cunmr2
#define cunmr3_ cunmr3
#define cunmrq_ cunmrq
#define cunmrz_ cunmrz
#define cunmtr_ cunmtr
#define cupgtr_ cupgtr
#define cupmtr_ cupmtr
#define dbdsdc_ dbdsdc
#define dbdsqr_ dbdsqr
#define ddisna_ ddisna
#define dgbbrd_ dgbbrd
#define dgbcon_ dgbcon
#define dgbequ_ dgbequ
#define dgbrfs_ dgbrfs
#define dgbsv_ dgbsv
#define dgbsvx_ dgbsvx
#define dgbtf2_ dgbtf2
#define dgbtrf_ dgbtrf
#define dgbtrs_ dgbtrs
#define dgebak_ dgebak
#define dgebal_ dgebal
#define dgebd2_ dgebd2
#define dgebrd_ dgebrd
#define dgecon_ dgecon
#define dgeequ_ dgeequ
#define dgees_ dgees
#define dgeesx_ dgeesx
#define dgeev_ dgeev
#define dgeevx_ dgeevx
#define dgegs_ dgegs
#define dgegv_ dgegv
#define dgehd2_ dgehd2
#define dgehrd_ dgehrd
#define dgelq2_ dgelq2
#define dgelqf_ dgelqf
#define dgels_ dgels
#define dgelsd_ dgelsd
#define dgelss_ dgelss
#define dgelsx_ dgelsx
#define dgelsy_ dgelsy
#define dgeql2_ dgeql2
#define dgeqlf_ dgeqlf
#define dgeqp3_ dgeqp3
#define dgeqpf_ dgeqpf
#define dgeqr2_ dgeqr2
#define dgeqrf_ dgeqrf
#define dgerfs_ dgerfs
#define dgerq2_ dgerq2
#define dgerqf_ dgerqf
#define dgesc2_ dgesc2
#define dgesdd_ dgesdd
#define dgesv_ dgesv
#define dgesvd_ dgesvd
#define dgesvx_ dgesvx
#define dgetc2_ dgetc2
#define dgetf2_ dgetf2
#define dgetrf_ dgetrf
#define dgetri_ dgetri
#define dgetrs_ dgetrs
#define dggbak_ dggbak
#define dggbal_ dggbal
#define dgges_ dgges
#define dggesx_ dggesx
#define dggev_ dggev
#define dggevx_ dggevx
#define dggglm_ dggglm
#define dgghrd_ dgghrd
#define dgglse_ dgglse
#define dggqrf_ dggqrf
#define dggrqf_ dggrqf
#define dggsvd_ dggsvd
#define dggsvp_ dggsvp
#define dgtcon_ dgtcon
#define dgtrfs_ dgtrfs
#define dgtsv_ dgtsv
#define dgtsvx_ dgtsvx
#define dgttrf_ dgttrf
#define dgttrs_ dgttrs
#define dgtts2_ dgtts2
#define dhgeqz_ dhgeqz
#define dhsein_ dhsein
#define dhseqr_ dhseqr
#define dlabad_ dlabad
#define dlabrd_ dlabrd
#define dlacon_ dlacon
#define dlacpy_ dlacpy
#define dladiv_ dladiv
#define dlae2_ dlae2
#define dlaebz_ dlaebz
#define dlaed0_ dlaed0
#define dlaed1_ dlaed1
#define dlaed2_ dlaed2
#define dlaed3_ dlaed3
#define dlaed4_ dlaed4
#define dlaed5_ dlaed5
#define dlaed6_ dlaed6
#define dlaed7_ dlaed7
#define dlaed8_ dlaed8
#define dlaed9_ dlaed9
#define dlaeda_ dlaeda
#define dlaein_ dlaein
#define dlaev2_ dlaev2
#define dlaexc_ dlaexc
#define dlag2_ dlag2
#define dlags2_ dlags2
#define dlagtf_ dlagtf
#define dlagtm_ dlagtm
#define dlagts_ dlagts
#define dlagv2_ dlagv2
#define dlahqr_ dlahqr
#define dlahrd_ dlahrd
#define dlaic1_ dlaic1
#define dlaln2_ dlaln2
#define dlals0_ dlals0
#define dlalsa_ dlalsa
#define dlalsd_ dlalsd
#define dlamch_ dlamch
#define dlamrg_ dlamrg
#define dlangb_ dlangb
#define dlange_ dlange
#define dlangt_ dlangt
#define dlanhs_ dlanhs
#define dlansb_ dlansb
#define dlansp_ dlansp
#define dlanst_ dlanst
#define dlansy_ dlansy
#define dlantb_ dlantb
#define dlantp_ dlantp
#define dlantr_ dlantr
#define dlanv2_ dlanv2
#define dlapll_ dlapll
#define dlapmt_ dlapmt
#define dlapy2_ dlapy2
#define dlapy3_ dlapy3
#define dlaqgb_ dlaqgb
#define dlaqge_ dlaqge
#define dlaqp2_ dlaqp2
#define dlaqps_ dlaqps
#define dlaqsb_ dlaqsb
#define dlaqsp_ dlaqsp
#define dlaqsy_ dlaqsy
#define dlaqtr_ dlaqtr
#define dlar1v_ dlar1v
#define dlar2v_ dlar2v
#define dlarf_ dlarf
#define dlarfb_ dlarfb
#define dlarfg_ dlarfg
#define dlarft_ dlarft
#define dlarfx_ dlarfx
#define dlargv_ dlargv
#define dlarnv_ dlarnv
#define dlarrb_ dlarrb
#define dlarre_ dlarre
#define dlarrf_ dlarrf
#define dlarrv_ dlarrv
#define dlartg_ dlartg
#define dlartv_ dlartv
#define dlaruv_ dlaruv
#define dlarz_ dlarz
#define dlarzb_ dlarzb
#define dlarzt_ dlarzt
#define dlas2_ dlas2
#define dlascl_ dlascl
#define dlasd0_ dlasd0
#define dlasd1_ dlasd1
#define dlasd2_ dlasd2
#define dlasd3_ dlasd3
#define dlasd4_ dlasd4
#define dlasd5_ dlasd5
#define dlasd6_ dlasd6
#define dlasd7_ dlasd7
#define dlasd8_ dlasd8
#define dlasd9_ dlasd9
#define dlasda_ dlasda
#define dlasdq_ dlasdq
#define dlasdt_ dlasdt
#define dlaset_ dlaset
#define dlasq1_ dlasq1
#define dlasq2_ dlasq2
#define dlasq3_ dlasq3
#define dlasq4_ dlasq4
#define dlasq5_ dlasq5
#define dlasq6_ dlasq6
#define dlasr_ dlasr
#define dlasrt_ dlasrt
#define dlassq_ dlassq
#define dlasv2_ dlasv2
#define dlaswp_ dlaswp
#define dlasy2_ dlasy2
#define dlasyf_ dlasyf
#define dlatbs_ dlatbs
#define dlatdf_ dlatdf
#define dlatps_ dlatps
#define dlatrd_ dlatrd
#define dlatrs_ dlatrs
#define dlatrz_ dlatrz
#define dlatzm_ dlatzm
#define dlauu2_ dlauu2
#define dlauum_ dlauum
#define dopgtr_ dopgtr
#define dopmtr_ dopmtr
#define dorg2l_ dorg2l
#define dorg2r_ dorg2r
#define dorgbr_ dorgbr
#define dorghr_ dorghr
#define dorgl2_ dorgl2
#define dorglq_ dorglq
#define dorgql_ dorgql
#define dorgqr_ dorgqr
#define dorgr2_ dorgr2
#define dorgrq_ dorgrq
#define dorgtr_ dorgtr
#define dorm2l_ dorm2l
#define dorm2r_ dorm2r
#define dormbr_ dormbr
#define dormhr_ dormhr
#define dorml2_ dorml2
#define dormlq_ dormlq
#define dormql_ dormql
#define dormqr_ dormqr
#define dormr2_ dormr2
#define dormr3_ dormr3
#define dormrq_ dormrq
#define dormrz_ dormrz
#define dormtr_ dormtr
#define dpbcon_ dpbcon
#define dpbequ_ dpbequ
#define dpbrfs_ dpbrfs
#define dpbstf_ dpbstf
#define dpbsv_ dpbsv
#define dpbsvx_ dpbsvx
#define dpbtf2_ dpbtf2
#define dpbtrf_ dpbtrf
#define dpbtrs_ dpbtrs
#define dpocon_ dpocon
#define dpoequ_ dpoequ
#define dporfs_ dporfs
#define dposv_ dposv
#define dposvx_ dposvx
#define dpotf2_ dpotf2
#define dpotrf_ dpotrf
#define dpotri_ dpotri
#define dpotrs_ dpotrs
#define dppcon_ dppcon
#define dppequ_ dppequ
#define dpprfs_ dpprfs
#define dppsv_ dppsv
#define dppsvx_ dppsvx
#define dpptrf_ dpptrf
#define dpptri_ dpptri
#define dpptrs_ dpptrs
#define dptcon_ dptcon
#define dpteqr_ dpteqr
#define dptrfs_ dptrfs
#define dptsv_ dptsv
#define dptsvx_ dptsvx
#define dpttrf_ dpttrf
#define dpttrs_ dpttrs
#define dptts2_ dptts2
#define drscl_ drscl
#define dsbev_ dsbev
#define dsbevd_ dsbevd
#define dsbevx_ dsbevx
#define dsbgst_ dsbgst
#define dsbgv_ dsbgv
#define dsbgvd_ dsbgvd
#define dsbgvx_ dsbgvx
#define dsbtrd_ dsbtrd
#define dsecnd_ dsecnd
#define dspcon_ dspcon
#define dspev_ dspev
#define dspevd_ dspevd
#define dspevx_ dspevx
#define dspgst_ dspgst
#define dspgv_ dspgv
#define dspgvd_ dspgvd
#define dspgvx_ dspgvx
#define dsprfs_ dsprfs
#define dspsv_ dspsv
#define dspsvx_ dspsvx
#define dsptrd_ dsptrd
#define dsptrf_ dsptrf
#define dsptri_ dsptri
#define dsptrs_ dsptrs
#define dstebz_ dstebz
#define dstedc_ dstedc
#define dstegr_ dstegr
#define dstein_ dstein
#define dsteqr_ dsteqr
#define dsterf_ dsterf
#define dstev_ dstev
#define dstevd_ dstevd
#define dstevr_ dstevr
#define dstevx_ dstevx
#define dsycon_ dsycon
#define dsyev_ dsyev
#define dsyevd_ dsyevd
#define dsyevr_ dsyevr
#define dsyevx_ dsyevx
#define dsygs2_ dsygs2
#define dsygst_ dsygst
#define dsygv_ dsygv
#define dsygvd_ dsygvd
#define dsygvx_ dsygvx
#define dsyrfs_ dsyrfs
#define dsysv_ dsysv
#define dsysvx_ dsysvx
#define dsytd2_ dsytd2
#define dsytf2_ dsytf2
#define dsytrd_ dsytrd
#define dsytrf_ dsytrf
#define dsytri_ dsytri
#define dsytrs_ dsytrs
#define dtbcon_ dtbcon
#define dtbrfs_ dtbrfs
#define dtbtrs_ dtbtrs
#define dtgevc_ dtgevc
#define dtgex2_ dtgex2
#define dtgexc_ dtgexc
#define dtgsen_ dtgsen
#define dtgsja_ dtgsja
#define dtgsna_ dtgsna
#define dtgsy2_ dtgsy2
#define dtgsyl_ dtgsyl
#define dtpcon_ dtpcon
#define dtprfs_ dtprfs
#define dtptri_ dtptri
#define dtptrs_ dtptrs
#define dtrcon_ dtrcon
#define dtrevc_ dtrevc
#define dtrexc_ dtrexc
#define dtrrfs_ dtrrfs
#define dtrsen_ dtrsen
#define dtrsna_ dtrsna
#define dtrsyl_ dtrsyl
#define dtrti2_ dtrti2
#define dtrtri_ dtrtri
#define dtrtrs_ dtrtrs
#define dtzrqf_ dtzrqf
#define dtzrzf_ dtzrzf
#define dzsum1_ dzsum1
#define f2c_ f2c
#define icmax1_ icmax1
#define ieeeck_ ieeeck
#define ilaenv_ ilaenv
#define izmax1_ izmax1
#define lsame_ lsame
#define lsamen_ lsamen
#define sbdsdc_ sbdsdc
#define sbdsqr_ sbdsqr
#define scsum1_ scsum1
#define sdisna_ sdisna
#define second_ second
#define sgbbrd_ sgbbrd
#define sgbcon_ sgbcon
#define sgbequ_ sgbequ
#define sgbrfs_ sgbrfs
#define sgbsv_ sgbsv
#define sgbsvx_ sgbsvx
#define sgbtf2_ sgbtf2
#define sgbtrf_ sgbtrf
#define sgbtrs_ sgbtrs
#define sgebak_ sgebak
#define sgebal_ sgebal
#define sgebd2_ sgebd2
#define sgebrd_ sgebrd
#define sgecon_ sgecon
#define sgeequ_ sgeequ
#define sgees_ sgees
#define sgeesx_ sgeesx
#define sgeev_ sgeev
#define sgeevx_ sgeevx
#define sgegs_ sgegs
#define sgegv_ sgegv
#define sgehd2_ sgehd2
#define sgehrd_ sgehrd
#define sgelq2_ sgelq2
#define sgelqf_ sgelqf
#define sgels_ sgels
#define sgelsd_ sgelsd
#define sgelss_ sgelss
#define sgelsx_ sgelsx
#define sgelsy_ sgelsy
#define sgeql2_ sgeql2
#define sgeqlf_ sgeqlf
#define sgeqp3_ sgeqp3
#define sgeqpf_ sgeqpf
#define sgeqr2_ sgeqr2
#define sgeqrf_ sgeqrf
#define sgerfs_ sgerfs
#define sgerq2_ sgerq2
#define sgerqf_ sgerqf
#define sgesc2_ sgesc2
#define sgesdd_ sgesdd
#define sgesv_ sgesv
#define sgesvd_ sgesvd
#define sgesvx_ sgesvx
#define sgetc2_ sgetc2
#define sgetf2_ sgetf2
#define sgetrf_ sgetrf
#define sgetri_ sgetri
#define sgetrs_ sgetrs
#define sggbak_ sggbak
#define sggbal_ sggbal
#define sgges_ sgges
#define sggesx_ sggesx
#define sggev_ sggev
#define sggevx_ sggevx
#define sggglm_ sggglm
#define sgghrd_ sgghrd
#define sgglse_ sgglse
#define sggqrf_ sggqrf
#define sggrqf_ sggrqf
#define sggsvd_ sggsvd
#define sggsvp_ sggsvp
#define sgtcon_ sgtcon
#define sgtrfs_ sgtrfs
#define sgtsv_ sgtsv
#define sgtsvx_ sgtsvx
#define sgttrf_ sgttrf
#define sgttrs_ sgttrs
#define sgtts2_ sgtts2
#define shgeqz_ shgeqz
#define shsein_ shsein
#define shseqr_ shseqr
#define slabad_ slabad
#define slabrd_ slabrd
#define slacon_ slacon
#define slacpy_ slacpy
#define sladiv_ sladiv
#define slae2_ slae2
#define slaebz_ slaebz
#define slaed0_ slaed0
#define slaed1_ slaed1
#define slaed2_ slaed2
#define slaed3_ slaed3
#define slaed4_ slaed4
#define slaed5_ slaed5
#define slaed6_ slaed6
#define slaed7_ slaed7
#define slaed8_ slaed8
#define slaed9_ slaed9
#define slaeda_ slaeda
#define slaein_ slaein
#define slaev2_ slaev2
#define slaexc_ slaexc
#define slag2_ slag2
#define slags2_ slags2
#define slagtf_ slagtf
#define slagtm_ slagtm
#define slagts_ slagts
#define slagv2_ slagv2
#define slahqr_ slahqr
#define slahrd_ slahrd
#define slaic1_ slaic1
#define slaln2_ slaln2
#define slals0_ slals0
#define slalsa_ slalsa
#define slalsd_ slalsd
#define slamch_ slamch
#define slamrg_ slamrg
#define slangb_ slangb
#define slange_ slange
#define slangt_ slangt
#define slanhs_ slanhs
#define slansb_ slansb
#define slansp_ slansp
#define slanst_ slanst
#define slansy_ slansy
#define slantb_ slantb
#define slantp_ slantp
#define slantr_ slantr
#define slanv2_ slanv2
#define slapll_ slapll
#define slapmt_ slapmt
#define slapy2_ slapy2
#define slapy3_ slapy3
#define slaqgb_ slaqgb
#define slaqge_ slaqge
#define slaqp2_ slaqp2
#define slaqps_ slaqps
#define slaqsb_ slaqsb
#define slaqsp_ slaqsp
#define slaqsy_ slaqsy
#define slaqtr_ slaqtr
#define slar1v_ slar1v
#define slar2v_ slar2v
#define slarf_ slarf
#define slarfb_ slarfb
#define slarfg_ slarfg
#define slarft_ slarft
#define slarfx_ slarfx
#define slargv_ slargv
#define slarnv_ slarnv
#define slarrb_ slarrb
#define slarre_ slarre
#define slarrf_ slarrf
#define slarrv_ slarrv
#define slartg_ slartg
#define slartv_ slartv
#define slaruv_ slaruv
#define slarz_ slarz
#define slarzb_ slarzb
#define slarzt_ slarzt
#define slas2_ slas2
#define slascl_ slascl
#define slasd0_ slasd0
#define slasd1_ slasd1
#define slasd2_ slasd2
#define slasd3_ slasd3
#define slasd4_ slasd4
#define slasd5_ slasd5
#define slasd6_ slasd6
#define slasd7_ slasd7
#define slasd8_ slasd8
#define slasd9_ slasd9
#define slasda_ slasda
#define slasdq_ slasdq
#define slasdt_ slasdt
#define slaset_ slaset
#define slasq1_ slasq1
#define slasq2_ slasq2
#define slasq3_ slasq3
#define slasq4_ slasq4
#define slasq5_ slasq5
#define slasq6_ slasq6
#define slasr_ slasr
#define slasrt_ slasrt
#define slassq_ slassq
#define slasv2_ slasv2
#define slaswp_ slaswp
#define slasy2_ slasy2
#define slasyf_ slasyf
#define slatbs_ slatbs
#define slatdf_ slatdf
#define slatps_ slatps
#define slatrd_ slatrd
#define slatrs_ slatrs
#define slatrz_ slatrz
#define slatzm_ slatzm
#define slauu2_ slauu2
#define slauum_ slauum
#define sopgtr_ sopgtr
#define sopmtr_ sopmtr
#define sorg2l_ sorg2l
#define sorg2r_ sorg2r
#define sorgbr_ sorgbr
#define sorghr_ sorghr
#define sorgl2_ sorgl2
#define sorglq_ sorglq
#define sorgql_ sorgql
#define sorgqr_ sorgqr
#define sorgr2_ sorgr2
#define sorgrq_ sorgrq
#define sorgtr_ sorgtr
#define sorm2l_ sorm2l
#define sorm2r_ sorm2r
#define sormbr_ sormbr
#define sormhr_ sormhr
#define sorml2_ sorml2
#define sormlq_ sormlq
#define sormql_ sormql
#define sormqr_ sormqr
#define sormr2_ sormr2
#define sormr3_ sormr3
#define sormrq_ sormrq
#define sormrz_ sormrz
#define sormtr_ sormtr
#define spbcon_ spbcon
#define spbequ_ spbequ
#define spbrfs_ spbrfs
#define spbstf_ spbstf
#define spbsv_ spbsv
#define spbsvx_ spbsvx
#define spbtf2_ spbtf2
#define spbtrf_ spbtrf
#define spbtrs_ spbtrs
#define spocon_ spocon
#define spoequ_ spoequ
#define sporfs_ sporfs
#define sposv_ sposv
#define sposvx_ sposvx
#define spotf2_ spotf2
#define spotrf_ spotrf
#define spotri_ spotri
#define spotrs_ spotrs
#define sppcon_ sppcon
#define sppequ_ sppequ
#define spprfs_ spprfs
#define sppsv_ sppsv
#define sppsvx_ sppsvx
#define spptrf_ spptrf
#define spptri_ spptri
#define spptrs_ spptrs
#define sptcon_ sptcon
#define spteqr_ spteqr
#define sptrfs_ sptrfs
#define sptsv_ sptsv
#define sptsvx_ sptsvx
#define spttrf_ spttrf
#define spttrs_ spttrs
#define sptts2_ sptts2
#define srscl_ srscl
#define ssbev_ ssbev
#define ssbevd_ ssbevd
#define ssbevx_ ssbevx
#define ssbgst_ ssbgst
#define ssbgv_ ssbgv
#define ssbgvd_ ssbgvd
#define ssbgvx_ ssbgvx
#define ssbtrd_ ssbtrd
#define sspcon_ sspcon
#define sspev_ sspev
#define sspevd_ sspevd
#define sspevx_ sspevx
#define sspgst_ sspgst
#define sspgv_ sspgv
#define sspgvd_ sspgvd
#define sspgvx_ sspgvx
#define ssprfs_ ssprfs
#define sspsv_ sspsv
#define sspsvx_ sspsvx
#define ssptrd_ ssptrd
#define ssptrf_ ssptrf
#define ssptri_ ssptri
#define ssptrs_ ssptrs
#define sstebz_ sstebz
#define sstedc_ sstedc
#define sstegr_ sstegr
#define sstein_ sstein
#define ssteqr_ ssteqr
#define ssterf_ ssterf
#define sstev_ sstev
#define sstevd_ sstevd
#define sstevr_ sstevr
#define sstevx_ sstevx
#define ssycon_ ssycon
#define ssyev_ ssyev
#define ssyevd_ ssyevd
#define ssyevr_ ssyevr
#define ssyevx_ ssyevx
#define ssygs2_ ssygs2
#define ssygst_ ssygst
#define ssygv_ ssygv
#define ssygvd_ ssygvd
#define ssygvx_ ssygvx
#define ssyrfs_ ssyrfs
#define ssysv_ ssysv
#define ssysvx_ ssysvx
#define ssytd2_ ssytd2
#define ssytf2_ ssytf2
#define ssytrd_ ssytrd
#define ssytrf_ ssytrf
#define ssytri_ ssytri
#define ssytrs_ ssytrs
#define stbcon_ stbcon
#define stbrfs_ stbrfs
#define stbtrs_ stbtrs
#define stgevc_ stgevc
#define stgex2_ stgex2
#define stgexc_ stgexc
#define stgsen_ stgsen
#define stgsja_ stgsja
#define stgsna_ stgsna
#define stgsy2_ stgsy2
#define stgsyl_ stgsyl
#define stpcon_ stpcon
#define stprfs_ stprfs
#define stptri_ stptri
#define stptrs_ stptrs
#define strcon_ strcon
#define strevc_ strevc
#define strexc_ strexc
#define strrfs_ strrfs
#define strsen_ strsen
#define strsna_ strsna
#define strsyl_ strsyl
#define strti2_ strti2
#define strtri_ strtri
#define strtrs_ strtrs
#define stzrqf_ stzrqf
#define stzrzf_ stzrzf
#define xerbla_ xerbla
#define zbdsqr_ zbdsqr
#define zdrot_ zdrot
#define zdrscl_ zdrscl
#define zgbbrd_ zgbbrd
#define zgbcon_ zgbcon
#define zgbequ_ zgbequ
#define zgbrfs_ zgbrfs
#define zgbsv_ zgbsv
#define zgbsvx_ zgbsvx
#define zgbtf2_ zgbtf2
#define zgbtrf_ zgbtrf
#define zgbtrs_ zgbtrs
#define zgebak_ zgebak
#define zgebal_ zgebal
#define zgebd2_ zgebd2
#define zgebrd_ zgebrd
#define zgecon_ zgecon
#define zgeequ_ zgeequ
#define zgees_ zgees
#define zgeesx_ zgeesx
#define zgeev_ zgeev
#define zgeevx_ zgeevx
#define zgegs_ zgegs
#define zgegv_ zgegv
#define zgehd2_ zgehd2
#define zgehrd_ zgehrd
#define zgelq2_ zgelq2
#define zgelqf_ zgelqf
#define zgels_ zgels
#define zgelsd_ zgelsd
#define zgelss_ zgelss
#define zgelsx_ zgelsx
#define zgelsy_ zgelsy
#define zgeql2_ zgeql2
#define zgeqlf_ zgeqlf
#define zgeqp3_ zgeqp3
#define zgeqpf_ zgeqpf
#define zgeqr2_ zgeqr2
#define zgeqrf_ zgeqrf
#define zgerfs_ zgerfs
#define zgerq2_ zgerq2
#define zgerqf_ zgerqf
#define zgesc2_ zgesc2
#define zgesdd_ zgesdd
#define zgesv_ zgesv
#define zgesvd_ zgesvd
#define zgesvx_ zgesvx
#define zgetc2_ zgetc2
#define zgetf2_ zgetf2
#define zgetrf_ zgetrf
#define zgetri_ zgetri
#define zgetrs_ zgetrs
#define zggbak_ zggbak
#define zggbal_ zggbal
#define zgges_ zgges
#define zggesx_ zggesx
#define zggev_ zggev
#define zggevx_ zggevx
#define zggglm_ zggglm
#define zgghrd_ zgghrd
#define zgglse_ zgglse
#define zggqrf_ zggqrf
#define zggrqf_ zggrqf
#define zggsvd_ zggsvd
#define zggsvp_ zggsvp
#define zgtcon_ zgtcon
#define zgtrfs_ zgtrfs
#define zgtsv_ zgtsv
#define zgtsvx_ zgtsvx
#define zgttrf_ zgttrf
#define zgttrs_ zgttrs
#define zgtts2_ zgtts2
#define zhbev_ zhbev
#define zhbevd_ zhbevd
#define zhbevx_ zhbevx
#define zhbgst_ zhbgst
#define zhbgv_ zhbgv
#define zhbgvd_ zhbgvd
#define zhbgvx_ zhbgvx
#define zhbtrd_ zhbtrd
#define zhecon_ zhecon
#define zheev_ zheev
#define zheevd_ zheevd
#define zheevr_ zheevr
#define zheevx_ zheevx
#define zhegs2_ zhegs2
#define zhegst_ zhegst
#define zhegv_ zhegv
#define zhegvd_ zhegvd
#define zhegvx_ zhegvx
#define zherfs_ zherfs
#define zhesv_ zhesv
#define zhesvx_ zhesvx
#define zhetd2_ zhetd2
#define zhetf2_ zhetf2
#define zhetrd_ zhetrd
#define zhetrf_ zhetrf
#define zhetri_ zhetri
#define zhetrs_ zhetrs
#define zhgeqz_ zhgeqz
#define zhpcon_ zhpcon
#define zhpev_ zhpev
#define zhpevd_ zhpevd
#define zhpevx_ zhpevx
#define zhpgst_ zhpgst
#define zhpgv_ zhpgv
#define zhpgvd_ zhpgvd
#define zhpgvx_ zhpgvx
#define zhprfs_ zhprfs
#define zhpsv_ zhpsv
#define zhpsvx_ zhpsvx
#define zhptrd_ zhptrd
#define zhptrf_ zhptrf
#define zhptri_ zhptri
#define zhptrs_ zhptrs
#define zhsein_ zhsein
#define zhseqr_ zhseqr
#define zlabrd_ zlabrd
#define zlacgv_ zlacgv
#define zlacon_ zlacon
#define zlacp2_ zlacp2
#define zlacpy_ zlacpy
#define zlacrm_ zlacrm
#define zlacrt_ zlacrt
#define zladiv_ zladiv
#define zlaed0_ zlaed0
#define zlaed7_ zlaed7
#define zlaed8_ zlaed8
#define zlaein_ zlaein
#define zlaesy_ zlaesy
#define zlaev2_ zlaev2
#define zlags2_ zlags2
#define zlagtm_ zlagtm
#define zlahef_ zlahef
#define zlahqr_ zlahqr
#define zlahrd_ zlahrd
#define zlaic1_ zlaic1
#define zlals0_ zlals0
#define zlalsa_ zlalsa
#define zlalsd_ zlalsd
#define zlangb_ zlangb
#define zlange_ zlange
#define zlangt_ zlangt
#define zlanhb_ zlanhb
#define zlanhe_ zlanhe
#define zlanhp_ zlanhp
#define zlanhs_ zlanhs
#define zlanht_ zlanht
#define zlansb_ zlansb
#define zlansp_ zlansp
#define zlansy_ zlansy
#define zlantb_ zlantb
#define zlantp_ zlantp
#define zlantr_ zlantr
#define zlapll_ zlapll
#define zlapmt_ zlapmt
#define zlaqgb_ zlaqgb
#define zlaqge_ zlaqge
#define zlaqhb_ zlaqhb
#define zlaqhe_ zlaqhe
#define zlaqhp_ zlaqhp
#define zlaqp2_ zlaqp2
#define zlaqps_ zlaqps
#define zlaqsb_ zlaqsb
#define zlaqsp_ zlaqsp
#define zlaqsy_ zlaqsy
#define zlar1v_ zlar1v
#define zlar2v_ zlar2v
#define zlarcm_ zlarcm
#define zlarf_ zlarf
#define zlarfb_ zlarfb
#define zlarfg_ zlarfg
#define zlarft_ zlarft
#define zlarfx_ zlarfx
#define zlargv_ zlargv
#define zlarnv_ zlarnv
#define zlarrv_ zlarrv
#define zlartg_ zlartg
#define zlartv_ zlartv
#define zlarz_ zlarz
#define zlarzb_ zlarzb
#define zlarzt_ zlarzt
#define zlascl_ zlascl
#define zlaset_ zlaset
#define zlasr_ zlasr
#define zlassq_ zlassq
#define zlaswp_ zlaswp
#define zlasyf_ zlasyf
#define zlatbs_ zlatbs
#define zlatdf_ zlatdf
#define zlatps_ zlatps
#define zlatrd_ zlatrd
#define zlatrs_ zlatrs
#define zlatrz_ zlatrz
#define zlatzm_ zlatzm
#define zlauu2_ zlauu2
#define zlauum_ zlauum
#define zpbcon_ zpbcon
#define zpbequ_ zpbequ
#define zpbrfs_ zpbrfs
#define zpbstf_ zpbstf
#define zpbsv_ zpbsv
#define zpbsvx_ zpbsvx
#define zpbtf2_ zpbtf2
#define zpbtrf_ zpbtrf
#define zpbtrs_ zpbtrs
#define zpocon_ zpocon
#define zpoequ_ zpoequ
#define zporfs_ zporfs
#define zposv_ zposv
#define zposvx_ zposvx
#define zpotf2_ zpotf2
#define zpotrf_ zpotrf
#define zpotri_ zpotri
#define zpotrs_ zpotrs
#define zppcon_ zppcon
#define zppequ_ zppequ
#define zpprfs_ zpprfs
#define zppsv_ zppsv
#define zppsvx_ zppsvx
#define zpptrf_ zpptrf
#define zpptri_ zpptri
#define zpptrs_ zpptrs
#define zptcon_ zptcon
#define zpteqr_ zpteqr
#define zptrfs_ zptrfs
#define zptsv_ zptsv
#define zptsvx_ zptsvx
#define zpttrf_ zpttrf
#define zpttrs_ zpttrs
#define zptts2_ zptts2
#define zrot_ zrot
#define zspcon_ zspcon
#define zspmv_ zspmv
#define zspr_ zspr
#define zsprfs_ zsprfs
#define zspsv_ zspsv
#define zspsvx_ zspsvx
#define zsptrf_ zsptrf
#define zsptri_ zsptri
#define zsptrs_ zsptrs
#define zstedc_ zstedc
#define zstegr_ zstegr
#define zstein_ zstein
#define zsteqr_ zsteqr
#define zsycon_ zsycon
#define zsymv_ zsymv
#define zsyr_ zsyr
#define zsyrfs_ zsyrfs
#define zsysv_ zsysv
#define zsysvx_ zsysvx
#define zsytf2_ zsytf2
#define zsytrf_ zsytrf
#define zsytri_ zsytri
#define zsytrs_ zsytrs
#define ztbcon_ ztbcon
#define ztbrfs_ ztbrfs
#define ztbtrs_ ztbtrs
#define ztgevc_ ztgevc
#define ztgex2_ ztgex2
#define ztgexc_ ztgexc
#define ztgsen_ ztgsen
#define ztgsja_ ztgsja
#define ztgsna_ ztgsna
#define ztgsy2_ ztgsy2
#define ztgsyl_ ztgsyl
#define ztpcon_ ztpcon
#define ztprfs_ ztprfs
#define ztptri_ ztptri
#define ztptrs_ ztptrs
#define ztrcon_ ztrcon
#define ztrevc_ ztrevc
#define ztrexc_ ztrexc
#define ztrrfs_ ztrrfs
#define ztrsen_ ztrsen
#define ztrsna_ ztrsna
#define ztrsyl_ ztrsyl
#define ztrti2_ ztrti2
#define ztrtri_ ztrtri
#define ztrtrs_ ztrtrs
#define ztzrqf_ ztzrqf
#define ztzrzf_ ztzrzf
#define zung2l_ zung2l
#define zung2r_ zung2r
#define zungbr_ zungbr
#define zunghr_ zunghr
#define zungl2_ zungl2
#define zunglq_ zunglq
#define zungql_ zungql
#define zungqr_ zungqr
#define zungr2_ zungr2
#define zungrq_ zungrq
#define zungtr_ zungtr
#define zunm2l_ zunm2l
#define zunm2r_ zunm2r
#define zunmbr_ zunmbr
#define zunmhr_ zunmhr
#define zunml2_ zunml2
#define zunmlq_ zunmlq
#define zunmql_ zunmql
#define zunmqr_ zunmqr
#define zunmr2_ zunmr2
#define zunmr3_ zunmr3
#define zunmrq_ zunmrq
#define zunmrz_ zunmrz
#define zunmtr_ zunmtr
#define zupgtr_ zupgtr
#define zupmtr_ zupmtr
#endif /* end of ifdef MATLAB_MEX_FILE */

/*Define the namespace*/
namespace ROPTLIB
{
	// #include <dasum.h>
	// #include <daxpy.h>
	void axpy_(integer *n, complexRopt *ca, complexRopt *cx, integer *incx, complexRopt *cy, integer *incy);
	void axpy_(integer *n, doublerealRopt *da, doublerealRopt *dx, integer *incx, doublerealRopt *dy, integer *incy);
	void axpy_(integer *n, realRopt *sa, realRopt *sx, integer *incx, realRopt *sy, integer *incy);
	void axpy_(integer *n, doublecomplexRopt *za, doublecomplexRopt *zx, integer *incx, doublecomplexRopt *zy, integer *incy);

	// #include <dcabs1.h>
	// #include <dcopy.h>
	void copy_(integer *n, complexRopt *cx, integer *incx, complexRopt *cy, integer *incy);
	void copy_(integer *n, doublerealRopt *cx, integer *incx, doublerealRopt *cy, integer *incy);
	void copy_(integer *n, realRopt *cx, integer *incx, realRopt *cy, integer *incy);
	void copy_(integer *n, doublecomplexRopt *cx, integer *incx, doublecomplexRopt *cy, integer *incy);

	// #include <ddot.h>
	doublerealRopt dot_(integer *n, doublerealRopt *dx, integer *incx, doublerealRopt *dy, integer *incy);
	realRopt dot_(integer *n, realRopt *dx, integer *incx, realRopt *dy, integer *incy);
	void dotu_(complexRopt *ret_val, integer *n, complexRopt *cx, integer *incx, complexRopt *cy, integer *incy);
	void dotu_(doublecomplexRopt *ret_val, integer *n, doublecomplexRopt *cx, integer *incx, doublecomplexRopt *cy, integer *incy);
	void dotc_(complexRopt *ret_val, integer *n, complexRopt *cx, integer *incx, complexRopt *cy, integer *incy);
	void dotc_(doublecomplexRopt *ret_val, integer *n, doublecomplexRopt *cx, integer *incx, doublecomplexRopt *cy, integer *incy);

	// #include <dgbmv.h>
	// #include <dgemm.h>
	void gemm_(char *transa, char *transb, integer *m, integer *n, integer *k, complexRopt *alpha, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, complexRopt *beta, complexRopt *c__, integer *ldc);
	void gemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublerealRopt *alpha, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, doublerealRopt *beta, doublerealRopt *c__, integer *ldc);
	void gemm_(char *transa, char *transb, integer *m, integer *n, integer *k, realRopt *alpha, realRopt *a, integer *lda, realRopt *b, integer *ldb, realRopt *beta, realRopt *c__, integer *ldc);
	void gemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublecomplexRopt *alpha, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *beta, doublecomplexRopt *c__, integer *ldc);

	// #include <dgemv.h>
	void gemv_(char *trans, integer *m, integer *n, complexRopt *alpha, complexRopt *a, integer *lda, complexRopt *x, integer *incx, complexRopt *beta, complexRopt *y, integer *incy);
	void gemv_(char *trans, integer *m, integer *n, doublerealRopt *alpha, doublerealRopt *a, integer *lda, doublerealRopt *x, integer *incx, doublerealRopt *beta, doublerealRopt *y, integer *incy);
	void gemv_(char *trans, integer *m, integer *n, realRopt *alpha, realRopt *a, integer *lda, realRopt *x, integer *incx, realRopt *beta, realRopt *y, integer *incy);
	void gemv_(char *trans, integer *m, integer *n, doublecomplexRopt *alpha, doublecomplexRopt *a, integer *lda, doublecomplexRopt *x, integer *incx, doublecomplexRopt *beta, doublecomplexRopt *y, integer *incy);

	// #include <dger.h>
	void ger_(integer *m, integer *n, doublerealRopt *alpha, doublerealRopt *x, integer *incx, doublerealRopt *y, integer *incy, doublerealRopt *a, integer *lda);
	void ger_(integer *m, integer *n, realRopt *alpha, realRopt *x, integer *incx, realRopt *y, integer *incy, realRopt *a, integer *lda);

	void gerc_(integer *m, integer *n, complexRopt *alpha, complexRopt *x, integer *incx, complexRopt *y, integer *incy, complexRopt *a, integer *lda);
	void gerc_(integer *m, integer *n, doublecomplexRopt *alpha, doublecomplexRopt *x, integer *incx, doublecomplexRopt *y, integer *incy, doublecomplexRopt *a, integer *lda);

	void ger_(integer *m, integer *n, complexRopt *alpha, complexRopt *x, integer *incx, complexRopt *y, integer *incy, complexRopt *a, integer *lda);
	void ger_(integer *m, integer *n, doublecomplexRopt *alpha, doublecomplexRopt *x, integer *incx, doublecomplexRopt *y, integer *incy, doublecomplexRopt *a, integer *lda);

	// #include <dnrm2.h>
	doublerealRopt nrm2_(integer *n, doublerealRopt *x, integer *incx);
	realRopt nrm2_(integer *n, realRopt *x, integer *incx);

	// #include <drot.h>
	// #include <drotg.h>
	// #include <dsbmv.h>
	// #include <dscal.h>
	void scal_(integer *n, complexRopt *ca, complexRopt *cx, integer *incx);
	void scal_(integer *n, doublerealRopt *ca, doublerealRopt *cx, integer *incx);
	void scal_(integer *n, realRopt *ca, realRopt *cx, integer *incx);
	void scal_(integer *n, doublecomplexRopt *ca, doublecomplexRopt *cx, integer *incx);

	// #include <dspmv.h>
	// #include <dspr.h>
	// #include <dspr2.h>
	// #include <dswap.h>
	// #include <dsymm.h>
	// #include <dsymv.h>
	void symv_(char *uplo, integer *n, complexRopt *alpha, complexRopt *a, integer *lda, complexRopt *x, integer *incx, complexRopt *beta, complexRopt *y, integer *incy);
	void symv_(char *uplo, integer *n, doublerealRopt *alpha, doublerealRopt *a, integer *lda, doublerealRopt *x, integer *incx, doublerealRopt *beta, doublerealRopt *y, integer *incy);
	void symv_(char *uplo, integer *n, realRopt *alpha, realRopt *a, integer *lda, realRopt *x, integer *incx, realRopt *beta, realRopt *y, integer *incy);
	void symv_(char *uplo, integer *n, doublecomplexRopt *alpha, doublecomplexRopt *a, integer *lda, doublecomplexRopt *x, integer *incx, doublecomplexRopt *beta, doublecomplexRopt *y, integer *incy);

	// #include <dsyr.h>
	// #include <dsyr2.h>
	// #include <dsyr2k.h>
	// #include <dsyrk.h>
	// #include <dtbmv.h>
	// #include <dtbsv.h>
	// #include <dtpmv.h>
	// #include <dtpsv.h>
	// #include <dtrmm.h>
	// #include <dtrmv.h>
	// #include <dtrsm.h>
	void trsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublerealRopt *alpha, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb);
	void trsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, realRopt *alpha, realRopt *a, integer *lda, realRopt *b, integer *ldb);
	void trsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, doublecomplexRopt *alpha, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb);
	void trsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, complexRopt *alpha, complexRopt *a, integer *lda, complexRopt *b, integer *ldb);
	// #include <dtrsv.h>
	// #include <dzasum.h>
	// #include <dznrm2.h>

	// #include <dbdsdc.h>
	// #include <dbdsqr.h>
	// #include <ddisna.h>
	// #include <dgbbrd.h>
	// #include <dgbcon.h>
	// #include <dgbequ.h>
	// #include <dgbrfs.h>
	// #include <dgbsv.h>
	// #include <dgbsvx.h>
	// #include <dgbtf2.h>
	// #include <dgbtrf.h>
	// #include <dgbtrs.h>
	// #include <dgebak.h>
	// #include <dgebal.h>
	// #include <dgebd2.h>
	// #include <dgebrd.h>
	// #include <dgecon.h>
	// #include <dgeequ.h>
	// #include <dgees.h>
	void gees_(char *jobvs, char *sort, L_fp select, integer *n, complexRopt *a, integer *lda, integer *sdim, complexRopt *w, complexRopt *vs, integer *ldvs, complexRopt *work, integer *lwork, realRopt *rwork, logical *bwork, integer *info);
	void gees_(char *jobvs, char *sort, L_fp select, integer *n, doublerealRopt *a, integer *lda, integer *sdim, doublerealRopt *wr, doublerealRopt *wi, doublerealRopt *vs, integer *ldvs, doublerealRopt *work, integer *lwork, logical *bwork, integer *info);
	void gees_(char *jobvs, char *sort, L_fp select, integer *n, realRopt *a, integer *lda, integer *sdim, realRopt *wr, realRopt *wi, realRopt *vs, integer *ldvs, realRopt *work, integer *lwork, logical *bwork, integer *info);
	void gees_(char *jobvs, char *sort, L_fp select, integer *n, doublecomplexRopt *a, integer *lda, integer *sdim, doublecomplexRopt *w, doublecomplexRopt *vs, integer *ldvs, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, logical *bwork, integer *info);

	// #include <dgeesx.h>
	// #include <dgeev.h>
	// #include <dgeevx.h>
	// #include <dgegs.h>
	// #include <dgegv.h>
	// #include <dgehd2.h>
	// #include <dgehrd.h>
	// #include <dgelq2.h>
	// #include <dgelqf.h>
	// #include <dgels.h>
	// #include <dgelsd.h>
	// #include <dgelss.h>
	// #include <dgelsx.h>
	// #include <dgelsy.h>
	// #include <dgeql2.h>
	// #include <dgeqlf.h>
	// #include <dgeqp3.h>
	void geqp3_(integer *m, integer *n, complexRopt *a, integer *lda, integer *jpvt, complexRopt *tau, complexRopt *work, integer *lwork, realRopt *rwork, integer *info);
	void geqp3_(integer *m, integer *n, doublerealRopt *a, integer *lda, integer *jpvt, doublerealRopt *tau, doublerealRopt *work, integer *lwork, integer *info);
	void geqp3_(integer *m, integer *n, realRopt *a, integer *lda, integer *jpvt, realRopt *tau, realRopt *work, integer *lwork, integer *info);
	void geqp3_(integer *m, integer *n, doublecomplexRopt *a, integer *lda, integer *jpvt, doublecomplexRopt *tau, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *info);
	// #include <dgeqpf.h>
	// #include <dgeqr2.h>
	// #include <dgeqrf.h>
	// #include <dgerfs.h>
	// #include <dgerq2.h>
	// #include <dgerqf.h>
	// #include <dgesc2.h>
	// #include <dgesdd.h>
	void gesdd_(char *jobz, integer *m, integer *n, complexRopt *a, integer *lda, realRopt *s, complexRopt *u, integer *ldu, complexRopt *vt, integer *ldvt, complexRopt *work, integer *lwork, realRopt *rwork, integer *iwork, integer *info);
	void gesdd_(char *jobz, integer *m, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *s, doublerealRopt *u, integer *ldu, doublerealRopt *vt, integer *ldvt, doublerealRopt *work, integer *lwork, integer *iwork, integer *info);
	void gesdd_(char *jobz, integer *m, integer *n, realRopt *a, integer *lda, realRopt *s, realRopt *u, integer *ldu, realRopt *vt, integer *ldvt, realRopt *work, integer *lwork, integer *iwork, integer *info);
	void gesdd_(char *jobz, integer *m, integer *n, doublecomplexRopt *a, integer *lda, doublerealRopt *s, doublecomplexRopt *u, integer *ldu, doublecomplexRopt *vt, integer *ldvt, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *iwork, integer *info);
	// #include <dgesv.h>
	void gesv_(integer *n, integer *nrhs, complexRopt *a, integer *lda, integer *ipiv, complexRopt *b, integer *ldb, integer *info);
	void gesv_(integer *n, integer *nrhs, doublerealRopt *a, integer *lda, integer *ipiv, doublerealRopt *b, integer *ldb, integer *info);
	void gesv_(integer *n, integer *nrhs, realRopt *a, integer *lda, integer *ipiv, realRopt *b, integer *ldb, integer *info);
	void gesv_(integer *n, integer *nrhs, doublecomplexRopt *a, integer *lda, integer *ipiv, doublecomplexRopt *b, integer *ldb, integer *info);

	// #include <dgesvd.h>
	void gesvd_(char *jobu, char *jobvt, integer *m, integer *n, complexRopt *a, integer *lda, realRopt *s, complexRopt *u, integer *ldu, complexRopt *vt, integer *ldvt, complexRopt *work, integer *lwork, realRopt *rwork, integer *info);
	void gesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *s, doublerealRopt *u, integer *ldu, doublerealRopt *vt, integer *ldvt, doublerealRopt *work, integer *lwork, integer *info);
	void gesvd_(char *jobu, char *jobvt, integer *m, integer *n, realRopt *a, integer *lda, realRopt *s, realRopt *u, integer *ldu, realRopt *vt, integer *ldvt, realRopt *work, integer *lwork, integer *info);
	void gesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublecomplexRopt *a, integer *lda, doublerealRopt *s, doublecomplexRopt *u, integer *ldu, doublecomplexRopt *vt, integer *ldvt, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *info);

	// #include <dgesvx.h>
	// #include <dgetc2.h>
	// #include <dgetf2.h>
	// #include <dgetrf.h>
	void getrf_(integer *m, integer *n, complexRopt *a, integer *lda, integer *ipiv, integer *info);
	void getrf_(integer *m, integer *n, doublerealRopt *a, integer *lda, integer *ipiv, integer *info);
	void getrf_(integer *m, integer *n, realRopt *a, integer *lda, integer *ipiv, integer *info);
	void getrf_(integer *m, integer *n, doublecomplexRopt *a, integer *lda, integer *ipiv, integer *info);

	// #include <dgetri.h>
	void getri_(integer *n, complexRopt *a, integer *lda, integer *ipiv, complexRopt *work, integer *lwork, integer *info);
	void getri_(integer *n, doublerealRopt *a, integer *lda, integer *ipiv, doublerealRopt *work, integer *lwork, integer *info);
	void getri_(integer *n, realRopt *a, integer *lda, integer *ipiv, realRopt *work, integer *lwork, integer *info);
	void getri_(integer *n, doublecomplexRopt *a, integer *lda, integer *ipiv, doublecomplexRopt *work, integer *lwork, integer *info);

	// #include <dgetrs.h>
	void getrs_(char *trans, integer *n, integer *nrhs, complexRopt *a, integer *lda, integer *ipiv, complexRopt *b, integer *ldb, integer *info);
	void getrs_(char *trans, integer *n, integer *nrhs, doublerealRopt *a, integer *lda, integer *ipiv, doublerealRopt *b, integer *ldb, integer *info);
	void getrs_(char *trans, integer *n, integer *nrhs, realRopt *a, integer *lda, integer *ipiv, realRopt *b, integer *ldb, integer *info);
	void getrs_(char *trans, integer *n, integer *nrhs, doublecomplexRopt *a, integer *lda, integer *ipiv, doublecomplexRopt *b, integer *ldb, integer *info);
	// #include <dggbak.h>
	// #include <dggbal.h>
	// #include <dgges.h>
	// #include <dggesx.h>
	// #include <dggev.h>
	// #include <dggevx.h>
	// #include <dggglm.h>
	// #include <dgghrd.h>
	// #include <dgglse.h>
	// #include <dggqrf.h>
	// #include <dggrqf.h>
	// #include <dggsvd.h>
	// #include <dggsvp.h>
	// #include <dgtcon.h>
	// #include <dgtrfs.h>
	// #include <dgtsv.h>
	// #include <dgtsvx.h>
	// #include <dgttrf.h>
	// #include <dgttrs.h>
	// #include <dgtts2.h>
	// #include <dhgeqz.h>
	// #include <dhsein.h>
	// #include <dhseqr.h>
	// #include <dlabad.h>
	// #include <dlabrd.h>
	// #include <dlacon.h>
	// #include <dlacpy.h>
	// #include <dladiv.h>
	// #include <dlae2.h>
	// #include <dlaebz.h>
	// #include <dlaed0.h>
	// #include <dlaed1.h>
	// #include <dlaed2.h>
	// #include <dlaed3.h>
	// #include <dlaed4.h>
	// #include <dlaed5.h>
	// #include <dlaed6.h>
	// #include <dlaed7.h>
	// #include <dlaed8.h>
	// #include <dlaed9.h>
	// #include <dlaeda.h>
	// #include <dlaein.h>
	// #include <dlaev2.h>
	// #include <dlaexc.h>
	// #include <dlag2.h>
	// #include <dlags2.h>
	// #include <dlagtf.h>
	// #include <dlagtm.h>
	// #include <dlagts.h>
	// #include <dlagv2.h>
	// #include <dlahqr.h>
	// #include <dlahrd.h>
	// #include <dlaic1.h>
	// #include <dlaln2.h>
	// #include <dlals0.h>
	// #include <dlalsa.h>
	// #include <dlalsd.h>
	// #include <dlamch.h>
	// #include <dlamrg.h>
	// #include <dlangb.h>
	// #include <dlange.h>
	// #include <dlangt.h>
	// #include <dlanhs.h>
	// #include <dlansb.h>
	// #include <dlansp.h>
	// #include <dlanst.h>
	// #include <dlansy.h>
	// #include <dlantb.h>
	// #include <dlantp.h>
	// #include <dlantr.h>
	// #include <dlanv2.h>
	// #include <dlapll.h>
	// #include <dlapmt.h>
	void lapmt_(logical *forwrd, integer *m, integer *n, complexRopt *x, integer *ldx, integer *k);
	void lapmt_(logical *forwrd, integer *m, integer *n, doublerealRopt *x, integer *ldx, integer *k);
	void lapmt_(logical *forwrd, integer *m, integer *n, realRopt *x, integer *ldx, integer *k);
	void lapmt_(logical *forwrd, integer *m, integer *n, doublecomplexRopt *x, integer *ldx, integer *k);

	// #include <dlapy2.h>
	// #include <dlapy3.h>
	// #include <dlaqgb.h>
	// #include <dlaqge.h>
	// #include <dlaqp2.h>
	// #include <dlaqps.h>
	// #include <dlaqsb.h>
	// #include <dlaqsp.h>
	// #include <dlaqsy.h>
	// #include <dlaqtr.h>
	// #include <dlar1v.h>
	// #include <dlar2v.h>
	// #include <dlarf.h>
	// #include <dlarfb.h>
	// #include <dlarfg.h>
	// #include <dlarft.h>
	// #include <dlarfx.h>
	void larfx_(char *side, integer *m, integer *n, complexRopt *v, complexRopt *tau, complexRopt *c__, integer *ldc, complexRopt *work);
	void larfx_(char *side, integer *m, integer *n, doublerealRopt *v, doublerealRopt *tau, doublerealRopt *c__, integer *ldc, doublerealRopt *work);
	void larfx_(char *side, integer *m, integer *n, realRopt *v, realRopt *tau, realRopt *c__, integer *ldc, realRopt *work);
	void larfx_(char *side, integer *m, integer *n, doublecomplexRopt *v, doublecomplexRopt *tau, doublecomplexRopt *c__, integer *ldc, doublecomplexRopt *work);
	// #include <dlargv.h>
	// #include <dlarnv.h>
	// #include <dlarrb.h>
	// #include <dlarre.h>
	// #include <dlarrf.h>
	// #include <dlarrv.h>
	// #include <dlartg.h>
	// #include <dlartv.h>
	// #include <dlaruv.h>
	// #include <dlarz.h>
	// #include <dlarzb.h>
	// #include <dlarzt.h>
	// #include <dlas2.h>
	// #include <dlascl.h>
	// #include <dlasd0.h>
	// #include <dlasd1.h>
	// #include <dlasd2.h>
	// #include <dlasd3.h>
	// #include <dlasd4.h>
	// #include <dlasd5.h>
	// #include <dlasd6.h>
	// #include <dlasd7.h>
	// #include <dlasd8.h>
	// #include <dlasd9.h>
	// #include <dlasda.h>
	// #include <dlasdq.h>
	// #include <dlasdt.h>
	// #include <dlaset.h>
	// #include <dlasq1.h>
	// #include <dlasq2.h>
	// #include <dlasq3.h>
	// #include <dlasq4.h>
	// #include <dlasq5.h>
	// #include <dlasq6.h>
	// #include <dlasr.h>
	// #include <dlasrt.h>
	// #include <dlassq.h>
	// #include <dlasv2.h>
	// #include <dlaswp.h>
	// #include <dlasy2.h>
	// #include <dlasyf.h>
	// #include <dlatbs.h>
	// #include <dlatdf.h>
	// #include <dlatps.h>
	// #include <dlatrd.h>
	// #include <dlatrs.h>
	// #include <dlatrz.h>
	// #include <dlatzm.h>
	// #include <dlauu2.h>
	// #include <dlauum.h>
	// #include <dopgtr.h>
	// #include <dopmtr.h>
	// #include <dorg2l.h>
	// #include <dorg2r.h>
	// #include <dorgbr.h>
	// #include <dorghr.h>
	// #include <dorgl2.h>
	// #include <dorglq.h>
	// #include <dorgql.h>
	// #include <dorgqr.h>
	void orgqr_(integer *m, integer *n, integer *k, doublerealRopt *a, integer *lda, doublerealRopt *tau, doublerealRopt *work, integer *lwork, integer *info);
	void orgqr_(integer *m, integer *n, integer *k, realRopt *a, integer *lda, realRopt *tau, realRopt *work, integer *lwork, integer *info);

	// #include <dorgr2.h>
	// #include <dorgrq.h>
	// #include <dorgtr.h>
	// #include <dorm2l.h>
	// #include <dorm2r.h>
	// #include <dormbr.h>
	// #include <dormhr.h>
	// #include <dorml2.h>
	// #include <dormlq.h>
	// #include <dormql.h>
	// #include <dormqr.h>
	void ormqr_(char *side, char *trans, integer *m, integer *n, integer *k, doublerealRopt *a, integer *lda, doublerealRopt *tau, doublerealRopt *c__, integer *ldc, doublerealRopt *work, integer *lwork, integer *info);
	void ormqr_(char *side, char *trans, integer *m, integer *n, integer *k, realRopt *a, integer *lda, realRopt *tau, realRopt *c__, integer *ldc, realRopt *work, integer *lwork, integer *info);

	// #include <dormr2.h>
	// #include <dormr3.h>
	// #include <dormrq.h>
	// #include <dormrz.h>
	// #include <dormtr.h>
	// #include <dpbcon.h>
	// #include <dpbequ.h>
	// #include <dpbrfs.h>
	// #include <dpbstf.h>
	// #include <dpbsv.h>
	// #include <dpbsvx.h>
	// #include <dpbtf2.h>
	// #include <dpbtrf.h>
	// #include <dpbtrs.h>
	// #include <dpocon.h>
	// #include <dpoequ.h>
	// #include <dporfs.h>
	// #include <dposv.h>
	// #include <dposvx.h>
	// #include <dpotf2.h>
	// #include <dpotrf.h>
	void potrf_(char *uplo, integer *n, complexRopt *a, integer *lda, integer *info);
	void potrf_(char *uplo, integer *n, doublerealRopt *a, integer *lda, integer *info);
	void potrf_(char *uplo, integer *n, realRopt *a, integer *lda, integer *info);
	void potrf_(char *uplo, integer *n, doublecomplexRopt *a, integer *lda, integer *info);

	// #include <dpotri.h>
	// #include <dpotrs.h>
	void potrs_(char *uplo, integer *n, integer *nrhs, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, integer *info);
	void potrs_(char *uplo, integer *n, integer *nrhs, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, integer *info);
	void potrs_(char *uplo, integer *n, integer *nrhs, realRopt *a, integer *lda, realRopt *b, integer *ldb, integer *info);
	void potrs_(char *uplo, integer *n, integer *nrhs, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, integer *info);

	// #include <dppcon.h>
	// #include <dppequ.h>
	// #include <dpprfs.h>
	// #include <dppsv.h>
	// #include <dppsvx.h>
	// #include <dpptrf.h>
	// #include <dpptri.h>
	// #include <dpptrs.h>
	// #include <dptcon.h>
	// #include <dpteqr.h>
	// #include <dptrfs.h>
	// #include <dptsv.h>
	// #include <dptsvx.h>
	// #include <dpttrf.h>
	// #include <dpttrs.h>
	// #include <dptts2.h>
	// #include <drscl.h>
	// #include <dsbev.h>
	// #include <dsbevd.h>
	// #include <dsbevx.h>
	// #include <dsbgst.h>
	// #include <dsbgv.h>
	// #include <dsbgvd.h>
	// #include <dsbgvx.h>
	// #include <dsbtrd.h>
	// #include <dsecnd.h>
	// #include <dspcon.h>
	// #include <dspev.h>
	// #include <dspevd.h>
	// #include <dspevx.h>
	// #include <dspgst.h>
	// #include <dspgv.h>
	// #include <dspgvd.h>
	// #include <dspgvx.h>
	// #include <dsprfs.h>
	// #include <dspsv.h>
	// #include <dspsvx.h>
	// #include <dsptrd.h>
	// #include <dsptrf.h>
	// #include <dsptri.h>
	// #include <dsptrs.h>
	// #include <dstebz.h>
	// #include <dstedc.h>
	// #include <dstegr.h>
	// #include <dstein.h>
	// #include <dsteqr.h>
	// #include <dsterf.h>
	// #include <dstev.h>
	// #include <dstevd.h>
	// #include <dstevr.h>
	// #include <dstevx.h>
	// #include <dsycon.h>
	// #include <dsyev.h>
	void syev_(char *jobz, char *uplo, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *w, doublerealRopt *work, integer *lwork, integer *info);
	void syev_(char *jobz, char *uplo, integer *n, realRopt *a, integer *lda, realRopt *w, realRopt *work, integer *lwork, integer *info);

	// #include <dsyevd.h>
	void syevd_(char *jobz, char *uplo, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *w, doublerealRopt *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
	void syevd_(char *jobz, char *uplo, integer *n, realRopt *a, integer *lda, realRopt *w, realRopt *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
	// #include <dsyevr.h>
	void syevr_(char *jobz, char *range, char *uplo, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *vl, doublerealRopt *vu, integer *il, integer *iu, doublerealRopt *abstol, integer *m, doublerealRopt *w, doublerealRopt *z__, integer *ldz, integer *isuppz, doublerealRopt *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
	void syevr_(char *jobz, char *range, char *uplo, integer *n, realRopt *a, integer *lda, realRopt *vl, realRopt *vu, integer *il, integer *iu, realRopt *abstol, integer *m, realRopt *w, realRopt *z__, integer *ldz, integer *isuppz, realRopt *work, integer *lwork, integer *iwork, integer *liwork, integer *info);
	// #include <dsyevx.h>
	void syevx_(char *jobz, char *range, char *uplo, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *vl, doublerealRopt *vu, integer *il, integer *iu, doublerealRopt *abstol, integer *m, doublerealRopt *w, doublerealRopt *z__, integer *ldz, doublerealRopt *work, integer *lwork, integer *iwork, integer *ifail, integer *info);
	void syevx_(char *jobz, char *range, char *uplo, integer *n, realRopt *a, integer *lda, realRopt *vl, realRopt *vu, integer *il, integer *iu, realRopt *abstol, integer *m, realRopt *w, realRopt *z__, integer *ldz, realRopt *work, integer *lwork, integer *iwork, integer *ifail, integer *info);
	// #include <dsygs2.h>
	// #include <dsygst.h>
	// #include <dsygv.h>
	// #include <dsygvd.h>
	// #include <dsygvx.h>
	// #include <dsyrfs.h>
	// #include <dsysv.h>
	// #include <dsysvx.h>
	// #include <dsytd2.h>
	// #include <dsytf2.h>
	// #include <dsytrd.h>
	// #include <dsytrf.h>
	// #include <dsytri.h>
	// #include <dsytrs.h>
	// #include <dtbcon.h>
	// #include <dtbrfs.h>
	// #include <dtbtrs.h>
	// #include <dtgevc.h>
	// #include <dtgex2.h>
	// #include <dtgexc.h>
	// #include <dtgsen.h>
	// #include <dtgsja.h>
	// #include <dtgsna.h>
	// #include <dtgsy2.h>
	// #include <dtgsyl.h>
	void tgsyl_(char *trans, integer *ijob, integer *m, integer *n, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, complexRopt *c__, integer *ldc, complexRopt *d__, integer *ldd, complexRopt *e, integer *lde, complexRopt *f, integer *ldf, realRopt *scale, realRopt *dif, complexRopt *work, integer *lwork, integer *iwork, integer *info);
	void tgsyl_(char *trans, integer *ijob, integer *m, integer *n, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, doublerealRopt *c__, integer *ldc, doublerealRopt *d__, integer *ldd, doublerealRopt *e, integer *lde, doublerealRopt *f, integer *ldf, doublerealRopt *scale, doublerealRopt *dif, doublerealRopt *work, integer *lwork, integer *iwork, integer *info);
	void tgsyl_(char *trans, integer *ijob, integer *m, integer *n, realRopt *a, integer *lda, realRopt *b, integer *ldb, realRopt *c__, integer *ldc, realRopt *d__, integer *ldd, realRopt *e, integer *lde, realRopt *f, integer *ldf, realRopt *scale, realRopt *dif, realRopt *work, integer *lwork, integer *iwork, integer *info);
	void tgsyl_(char *trans, integer *ijob, integer *m, integer *n, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *c__, integer *ldc, doublecomplexRopt *d__, integer *ldd, doublecomplexRopt *e, integer *lde, doublecomplexRopt *f, integer *ldf, doublerealRopt *scale, doublerealRopt *dif, doublecomplexRopt *work, integer *lwork, integer *iwork, integer *info);
	// #include <dtpcon.h>
	// #include <dtprfs.h>
	// #include <dtptri.h>
	// #include <dtptrs.h>
	// #include <dtrcon.h>
	// #include <dtrevc.h>
	// #include <dtrexc.h>
	// #include <dtrrfs.h>
	// #include <dtrsen.h>
	// #include <dtrsna.h>
	// #include <dtrsyl.h>
	// #include <dtrti2.h>
	// #include <dtrtri.h>
	// #include <dtrtrs.h>
	void trtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, integer *info);
	void trtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublerealRopt *a, integer *lda, doublerealRopt *b, integer *ldb, integer *info);
	void trtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, realRopt *a, integer *lda, realRopt *b, integer *ldb, integer *info);
	void trtrs_(char *uplo, char *trans, char *diag, integer *n, integer *nrhs, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, integer *info);
	// #include <dtzrqf.h>
	// #include <dtzrzf.h>
	// #include <dzsum1.h>

	// zgegs_
	void gegs_(char *jobvsl, char *jobvsr, integer *n, complexRopt *a, integer *lda, complexRopt *b, integer *ldb, complexRopt *alpha, complexRopt *beta, complexRopt *vsl, integer *ldvsl, complexRopt *vsr, integer *ldvsr, complexRopt *work, integer *lwork, realRopt *rwork, integer *info);
	void gegs_(char *jobvsl, char *jobvsr, integer *n, doublecomplexRopt *a, integer *lda, doublecomplexRopt *b, integer *ldb, doublecomplexRopt *alpha, doublecomplexRopt *beta, doublecomplexRopt *vsl, integer *ldvsl, doublecomplexRopt *vsr, integer *ldvsr, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *info);

	// zunmqr_
	void unmqr_(char *side, char *trans, integer *m, integer *n, integer *k, complexRopt *a, integer *lda, complexRopt *tau, complexRopt *c__, integer *ldc, complexRopt *work, integer *lwork, integer *info);
	void unmqr_(char *side, char *trans, integer *m, integer *n, integer *k, doublecomplexRopt *a, integer *lda, doublecomplexRopt *tau, doublecomplexRopt *c__, integer *ldc, doublecomplexRopt *work, integer *lwork, integer *info);

	// zpotri_
	void potri_(char *uplo, integer *n, complexRopt *a, integer *lda, integer *info);
	void potri_(char *uplo, integer *n, doublecomplexRopt *a, integer *lda, integer *info);

	// zheev_
	void heev_(char *jobz, char *uplo, integer *n, complexRopt *A, integer *LDA, realRopt *W, complexRopt *work, integer *lwork, realRopt *rwork, integer *info);
	void heev_(char *jobz, char *uplo, integer *n, doublecomplexRopt *A, integer *LDA, doublerealRopt *W, doublecomplexRopt *work, integer *lwork, doublerealRopt *rwork, integer *info);

	// zungqr_
	void ungqr_(integer *m, integer *n, integer *k, complexRopt *a, integer *lda, complexRopt *tau, complexRopt *work, integer *lwork, integer *info);
	void ungqr_(integer *m, integer *n, integer *k, doublecomplexRopt *a, integer *lda, doublecomplexRopt *tau, doublecomplexRopt *work, integer *lwork, integer *info);

	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, int nrhs, float alpha, blas_sparse_matrix A, const float *b, int ldb, float *c, int ldc);
	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, int nrhs, double alpha, blas_sparse_matrix A, const double *b, int ldb, double *c, int ldc);
	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, int nrhs, const complexRopt *alpha, blas_sparse_matrix A, const complexRopt *b, int ldb, complexRopt *c, int ldc);
	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, int nrhs, const doublecomplexRopt *alpha, blas_sparse_matrix A, const doublecomplexRopt *b, int ldb, doublecomplexRopt *c, int ldc);

	//	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, integer nrhs, float alpha, blas_sparse_matrix A, const float *b, integer ldb, float *c, integer ldc);
	//	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, integer nrhs, double alpha, blas_sparse_matrix A, const double *b, integer ldb, double *c, integer ldc);
	//	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, integer nrhs, const complexRopt *alpha, blas_sparse_matrix A, const complexRopt *b, integer ldb, complexRopt *c, integer ldc);
	//	void BLAS_usmm(enum blas_order_type order, enum blas_trans_type transa, integer nrhs, const doublecomplexRopt *alpha, blas_sparse_matrix A, const doublecomplexRopt *b, integer ldb, doublecomplexRopt *c, integer ldc);

	void BLAS_uscr_insert_entries(blas_sparse_matrix A, integer nz, const float *val, const integer *indx, const integer *jndx);
	void BLAS_uscr_insert_entries(blas_sparse_matrix A, integer nz, const double *val, const integer *indx, const integer *jndx);
	void BLAS_uscr_insert_entries(blas_sparse_matrix A, integer nz, const complexRopt *val, const integer *indx, const integer *jndx);
	void BLAS_uscr_insert_entries(blas_sparse_matrix A, integer nz, const doublecomplexRopt *val, const integer *indx, const integer *jndx);

	//	void BLAS_uscr_insert_entries(blas_sparse_matrix A, int nz, const float *val, const int *indx, const int *jndx);
	//	void BLAS_uscr_insert_entries(blas_sparse_matrix A, int nz, const double *val, const int *indx, const int *jndx);
	//	void BLAS_uscr_insert_entries(blas_sparse_matrix A, int nz, const complexRopt *val, const int *indx, const int *jndx);
	//	void BLAS_uscr_insert_entries(blas_sparse_matrix A, int nz, const doublecomplexRopt *val, const int *indx, const int *jndx);

}; /*end of ROPTLIB namespace*/

#endif
