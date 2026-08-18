#include <ssom_procrustes.h>

SsomProcrustes::SsomProcrustes() {}

SsomProcrustes::SsomProcrustes(SomUtils::SomSize somSz, SomUtils::MatD &tijs, Eigen::MatrixXi &edges)
{
    sz_ = somSz;
    tijs_ = tijs;
    edges_ = edges;
    numEdges_ = tijs_.cols();
    fullSz_ = sz_.d_ * sz_.p_ * sz_.n_ + sz_.p_ * sz_.n_;

    Rgt_.resize(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
    Tgt_ = SomUtils::MatD::Zero(sz_.d_, sz_.n_);

    Rcurr_.resize(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
    Tcurr_ = SomUtils::MatD::Zero(sz_.d_, sz_.n_);

    Rout_.resize(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
    Tout_ = SomUtils::MatD::Zero(sz_.d_, sz_.n_);

    src_ = 0; // TODO: src_ VS src (for sure in globalize, maybe also in other places)

    costCurr_ = 1e+10;

    transfEndThresh_ = 1e-3;
    maxNumIterations_ = 10;
}

SsomProcrustes::SsomProcrustes(const SomUtils::SomSize somSz, const SomUtils::MatD &tijs, const Eigen::MatrixXi &edges)
{
    sz_ = somSz;
    tijs_ = tijs;
    edges_ = edges;
    numEdges_ = tijs_.cols();
    fullSz_ = sz_.d_ * sz_.p_ * sz_.n_ + sz_.p_ * sz_.n_;

    Rgt_.resize(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
    Tgt_ = SomUtils::MatD::Zero(sz_.d_, sz_.n_);

    Rcurr_.resize(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
    Tcurr_ = SomUtils::MatD::Zero(sz_.d_, sz_.n_);

    Rout_.resize(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
    Tout_ = SomUtils::MatD::Zero(sz_.d_, sz_.n_);

    src_ = 0; // TODO: src_ VS src (for sure in globalize, maybe also in other places)

    costCurr_ = 1e+10;

    transfEndThresh_ = 1e-3;
    maxNumIterations_ = 10;
}

SsomProcrustes::~SsomProcrustes() {}

void SsomProcrustes::getMiFromMmat(const SomUtils::MatD &Mmat, SomUtils::MatD &Mi, int idx) const
{
    // for edge_id = 1:size(edges, 1)
    //     if edges(edge_id, 1) == idx
    //         jj = edges(edge_id, 2);
    //         m_i(:, jj) = Mmat(:, edge_id);
    //     end
    // end
    for (int edgeIdx = 0; edgeIdx < numEdges_; ++edgeIdx)
    {
        if (edges_(edgeIdx, 0) - 1 == idx)
        {
            int jj = edges_(edgeIdx, 1) - 1;
            Mi.col(jj) = Mmat.col(edgeIdx);
        }
    }
}

void SsomProcrustes::stepOne(const SomUtils::MatD &T, const SomUtils::MatD &Lambdas)
{
    // make M_i, N_i matrices
    // Mmat = zeros(d, num_edges);
    // Nmat = zeros(d, num_edges);
    SomUtils::MatD Mmat = SomUtils::MatD::Zero(sz_.d_, numEdges_);
    SomUtils::MatD Nmat = SomUtils::MatD::Zero(sz_.d_, numEdges_);
    SomUtils::MatD TijsScaled = SomUtils::MatD::Zero(sz_.d_, numEdges_);
    makeTijsScaled(tijs_, Lambdas, TijsScaled);
    makeMN(T, Mmat, Nmat, TijsScaled);

    // for ii = 1:N
    //     M_i = get_Mi_from_Mmat(Mmat, edges, ii, d, N);
    //     N_i = -get_Ni_from_Nmat(Nmat, edges, ii, d, N);
    //     %compute R_ii
    //     [U,~,V] = svd(M_i*(N_i)'); % ~ would be S
    //     %when d=3 -> quasi_diag = diag(1,1,det(U*V')
    //     R_ii = U * diag([1,1,det(U*V')]) * V';

    //     %fill retval matrix with R_ii at corresponding index
    //     R(:,:,ii) = R_ii';

    for (int i = 0; i < sz_.n_; ++i)
    {
        SomUtils::MatD Mi = SomUtils::MatD::Zero(sz_.d_, sz_.n_);
        SomUtils::MatD Ni = SomUtils::MatD::Zero(sz_.d_, sz_.n_);
        getMiFromMmat(Mmat, Mi, i);
        getMiFromMmat(Nmat, Ni, i);
        Ni *= -1;
        // ROFL_VAR2(i, Mi)
        // ROFL_VAR2(i, Ni)
        Eigen::JacobiSVD<SomUtils::MatD> svd(Mi * Ni.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
        SomUtils::MatD U = svd.matrixU();
        SomUtils::MatD V = svd.matrixV();
        SomUtils::MatD tmp = SomUtils::MatD::Identity(sz_.d_, sz_.d_);
        tmp(sz_.d_ - 1, sz_.d_ - 1) = U.determinant() * V.determinant();
        SomUtils::MatD Rii = U * tmp * V.transpose();
        Rcurr_[i] = Rii.transpose();
    }
}

void SsomProcrustes::stepTwo(const SomUtils::VecMatD &R, const SomUtils::MatD &Lambdas)
{
    // [A,b] = make_A_b(R, T_globalframe, tijs_vec, edges, params);
    // A=zeros(d*num_edges,d*N);
    SomUtils::MatD A = SomUtils::MatD::Zero(sz_.d_ * numEdges_, sz_.d_ * sz_.n_);
    // b=zeros(d*num_edges, 1);
    SomUtils::MatD b = SomUtils::MatD::Zero(sz_.d_ * numEdges_, 1);
    SomUtils::MatD TijsScaled = SomUtils::MatD::Zero(sz_.d_, numEdges_);
    makeTijsScaled(tijs_, Lambdas, TijsScaled);
    makeAb(R, A, b, TijsScaled);

    // transl_out = A\(-b);

    // ROFL_VAR2(A, b)
    SomUtils::MatD TcurrVec = A.colPivHouseholderQr().solve(-b);
    // ROFL_VAR3(A, b, TcurrVec)
    // ROFL_VAR2("Residuals:", A * TcurrVec + b)
    Tcurr_ = TcurrVec.reshaped<Eigen::ColMajor>(sz_.d_, sz_.n_); // HP) this resize does the intended job
    // ROFL_VAR1(Tcurr_)
    // ROFL_VAR1(Tgt_)
}

void SsomProcrustes::stepThree(const SomUtils::VecMatD &R,
                               const SomUtils::MatD &T)
{
    // num_edges = size(edges, 1);
    // lambdas_out = zeros(num_edges, 1);
    int numEdges = edges_.rows();
    // for ee = 1 : num_edges

    //      tij = tijs( :, ee);

    //      ii = edges(ee, 1);
    //      jj = edges(ee, 2);
    //      Ri = R( :, :, ii);
    //      % Rj = R( :, :, jj);
    //      Ti = T( :, ii);
    //      Tj = T( :, jj);

    //      a = tij' * tij; b1 = tij ' * Ri' * Ti;
    //      b2 = -tij ' * Ri' * Tj;
    //      b3 = Ti' * Ri * tij; b4 = -Tj' * Ri * tij; b = b1 + b2 + b3 + b4;
    //      c1 = Ti' * Ti; c2 = -Ti' * Tj; c3 = -Tj' * Ti; % OBS: should be equal to c2 c4 = Tj' * Tj; c = c1 + c2 + c3 + c4;

    //      x = solve_quadratic_lambdas(a, b, c);

    //      lambdas_out(ee) = x;

    for (int e = 0; e < numEdges; ++e)
    {
        auto tij = tijs_.col(e);
        int ii = edges_(e, 0) - 1;
        int jj = edges_(e, 1) - 1;
        auto Ri = R[ii];
        auto Ti = T.col(ii);
        auto Tj = T.col(jj);

        double a = tij.transpose() * tij;
        double b1 = tij.transpose() * Ri.transpose() * Ti;
        double b2 = -tij.transpose() * Ri.transpose() * Tj;
        double b3 = Ti.transpose() * Ri * tij;
        double b4 = -Tj.transpose() * Ri * tij;
        double b = b1 + b2 + b3 + b4;

        double c1 = Ti.transpose() * Ti;
        double c2 = -Ti.transpose() * Tj;
        double c3 = -Tj.transpose() * Ti;
        double c4 = Tj.transpose() * Tj;
        double c = c1 + c2 + c3 + c4;

        LambdasCurr_(e) = solveQuadraticLambdas(a, b, c);
    }
}

double SsomProcrustes::solveQuadraticLambdas(const double a, const double b, const double c)
{
    // D = b.^2 - 4.*a.*c;                    % discriminant

    double D = b * b - 4.0 * a * c;

    // if D < 0
    //     x = 1;
    //     return;
    if (D < 0)
    {
        return 1.0;
    }
    else
    {
        double sqrtD = sqrt(D);
        double x1 = (-b + sqrtD) / (2.0 * a);
        double x2 = (-b - sqrtD) / (2.0 * a);
        if (x1 < 1.0 && x2 >= 1.0)
        {
            return x2;
        }
        else if (x2 < 1.0 && x1 >= 1.0)
        {
            return x1;
        }
        else if (x1 >= 1.0 && x2 >= 1.0)
        {
            return std::max(x1, x2);
        }
        else
        {
            // which case are we in?
            return 1.0;
        }
    }
    // else
    //     sqrtD = sqrt(D);
    //     x1 = (-b + sqrtD) / (2.*a);
    //     x2 = (-b - sqrtD) / (2.*a);
    //     if x1 < 1 && x2 >=1
    //         x = x2;
    //     elseif  x2 < 1 && x1 >=1
    //         x = x1;
    //     elseif x1 >= 1 && x2 >= 1
    //         x = max(x1, x2);
    //     else
    //         %which case are we in?
    //         x = 1;
}

void SsomProcrustes::makeMN(const SomUtils::MatD &T, SomUtils::MatD &Mmat, SomUtils::MatD &Nmat, SomUtils::MatD &TijsScaled) const
{
    // for edge_id = 1:num_edges
    //     ii = edges(edge_id, 1);
    //     jj = edges(edge_id, 2);
    //     Mmat(:,edge_id) = tijs_vec(:, edge_id);
    //     Nmat(:,edge_id) = T_globalframe(:,ii) - T_globalframe(:,jj);
    // end
    for (int e = 0; e < numEdges_; ++e)
    {
        int ii = edges_(e, 0) - 1;
        int jj = edges_(e, 1) - 1;
        Mmat.col(e) = TijsScaled.col(e);
        Nmat.col(e) = T.col(ii) - T.col(jj);
    }

    // ROFL_VAR2(Mmat, Nmat)
}

void SsomProcrustes::makeTijsScaled(const SomUtils::MatD &Tijs, const SomUtils::MatD &Lambdas, SomUtils::MatD &TijsScaled) const
{
    ROFL_ASSERT_VAR5(Tijs.rows() == TijsScaled.rows() && Tijs.cols() == TijsScaled.cols() && Lambdas.rows() == TijsScaled.cols(), Tijs.rows(), TijsScaled.rows(), Tijs.cols(), TijsScaled.cols(), Lambdas.rows());
    TijsScaled = Tijs;
    for (int e = 0; e < numEdges_; ++e)
    {
        int i = edges_(e, 0) - 1; // !! -1
        int j = edges_(e, 1) - 1; // !! -1

        auto tij = Tijs.col(e);
        double lambdaIJ = Lambdas(e, 0);

        TijsScaled.col(e) = tij * lambdaIJ;
    }
}

void SsomProcrustes::makeAb(const SomUtils::VecMatD &Rgf, SomUtils::MatD &A, SomUtils::MatD &b, SomUtils::MatD &TijsScaled) const
{
    // idxEdges=reshape(1:d*num_edges,d,num_edges);
    // idxNodes=reshape(1:d*N,d,N);

    // for edge_id = 1:num_edges
    //     ii = edges(edge_id, 1);
    //     jj = edges(edge_id, 2);
    //     %             A(col_id*d-(d-1):col_id*d, (ii*d)-(d-1):ii*d) = R(:,:,ii)';
    //     %             A(col_id*d-(d-1):col_id*d, (jj*d)-(d-1):jj*d) = -R(:,:,ii)';
    //     A(idxEdges(:,edge_id),idxNodes(:,ii)) = R_gf(:,:,ii)';
    //     A(idxEdges(:,edge_id),idxNodes(:,jj)) = -R_gf(:,:,ii)';
    //     b(idxEdges(:,edge_id)) = tijs_vec(:, edge_id);
    // end
    for (int e = 0; e < numEdges_; ++e)
    {
        int ii = edges_(e, 0) - 1;
        int jj = edges_(e, 1) - 1;
        A.block(sz_.d_ * e, sz_.d_ * ii, sz_.d_, sz_.d_) = Rgf[ii].transpose();  // HP) copilot made the indices correct
        A.block(sz_.d_ * e, sz_.d_ * jj, sz_.d_, sz_.d_) = -Rgf[ii].transpose(); // HP) copilot made the indices correct
        b.block(sz_.d_ * e, 0, sz_.d_, 1) = TijsScaled.col(e);                        // HP) copilot made the indices correct
    }

    // ROFL_VAR2(A, b.transpose())
}

void SsomProcrustes::run()
{
    /*iterate!*/
    // num_iterations = 0;
    // %COORD DESC - step 3: iterate until convergence
    // while (norm(transf_prev - transf_curr)>= transf_end_thresh && num_iterations<max_icp_iterations)
    // %     rot_prev =  rot_curr;
    // %     transl_prev = transl_curr;
    //     transf_prev = transf_curr;

    //     rot_curr = som_stepone_procrustes(T_globalframe_nois, tijs_vec, edges, params);

    //     %COORD DESC - step 2
    //     transl_curr = som_steptwo_procrustes(rot_curr, T_globalframe_nois, tijs_vec, edges, params);
    //     % T = reshape(T, d, []);

    //     transf_curr = make_transf(rot_curr,transl_curr);
    // end
    SomUtils::MatD Tprev = SomUtils::MatD::Zero(sz_.d_, sz_.n_);
    SomUtils::VecMatD Rprev(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
    SomUtils::MatD LambdasPrev = SomUtils::MatD::Zero(numEdges_, 1);
    double normDiff = 1e+10;
    int numIter = 0;
    do
    {
        ROFL_VAR2("Iteration SOM Procrustes run()", numIter)
        stepOne(Tcurr_, LambdasCurr_);
        stepTwo(Rcurr_, LambdasCurr_);
        stepThree(Rcurr_, Tcurr_);

        int src = 0;

        auto Rglobal = Rcurr_[src] * Rgt_[src].transpose();
        auto Tglobal = Rglobal * Tcurr_.col(src) - Tgt_.col(src);
        double lambdaFactor = LambdasGt_(src) / LambdasCurr_(src); // should be the same for all edges
        auto LambdasGlobal = lambdaFactor * LambdasCurr_;
        // ROFL_VAR4(Rglobal, (Rglobal * Tsedn.col(src)).transpose(), Tsedn.col(src).transpose(), Tgt_.col(src).transpose());
        // ROFL_VAR2(Tglobal.transpose(), Tsedn.col(src).transpose());
        // ROFL_VAR3(Rglobal.transpose(), Tcurr_, Tglobal)
        SomUtils::MatD TglobalRepmat(SomUtils::MatD::Zero(sz_.d_, sz_.n_));
        for (int i = 0; i < sz_.n_; ++i)
        {
            TglobalRepmat.col(i) = Tglobal; // TODO: maybe use some other adv init
        }
        auto TrecoveredGlobal = Rglobal.transpose() * Tcurr_ - TglobalRepmat;
        Tcurr_ = TrecoveredGlobal;
        ROFL_VAR3(Tcurr_, Tcurr_.rows(), Tcurr_.cols())

        for (int i = 0; i < sz_.n_; ++i)
        {
            ROFL_VAR1(Rcurr_[i])
        }

        ROFL_VAR1(LambdasCurr_.transpose())

        normDiff = norm(Rprev, Tprev, LambdasPrev, Rcurr_, Tcurr_, LambdasCurr_);

        ROFL_VAR1(normDiff)
        ROFL_VAR1("----------------------------------------------------------\n")

        Rprev = Rcurr_;
        Tprev = Tcurr_;
        LambdasPrev = LambdasCurr_;

        numIter++;
    } while (normDiff >= transfEndThresh_ && numIter < maxNumIterations_);

    Rout_ = Rcurr_;
    Tout_ = Tcurr_;
    LambdasOut_ = LambdasCurr_;

    costCurr_ = costEigen(Rcurr_, Tcurr_, LambdasCurr_);
}

double SsomProcrustes::getCost() const
{
    return costCurr_;
}

double SsomProcrustes::costEigen(const SomUtils::VecMatD &Reigen, const SomUtils::MatD &Teigen, const SomUtils::MatD &LambdasEigen) const
{
    double cost = 0.0f;
    for (int e = 0; e < numEdges_; ++e)
    {
        int i = edges_(e, 0) - 1; // !! -1
        int j = edges_(e, 1) - 1; // !! -1

        auto Ri = Reigen[i];
        auto Ti = Teigen.col(i);
        auto Tj = Teigen.col(j);

        SomUtils::VecD tij(SomUtils::VecD::Zero(sz_.d_));
        tij = tijs_.col(e);

        // ROFL_VAR3(i, j, e);
        // ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

        auto lambdaE = LambdasEigen(e, 0);

        double costEsq = (Ri * tij * lambdaE - Tj + Ti).squaredNorm(); // TODO: use squaredNorm() here directly
        // ROFL_VAR1(costE);

        cost += costEsq;
    }
    return cost;
}

void SsomProcrustes::setTstart(const SomUtils::MatD &Tstart)
{
    Tcurr_ = Tstart;
}

void SsomProcrustes::setLambdasStart(const SomUtils::MatD &LambdasStart)
{
    LambdasCurr_ = LambdasStart;
}

void SsomProcrustes::setRgt(const SomUtils::VecMatD &Rgt)
{
    Rgt_ = Rgt;
}

void SsomProcrustes::setTgt(const SomUtils::MatD &Tgt)
{
    Tgt_ = Tgt;
}

void SsomProcrustes::setLambdasGt(const SomUtils::MatD &LambdasGt)
{
    LambdasGt_ = LambdasGt;
}

SomUtils::MatD SsomProcrustes::getTout() const
{
    return Tout_;
}

SomUtils::VecMatD SsomProcrustes::getRout() const
{
    return Rout_;
}

SomUtils::MatD SsomProcrustes::getLambdasOut() const
{
    return LambdasOut_;
}

void SsomProcrustes::vectorizeRTLambdas(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas, SomUtils::MatD &XvecOut) const
{
    // int fullRotsSz = sz_.p_ * sz_.d_ * sz_.n_;

    // for (int i=0; i<fullRotsSz; ++i) {
    // }

    int n = R.size();

    // for (int i = 0; i < n; ++i)
    // {
    //     ROFL_VAR1(R[i])
    //     ROFL_VAR1(T.col(i))
    // }
    ROFL_ASSERT(T.cols() == n)

    int fullIdx = 0;
    for (int i = 0; i < n; ++i)
    {
        int ric = R[i].cols();
        int rir = R[i].rows();
        for (int j = 0; j < ric; ++j)
        {
            for (int k = 0; k < rir; ++k)
            {
                XvecOut(fullIdx, 0) = R[i](k, j);
                fullIdx++;
                // ROFL_VAR4(i, j, k, fullIdx);
            }
        }
    }

    for (int i = 0; i < T.cols(); ++i)
    {
        for (int j = 0; j < T.rows(); ++j)
        {
            XvecOut(fullIdx, 0) = T(j, i);
            fullIdx++;
            // ROFL_VAR4(i, j, k, fullIdx);
        }
    }

    for (int i = 0; i < Lambdas.rows(); ++i)
    {
        XvecOut(fullIdx, 0) = Lambdas(i, 0);
        fullIdx++;
    }

    // TODO: more asserts may be added

    ROFL_ASSERT(fullIdx == XvecOut.rows())
}

double SsomProcrustes::norm(const SomUtils::VecMatD &R1, const SomUtils::MatD &T1, const SomUtils::MatD &Lambdas1, const SomUtils::VecMatD &R2, const SomUtils::MatD &T2, const SomUtils::MatD &Lambdas2) const
{
    SomUtils::MatD X1(sz_.d_ * sz_.d_ * sz_.n_ + sz_.d_ * sz_.n_ + Lambdas1.rows(), 1);
    vectorizeRTLambdas(R1, T1, Lambdas1, X1);
    SomUtils::MatD X2(sz_.d_ * sz_.d_ * sz_.n_ + sz_.d_ * sz_.n_ + Lambdas2.rows(), 1);
    vectorizeRTLambdas(R2, T2, Lambdas2, X2);
    return (X1 - X2).norm();
}

double SsomProcrustes::norm(const SomUtils::VecMatD &R, const SomUtils::MatD &T, const SomUtils::MatD &Lambdas) const
{
    SomUtils::MatD X(sz_.d_ * sz_.d_ * sz_.n_ + sz_.d_ * sz_.n_ + Lambdas.rows(), 1);
    vectorizeRTLambdas(R, T, Lambdas, X);
    return X.norm();
}

void SsomProcrustes::setTcurr(const SomUtils::MatD &Tcurr)
{
    Tcurr_ = Tcurr;
}
