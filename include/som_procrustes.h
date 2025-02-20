#ifndef SOM_PROCRUSTES_H_
#define SOM_PROCRUSTES_H_

#include "som_utils.h"

class SomProcrustes
{
public:
    SomProcrustes() {}

    /**
     * Standard constructor with problem data inputs (most commonly used)
     */
    SomProcrustes(SomUtils::SomSize somSz, SomUtils::MatD &Tijs, Eigen::MatrixXi &edges)
    {
        sz_ = somSz;
        Tijs_ = Tijs;
        edges_ = edges;
        numEdges_ = Tijs_.cols();
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

    /**
     * Standard constructor with problem data inputs (most commonly used)
     */
    SomProcrustes(const SomUtils::SomSize somSz, const SomUtils::MatD &Tijs, const Eigen::MatrixXi &edges)
    {
        sz_ = somSz;
        Tijs_ = Tijs;
        edges_ = edges;
        numEdges_ = Tijs_.cols();
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

    ~SomProcrustes() {}

    void stepOne(const SomUtils::MatD &T)
    {
        // make M_i, N_i matrices
        // Mmat = zeros(d, num_edges);
        // Nmat = zeros(d, num_edges);
        SomUtils::MatD Mmat = SomUtils::MatD::Zero(sz_.d_, numEdges_);
        SomUtils::MatD Nmat = SomUtils::MatD::Zero(sz_.d_, numEdges_);
        makeMN(T, Mmat, Nmat);

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
            SomUtils::MatD Mi = Mmat.block(0, i, sz_.d_, 1);
            SomUtils::MatD Ni = -Nmat.block(0, i, sz_.d_, 1);
            Eigen::JacobiSVD<SomUtils::MatD> svd(Mi * Ni.transpose(), Eigen::ComputeFullU | Eigen::ComputeFullV);
            SomUtils::MatD U = svd.matrixU();
            SomUtils::MatD V = svd.matrixV();
            SomUtils::MatD tmp = SomUtils::MatD::Identity(sz_.d_, sz_.d_);
            tmp(sz_.d_ - 1, sz_.d_ - 1) = U.determinant() * V.determinant();
            SomUtils::MatD Rii = U * tmp * V.transpose();
            Rcurr_[i] = Rii.transpose();
        }
    }

    void stepTwo(const SomUtils::VecMatD &R)
    {
        // [A,b] = make_A_b(R, T_globalframe, Tijs_vec, edges, params);
        // A=zeros(d*num_edges,d*N);
        SomUtils::MatD A = SomUtils::MatD::Zero(sz_.d_ * numEdges_, sz_.d_ * sz_.n_);
        // b=zeros(d*num_edges, 1);
        SomUtils::MatD b = SomUtils::MatD::Zero(sz_.d_ * numEdges_, 1);
        makeAb(R, A, b);

        // transl_out = A\(-b);
        Tcurr_ = A.colPivHouseholderQr().solve(-b);
    }

    void makeMN(const SomUtils::MatD &T, SomUtils::MatD &Mmat, SomUtils::MatD &Nmat)
    {
        // for edge_id = 1:num_edges
        //     ii = edges(edge_id, 1);
        //     jj = edges(edge_id, 2);
        //     Mmat(:,edge_id) = Tijs_vec(:, edge_id);
        //     Nmat(:,edge_id) = T_globalframe(:,ii) - T_globalframe(:,jj);
        // end
        for (int e = 0; e < numEdges_; ++e)
        {
            int ii = edges_(e, 0) - 1;
            int jj = edges_(e, 1) - 1;
            Mmat.col(e) = Tijs_.col(e);
            Nmat.col(e) = T.col(ii) - T.col(jj);
        }
    }

    void makeAb(const SomUtils::VecMatD &Rgf, SomUtils::MatD &A, SomUtils::MatD &b)
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
        //     b(idxEdges(:,edge_id)) = Tijs_vec(:, edge_id);
        // end
        for (int e = 0; e < numEdges_; ++e)
        {
            int ii = edges_(e, 0) - 1;
            int jj = edges_(e, 1) - 1;
            A.block(sz_.d_ * e, sz_.d_ * ii, sz_.d_, sz_.d_) = Rgf[ii].transpose();
            A.block(sz_.d_ * e, sz_.d_ * jj, sz_.d_, sz_.d_) = -Rgf[ii].transpose();
            b.block(sz_.d_ * e, 0, sz_.d_, 1) = Tijs_.col(e);
        }
    }

    void run()
    {
        /*iterate!*/
        // num_iterations = 0;
        // %COORD DESC - step 3: iterate until convergence
        // while (norm(transf_prev - transf_curr)>= transf_end_thresh && num_iterations<max_icp_iterations)
        // %     rot_prev =  rot_curr;
        // %     transl_prev = transl_curr;
        //     transf_prev = transf_curr;

        //     rot_curr = som_stepone_procrustes(T_globalframe_nois, Tijs_vec, edges, params);

        //     %COORD DESC - step 2
        //     transl_curr = som_steptwo_procrustes(rot_curr, T_globalframe_nois, Tijs_vec, edges, params);
        //     % T = reshape(T, d, []);

        //     transf_curr = make_transf(rot_curr,transl_curr);
        // end
        SomUtils::MatD Tprev = SomUtils::MatD::Zero(sz_.d_, sz_.n_);
        SomUtils::VecMatD Rprev(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
        double normDiff = 1e+10;
        int numIterations = 0;
        do
        {
            ROFL_VAR2("Iteration SOM Procrustes run()", numIterations)
            stepOne(Tcurr_);
            stepTwo(Rcurr_);

            ROFL_VAR3(Tcurr_, Tcurr_.rows(), Tcurr_.cols())
            for (int i = 0; i < sz_.n_; ++i)
            {
                ROFL_VAR1(Rcurr_[i])
            }
            double normDiff = norm(Rprev, Tprev, Rcurr_, Tcurr_);

            Rprev = Rcurr_;
            Tprev = Tcurr_;

            numIterations++;
        } while (normDiff >= transfEndThresh_ && numIterations < maxNumIterations_);

        Rout_ = Rcurr_;
        Tout_ = Tcurr_;
    }

    void vectorizeRT(const SomUtils::VecMatD &R, const SomUtils::MatD &T, SomUtils::MatD &XvecOut) const
    {
        // int fullRotsSz = sz_.p_ * sz_.d_ * sz_.n_;

        // for (int i=0; i<fullRotsSz; ++i) {
        // }

        int n = R.size();

        for (int i = 0; i < n; ++i)
        {
            ROFL_VAR1(R[i])
            ROFL_VAR1(T.col(i))
        }
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

        // TODO: more asserts may be added

        ROFL_ASSERT(fullIdx == XvecOut.rows())
    }

    double norm(const SomUtils::VecMatD &R1, const SomUtils::MatD &T1, const SomUtils::VecMatD &R2, const SomUtils::MatD &T2)
    {
        SomUtils::MatD X1(sz_.d_ * sz_.d_ * sz_.n_ + sz_.d_ * sz_.n_, 1);
        vectorizeRT(R1, T1, X1);
        SomUtils::MatD X2(sz_.d_ * sz_.d_ * sz_.n_ + sz_.d_ * sz_.n_, 1);
        vectorizeRT(R2, T2, X2);
        return (X1 - X2).norm();
    }

    double norm(const SomUtils::VecMatD &R, const SomUtils::MatD &T)
    {
        SomUtils::MatD X(sz_.d_ * sz_.d_ * sz_.n_ + sz_.d_ * sz_.n_, 1);
        vectorizeRT(R, T, X);
        return X.norm();
    }

    void setTcurr(const SomUtils::MatD &Tcurr)
    {
        Tcurr_ = Tcurr;
    }

    SomUtils::MatD getTout()
    {
        return Tout_;
    }

    SomUtils::VecMatD getRout()
    {
        return Rout_;
    }

    /**
     * CLASS MEMBERS
     */

    /**
     * Struct that contains problem size info
     */
    SomUtils::SomSize sz_;

    /**
     * (Estimated) Relative translations between nodes
     * Size: d x e
     */
    SomUtils::MatD Tijs_;

    /**
     * Edges
     * Size: e x 2
     */
    Eigen::MatrixXi edges_;

    /**
     * Basically just edges_.rows()
     * Saved separately for ease
     */
    int numEdges_; // num edges

    /**
     * Full size of problem: pxdxn (rotations) + pxn (translations)
     */
    int fullSz_;

    /**
     * Ground truth information for R
     */
    SomUtils::VecMatD Rgt_;

    /**
     * Ground truth information for T
     */
    SomUtils::MatD Tgt_;

    /**
     * "Global" reference node id (generally 0)
     */
    int src_;

    /**
     * @brief Current cost value (useful at some points e.g. linesearches)
     */
    double costCurr_;

    /**
     * Output of RSOM for R
     */
    SomUtils::VecMatD Rcurr_;

    /**
     * Output of RSOM for T
     */
    SomUtils::MatD Tcurr_;

    /**
     * Output of RSOM for R
     */
    SomUtils::VecMatD Rout_;

    /**
     * Output of RSOM for T
     */
    SomUtils::MatD Tout_;

    /**
     * Threshold for early stopping
     */
    double transfEndThresh_;

    /**
     * Maximum number of iterations
     */
    int maxNumIterations_;
};

#endif /*SOM_PROCRUSTES_H_*/