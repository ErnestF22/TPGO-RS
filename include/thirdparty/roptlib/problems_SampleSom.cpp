#include "problems_SampleSom.h"

namespace ROPTLIB
{

    SampleSomProblem::SampleSomProblem(std::vector<double> &rshSrc, std::vector<double> &rshDst, int arsLmax)
    {
        rshSrc_ = rshSrc;
        rshDst_ = rshDst;
        lmax_ = arsLmax;
    };

    SampleSomProblem::~SampleSomProblem(){};

    realdp SampleSomProblem::f(const Variable &x) const
    {
        // Check whether rshSrc, rshDst are actually well ok

        // MATLAB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // eul_angles = rotm2eul(x, 'ZYZ');
        // alpha = eul_angles(1);
        // beta = eul_angles(2);
        // gamma = eul_angles(3);
        // c = - rsh_corr_grad(rsh_src, rsh_dst, alpha, beta, gamma, lmax);

        // ars::RshRotation arsRshRot;
        // arsRshRot.create(lmax_, xEigen);
        // double corr = arsRshRot.correlation(rshSrc_, rshDst_);

        Eigen::Matrix3d xEigen;
        ROFL_VAR1(x);
        ROFL_VAR1(xEigen);
        RoptToEig(x, xEigen);

        const Eigen::Vector3d ea = xEigen.eulerAngles(2, 1, 2);
        ROFL_VAR1(ea.transpose());
        double corr;
        Eigen::Vector3d grad;    // unused here
        Eigen::Matrix3d Hessian; // unused (and unimplemented)
        rsh_corr_grad(ea, corr, grad, Hessian);

        ROFL_VAR1(corr);
        ROFL_ASSERT(!std::isnan(corr));
        // system("pause");
        return -corr; // checked -> the - here should be OK
    };

    Vector &SampleSomProblem::EucGrad(const Variable &x, Vector *result) const
    {
        // eul_angles = rotm2eul(x, 'ZYZ');
        // alpha = eul_angles(1);
        // beta = eul_angles(2);
        // gamma = eul_angles(3);
        Eigen::Matrix3d xEigen;
        RoptToEig(x, xEigen);
        Eigen::Vector3d xEulAngles = xEigen.eulerAngles(2, 1, 2);

        // [~, grad, ~, ~] = rsh_corr_grad(rsh_src, rsh_dst, alpha, beta, gamma, lmax);
        double corr; // unused here
        Eigen::Vector3d grad;
        Eigen::Matrix3d Hessian; // unused (and unimplemented)
        // std::vector<Eigen::MatrixXd> Mf(lmax_+1); // unused
        // Mf.reserve(lmax_+1);
        // rsh_corr_grad(xEulAngles, corr, grad, Hessian, Mf);
        rsh_corr_grad(xEulAngles, corr, grad, Hessian);

        // J = euler_angles_jacobian_zyz(alpha,beta,gamma);
        // pinvJ = pinv(J);
        Eigen::MatrixXd J(9, 3);
        J.setZero();
        // eulerAnglesJacobianZYZ(lmax_, xEulAngles, J);
        Eigen::MatrixXd pinvJ = J.completeOrthogonalDecomposition().pseudoInverse();

        // Gvec = g * pinvJ;
        // G = - reshape(Gvec, 3,3);
        // // ROFL_VAR3(grad.rows(), grad.cols(), grad.transpose());
        // // ROFL_VAR3(J.rows(), J.cols(), J);
        // // ROFL_VAR3(pinvJ.rows(), pinvJ.cols(), pinvJ);
        Eigen::RowVectorXd Gvec(9);
        Gvec = -grad.transpose() * pinvJ; // grad.transpose() is NOT NEEDED in Matlab as rsh_corr_grad already outputs a row vector
        // %     G = zeros(3,3);
        Eigen::Matrix3d G(Gvec.data());
        G.resize(3, 3);
        G.transposeInPlace();

        // back to Roptlib
        Vector Groptlib = Rotations(3).RandominManifold();
        Groptlib.ObtainWriteEntireData();
        realdp *GroptlibWriteArray = Groptlib.ObtainWriteEntireData(); //!! assignment in col-major order
        // Groptlib
        GroptlibWriteArray[0] = G(0, 0);
        GroptlibWriteArray[1] = G(0, 1);
        GroptlibWriteArray[2] = G(0, 2);
        GroptlibWriteArray[3] = G(1, 0);
        GroptlibWriteArray[4] = G(1, 1);
        GroptlibWriteArray[5] = G(1, 2);
        GroptlibWriteArray[6] = G(2, 0);
        GroptlibWriteArray[7] = G(2, 1);
        GroptlibWriteArray[8] = G(2, 2);
        // Groptlib.Print("Groptlib"); //This was causing a loop
        ROFL_VAR1(xEigen);
        ROFL_VAR1(xEulAngles.transpose());
        ROFL_VAR1(G);

        result->NewMemoryOnWrite();
        result->SetToZeros();
        *result = Groptlib;
        return *result;
    };

    void SampleSomProblem::RoptToEig(Vector x, Eigen::Matrix3d &xEigen) const
    {
        Vector xT = x.GetTranspose(); // Eigen ADV init is row-major!!
        // Eigen::Matrix3d xEigen;

        const realdp *xArr = xT.ObtainWriteEntireData();
        xEigen << xArr[0], xArr[1], xArr[2], xArr[3], xArr[4], xArr[5], xArr[6], xArr[7], xArr[8];
    }

    // Vector &SampleSomProblem::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
    // {
    //     //TODO: implement
    //     return *result;
    // };

    void SampleSomProblem::rsh_corr_grad(const Eigen::Vector3d &ea,
                                         double &corr,
                                         Eigen::Vector3d &grad,
                                         Eigen::Matrix3d &Hessian,
                                         std::vector<Eigen::MatrixXd> &Mf) const
    {
        ea.trace();
        // ...
        // function [corr, grad, Hessian, Mf] = rsh_corr_grad(rsh_src,rsh_dst,alpha,beta,gamma,lmax)
        // %
        // % Note: the function assumes that ZYZ Euler angles are used:
        // %  Rot = RotZ(gamma) * RotY(beta) * RotZ(alpha)
        //   corr = rsh_src(1) * rsh_dst(1);
        //   grad = zeros(1,3);
        //   Hessian =zeros(3,3);
        corr = rshSrc_[0] * rshDst_[0];
        Eigen::MatrixXd corrEig(1, 1);
        corrEig(0, 0) = corr;
        grad = Eigen::Vector3d::Zero();
        Eigen::MatrixXd grad0(1, 1), grad1(1, 1), grad2(1, 1);
        grad0.setZero();
        grad1.setZero();
        grad2.setZero();
        Hessian = Eigen::Matrix3d::Zero();
        //   [M_alpha, dM_alpha, ddM_alpha] = rsh_rot_z(alpha,lmax);
        //   [M_beta, dM_beta, ddM_beta] = rsh_rot_z(beta,lmax);
        //   [M_gamma, dM_gamma, ddM_gamma] = rsh_rot_z(gamma,lmax);
        double alpha = ea(0);
        double beta = ea(1);
        double gamma = ea(2);
        std::vector<Eigen::MatrixXd> M_alpha(lmax_ + 1), dM_alpha(lmax_ + 1), ddM_alpha(lmax_ + 1);
        std::vector<Eigen::MatrixXd> M_beta(lmax_ + 1), dM_beta(lmax_ + 1), ddM_beta(lmax_ + 1);
        std::vector<Eigen::MatrixXd> M_gamma(lmax_ + 1), dM_gamma(lmax_ + 1), ddM_gamma(lmax_ + 1);
        // rsh_rot_z(lmax_, alpha, M_alpha, dM_alpha, ddM_alpha);
        // rsh_rot_z(lmax_,beta, M_beta, dM_beta, ddM_beta);
        // rsh_rot_z(lmax_, gamma, M_gamma, dM_gamma, ddM_gamma);
        //   E = [0 0 1;
        //         1 0 0;
        //         0 1 0];
        //   Me = rsh_rot_matrix(E, lmax);
        Eigen::Matrix3d E;
        E << 0, 0, 1,
            1, 0, 0,
            0, 1, 0;
        std::vector<Eigen::MatrixXd> Me(lmax_ + 1);
        std::vector<Eigen::MatrixXd> Ue(lmax_ + 1), Ve(lmax_ + 1), We(lmax_ + 1); // will be unused in this case, but rsh_rot_matrix "needs" them
        // rsh_rot_matrix(lmax_, E, Me, Ue, Ve, We);
        //   Mf{1} = 1;
        Mf.clear();
        Eigen::MatrixXd oneMat(1, 1);
        oneMat(0, 0) = 1;
        Mf.push_back(oneMat);
        //
        //   for l=1:lmax
        for (int l = 1; l <= lmax_; ++l)
        {
            // ROFL_VAR3("FOR INFINITE LOOP?", l, lmax_);
            // // ROFL_VAR1(l);
            //       rsh_src_l = rsh_src(sh_lm_to_index(l,-l):sh_lm_to_index(l,l));
            //       rsh_dst_l = rsh_dst(sh_lm_to_index(l,-l):sh_lm_to_index(l,l));
            // int start_id = sh_lm_to_index(l, -l) - 1; //-1 due to C++ VS MATLAB indicization starting from 0 VS 1
            // int end_id = sh_lm_to_index(l, l);        //+1 due to a:b in MATLAB including both endpoints
            int start_id, end_id;
            std::vector<double> rsh_src_l(rshSrc_.begin() + start_id, rshSrc_.begin() + end_id);
            std::vector<double> rsh_dst_l(rshDst_.begin() + start_id, rshDst_.begin() + end_id);
            Eigen::Map<Eigen::MatrixXd> rsh_src_l_eig(rsh_src_l.data(), rsh_src_l.size(), 1);
            Eigen::Map<Eigen::MatrixXd> rsh_dst_l_eig(rsh_dst_l.data(), rsh_dst_l.size(), 1);
            //       corr = corr + rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            // // ROFL_VAR6(rsh_src_l_eig.transpose(), rsh_dst_l_eig.transpose(), Me[l].transpose(), M_alpha[l].transpose(), M_beta[l].transpose(), M_gamma[l].transpose());
            // // ROFL_VAR6(rsh_src_l_eig.size(), rsh_dst_l_eig.size(), Me[l].size(), M_alpha[l].size(), M_beta[l].size(), M_gamma[l].size());
            // // ROFL_VAR2(start_id, end_id);
            // // for (int l = 0; l <= lmax_; ++l)
            // // {
            // // ROFL_VAR5(l, Me[l].size(), M_alpha[l].size(), M_beta[l].size(), M_gamma[l].size());
            // // }
            corrEig += rsh_dst_l_eig.transpose() * M_alpha[l] * Me[l].transpose() * M_beta[l] * Me[l] * M_gamma[l] * rsh_src_l_eig;
            //       grad(1) = grad(1) + rsh_dst_l' * dM_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            //       grad(2) = grad(2) + rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * dM_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            //       grad(3) = grad(3) + rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * dM_gamma{l+1} * rsh_src_l;
            grad0 += rsh_dst_l_eig.transpose() * dM_alpha[l] * Me[l].transpose() * M_beta[l] * Me[l] * M_gamma[l] * rsh_src_l_eig;
            grad1 += rsh_dst_l_eig.transpose() * M_alpha[l] * Me[l].transpose() * dM_beta[l] * Me[l] * M_gamma[l] * rsh_src_l_eig;
            grad2 += rsh_dst_l_eig.transpose() * M_alpha[l] * Me[l].transpose() * M_beta[l] * Me[l] * dM_gamma[l] * rsh_src_l_eig;
            //       Hessian(1,1) = Hessian(1,1) +  rsh_dst_l' * ddM_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            //       Hessian(1,2) = Hessian(1,2) +  rsh_dst_l' * dM_alpha{l+1} * ...
            //                Me{l+1}' * dM_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            //       Hessian(1,3) = Hessian(1,3) +  rsh_dst_l' * dM_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * dM_gamma{l+1} * rsh_src_l;
            //       Hessian(2,1) = Hessian(2,1) +  rsh_dst_l' * dM_alpha{l+1} * ...
            //                Me{l+1}' * dM_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            //       Hessian(2,2) = Hessian(2,2) +  rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * ddM_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            //       Hessian(2,3) = Hessian(2,3) +  rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * dM_beta{l+1} * Me{l+1} * dM_gamma{l+1} * rsh_src_l;
            //       Hessian(3,1) = Hessian(3,1) +  rsh_dst_l' * dM_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * dM_gamma{l+1} * rsh_src_l;
            //       Hessian(3,2) = Hessian(3,2) +  rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * dM_beta{l+1} * Me{l+1} * dM_gamma{l+1} * rsh_src_l;
            //       Hessian(3,3) = Hessian(3,3) +  rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * ddM_gamma{l+1} * rsh_src_l;
            //  TODO: implement in C++ after fixing it
            // Mf
            // Mf{l+1} = M_alpha{l+1} * Me{l+1}' * M_beta{l+1} * Me{l+1} * M_gamma{l+1};
            Eigen::MatrixXd MfNew = M_alpha[l] * Me[l].transpose() * M_beta[l] * Me[l] * M_gamma[l];
            Mf.push_back(MfNew);
        }
        corr = corrEig(0, 0);
        grad << grad0, grad1, grad2;
    };

    void SampleSomProblem::rsh_corr_grad(const Eigen::Vector3d &ea,
                                         double &corr,
                                         Eigen::Vector3d &grad,
                                         Eigen::Matrix3d &Hessian) const
    {
        ea.trace();
        // ...
        // function [corr, grad, Hessian, Mf] = rsh_corr_grad(rsh_src,rsh_dst,alpha,beta,gamma,lmax)
        // %
        // % Note: the function assumes that ZYZ Euler angles are used:
        // %  Rot = RotZ(gamma) * RotY(beta) * RotZ(alpha)
        //   corr = rsh_src(1) * rsh_dst(1);
        //   grad = zeros(1,3);
        //   Hessian =zeros(3,3);
        corr = rshSrc_[0] * rshDst_[0];
        Eigen::MatrixXd corrEig(1, 1);
        corrEig(0, 0) = corr;
        grad = Eigen::Vector3d::Zero();
        Eigen::MatrixXd grad0(1, 1), grad1(1, 1), grad2(1, 1);
        grad0.setZero();
        grad1.setZero();
        grad2.setZero();
        Hessian = Eigen::Matrix3d::Zero();
        //   [M_alpha, dM_alpha, ddM_alpha] = rsh_rot_z(alpha,lmax);
        //   [M_beta, dM_beta, ddM_beta] = rsh_rot_z(beta,lmax);
        //   [M_gamma, dM_gamma, ddM_gamma] = rsh_rot_z(gamma,lmax);
        double alpha = ea(0);
        double beta = ea(1);
        double gamma = ea(2);
        std::vector<Eigen::MatrixXd> M_alpha(lmax_ + 1), dM_alpha(lmax_ + 1), ddM_alpha(lmax_ + 1);
        std::vector<Eigen::MatrixXd> M_beta(lmax_ + 1), dM_beta(lmax_ + 1), ddM_beta(lmax_ + 1);
        std::vector<Eigen::MatrixXd> M_gamma(lmax_ + 1), dM_gamma(lmax_ + 1), ddM_gamma(lmax_ + 1);
        // rsh_rot_z(lmax_, alpha, M_alpha, dM_alpha, ddM_alpha);
        // rsh_rot_z(lmax_, beta, M_beta, dM_beta, ddM_beta);
        // rsh_rot_z(lmax_, gamma, M_gamma, dM_gamma, ddM_gamma);
        //   E = [0 0 1;
        //         1 0 0;
        //         0 1 0];
        //   Me = rsh_rot_matrix(E, lmax);
        Eigen::Matrix3d E;
        E << 0, 0, 1,
            1, 0, 0,
            0, 1, 0;
        std::vector<Eigen::MatrixXd> Me(lmax_ + 1);
        std::vector<Eigen::MatrixXd> Ue(lmax_ + 1), Ve(lmax_ + 1), We(lmax_ + 1); // will be unused in this case, but rsh_rot_matrix "needs" them
        // rsh_rot_matrix(lmax_, E, Me, Ue, Ve, We);
        //
        //   for l=1:lmax
        for (int l = 1; l <= lmax_; ++l)
        {
            // ROFL_VAR2("FOR INFINITE LOOP?",l);
            // // ROFL_VAR1(l);
            //       rsh_src_l = rsh_src(sh_lm_to_index(l,-l):sh_lm_to_index(l,l));
            //       rsh_dst_l = rsh_dst(sh_lm_to_index(l,-l):sh_lm_to_index(l,l));
            // int start_id = sh_lm_to_index(l, -l) - 1; //-1 due to C++ VS MATLAB indicization starting from 0 VS 1
            // int end_id = sh_lm_to_index(l, l);        //+1 due to a:b in MATLAB including both endpoints
            int start_id, end_id;
            std::vector<double> rsh_src_l(rshSrc_.begin() + start_id, rshSrc_.begin() + end_id);
            std::vector<double> rsh_dst_l(rshDst_.begin() + start_id, rshDst_.begin() + end_id);
            Eigen::Map<Eigen::MatrixXd> rsh_src_l_eig(rsh_src_l.data(), rsh_src_l.size(), 1);
            Eigen::Map<Eigen::MatrixXd> rsh_dst_l_eig(rsh_dst_l.data(), rsh_dst_l.size(), 1);
            //       corr = corr + rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            // // ROFL_VAR6(rsh_src_l_eig.transpose(), rsh_dst_l_eig.transpose(), Me[l].transpose(), M_alpha[l].transpose(), M_beta[l].transpose(), M_gamma[l].transpose());
            // // ROFL_VAR6(rsh_src_l_eig.size(), rsh_dst_l_eig.size(), Me[l].size(), M_alpha[l].size(), M_beta[l].size(), M_gamma[l].size());
            // // ROFL_VAR2(start_id, end_id);
            // // for (int l = 0; l <= lmax_; ++l)
            // // {
            // // ROFL_VAR5(l, Me[l].size(), M_alpha[l].size(), M_beta[l].size(), M_gamma[l].size());
            // // }
            corrEig += rsh_dst_l_eig.transpose() * M_alpha[l] * Me[l].transpose() * M_beta[l] * Me[l] * M_gamma[l] * rsh_src_l_eig;
            //       grad(1) = grad(1) + rsh_dst_l' * dM_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            //       grad(2) = grad(2) + rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * dM_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            //       grad(3) = grad(3) + rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * dM_gamma{l+1} * rsh_src_l;
            grad0 += rsh_dst_l_eig.transpose() * dM_alpha[l] * Me[l].transpose() * M_beta[l] * Me[l] * M_gamma[l] * rsh_src_l_eig;
            grad1 += rsh_dst_l_eig.transpose() * M_alpha[l] * Me[l].transpose() * dM_beta[l] * Me[l] * M_gamma[l] * rsh_src_l_eig;
            grad2 += rsh_dst_l_eig.transpose() * M_alpha[l] * Me[l].transpose() * M_beta[l] * Me[l] * dM_gamma[l] * rsh_src_l_eig;
            //       Hessian(1,1) = Hessian(1,1) +  rsh_dst_l' * ddM_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            //       Hessian(1,2) = Hessian(1,2) +  rsh_dst_l' * dM_alpha{l+1} * ...
            //                Me{l+1}' * dM_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            //       Hessian(1,3) = Hessian(1,3) +  rsh_dst_l' * dM_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * dM_gamma{l+1} * rsh_src_l;
            //       Hessian(2,1) = Hessian(2,1) +  rsh_dst_l' * dM_alpha{l+1} * ...
            //                Me{l+1}' * dM_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            //       Hessian(2,2) = Hessian(2,2) +  rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * ddM_beta{l+1} * Me{l+1} * M_gamma{l+1} * rsh_src_l;
            //       Hessian(2,3) = Hessian(2,3) +  rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * dM_beta{l+1} * Me{l+1} * dM_gamma{l+1} * rsh_src_l;
            //       Hessian(3,1) = Hessian(3,1) +  rsh_dst_l' * dM_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * dM_gamma{l+1} * rsh_src_l;
            //       Hessian(3,2) = Hessian(3,2) +  rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * dM_beta{l+1} * Me{l+1} * dM_gamma{l+1} * rsh_src_l;
            //       Hessian(3,3) = Hessian(3,3) +  rsh_dst_l' * M_alpha{l+1} * ...
            //                Me{l+1}' * M_beta{l+1} * Me{l+1} * ddM_gamma{l+1} * rsh_src_l;
            //  TODO: implement in C++ after fixing it
        }
        corr = corrEig(0, 0);
        grad << grad0, grad1, grad2;
        // ROFL_VAR2(corr, grad.transpose());
    };

}
