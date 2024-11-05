#include "problems_SampleSom.h"

namespace ROPTLIB
{
    SampleSomProblem::SampleSomProblem() {}

    SampleSomProblem::SampleSomProblem(SomSize somSz, Eigen::MatrixXf &Tijs, Eigen::MatrixXi &edges)
    {
        sz_ = somSz;
        Tijs_ = Tijs;
        edges_ = edges;
        e_ = Tijs_.cols();
    }

    SampleSomProblem::~SampleSomProblem() {};

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
        // arsRshRot.create(e_, xEigen);
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
        // std::vector<Eigen::MatrixXd> Mf(e_+1); // unused
        // Mf.reserve(e_+1);
        // rsh_corr_grad(xEulAngles, corr, grad, Hessian, Mf);
        rsh_corr_grad(xEulAngles, corr, grad, Hessian);

        // J = euler_angles_jacobian_zyz(alpha,beta,gamma);
        // pinvJ = pinv(J);
        Eigen::MatrixXd J(9, 3);
        J.setZero();
        // eulerAnglesJacobianZYZ(e_, xEulAngles, J);
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
        
        // corr = rshSrc_[0] * rshDst_[0];
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
        std::vector<Eigen::MatrixXd> M_alpha(e_ + 1), dM_alpha(e_ + 1), ddM_alpha(e_ + 1);
        std::vector<Eigen::MatrixXd> M_beta(e_ + 1), dM_beta(e_ + 1), ddM_beta(e_ + 1);
        std::vector<Eigen::MatrixXd> M_gamma(e_ + 1), dM_gamma(e_ + 1), ddM_gamma(e_ + 1);
        
        Eigen::Matrix3d E;
        E << 0, 0, 1,
            1, 0, 0,
            0, 1, 0;
        std::vector<Eigen::MatrixXd> Me(e_ + 1);
        std::vector<Eigen::MatrixXd> Ue(e_ + 1), Ve(e_ + 1), We(e_ + 1); // will be unused in this case, but rsh_rot_matrix "needs" them
        // rsh_rot_matrix(e_, E, Me, Ue, Ve, We);
        //   Mf{1} = 1;
        Mf.clear();
        Eigen::MatrixXd oneMat(1, 1);
        oneMat(0, 0) = 1;
        Mf.push_back(oneMat);
        //
        //   for l=1:lmax
        for (int l = 1; l <= e_; ++l)
        {
            //
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
        corr = edges_(0,0) * edges_(0,0);
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
        std::vector<Eigen::MatrixXd> M_alpha(e_ + 1), dM_alpha(e_ + 1), ddM_alpha(e_ + 1);
        std::vector<Eigen::MatrixXd> M_beta(e_ + 1), dM_beta(e_ + 1), ddM_beta(e_ + 1);
        std::vector<Eigen::MatrixXd> M_gamma(e_ + 1), dM_gamma(e_ + 1), ddM_gamma(e_ + 1);
        // rsh_rot_z(e_, alpha, M_alpha, dM_alpha, ddM_alpha);
        // rsh_rot_z(e_, beta, M_beta, dM_beta, ddM_beta);
        // rsh_rot_z(e_, gamma, M_gamma, dM_gamma, ddM_gamma);
        //   E = [0 0 1;
        //         1 0 0;
        //         0 1 0];
        //   Me = rsh_rot_matrix(E, lmax);
        Eigen::Matrix3d E;
        E << 0, 0, 1,
            1, 0, 0,
            0, 1, 0;
        std::vector<Eigen::MatrixXd> Me(e_ + 1);
        std::vector<Eigen::MatrixXd> Ue(e_ + 1), Ve(e_ + 1), We(e_ + 1); // will be unused in this case, but rsh_rot_matrix "needs" them
        // rsh_rot_matrix(e_, E, Me, Ue, Ve, We);
        //
        //   for l=1:lmax
        for (int l = 1; l <= e_; ++l)
        {
            //
        }
        corr = corrEig(0, 0);
        grad << grad0, grad1, grad2;
        // ROFL_VAR2(corr, grad.transpose());
    };

}
