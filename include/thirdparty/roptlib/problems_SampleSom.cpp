#include "problems_SampleSom.h"

namespace ROPTLIB
{
    SampleSomProblem::SampleSomProblem() {}

    SampleSomProblem::SampleSomProblem(SomSize somSz, Eigen::MatrixXd &Tijs, Eigen::MatrixXi &edges)
    {
        sz_ = somSz;
        Tijs_ = Tijs;
        edges_ = edges;
        numEdges_ = Tijs_.cols();
        fullSz_ = sz_.d_ * sz_.p_ * sz_.n_ + sz_.p_ * sz_.n_;
    }

    SampleSomProblem::~SampleSomProblem() {};

    realdp SampleSomProblem::f(const Variable &x) const
    {
        Eigen::MatrixXd xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);
        ROFL_VAR1(x);

        realdp cost = 0.0;

        for (int e = 0; e < numEdges_; ++e)
        {
            Eigen::MatrixXd Ri(sz_.p_, sz_.d_);
            Eigen::MatrixXd Ti(sz_.p_, 1);
            Eigen::MatrixXd Tj(sz_.p_, 1);

            Eigen::VectorXd tij(sz_.d_);
            tij = Tijs_.col(e);

            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1
            getRi(xEigen, Ri, i);
            getTi(xEigen, Ti, i);
            getTi(xEigen, Tj, j);

            ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

            double cost_e = (Ri * tij - Tj + Ti).squaredNorm();

            cost += cost_e;
        }

        ROFL_VAR1(cost);
        // ROFL_ASSERT(!std::isnan(corr));

        return cost; // checked -> the - here should be OK
    };

    // Vector &SampleSomProblem::EucGrad(const Variable &x, Vector *result) const
    // {
    //     // eul_angles = rotm2eul(x, 'ZYZ');
    //     // alpha = eul_angles(1);
    //     // beta = eul_angles(2);
    //     // gamma = eul_angles(3);
    //     Eigen::MatrixXd xEigen(fullSz_, 1);
    //     RoptToEig(x, xEigen);

    //     Eigen::RowVectorXd Gvec(9);
    //     // %     G = zeros(3,3);
    //     Eigen::Matrix3d G(Gvec.data());
    //     G.resize(3, 3);
    //     G.transposeInPlace();
        

    //     // back to Roptlib
    //     Vector Groptlib = Rotations(3).RandominManifold();
    //     Groptlib.ObtainWriteEntireData();
    //     realdp *GroptlibWriteArray = Groptlib.ObtainWriteEntireData(); //!! assignment in col-major order
    //     // Groptlib
    //     GroptlibWriteArray[0] = G(0, 0);
    //     GroptlibWriteArray[1] = G(0, 1);
    //     GroptlibWriteArray[2] = G(0, 2);
    //     GroptlibWriteArray[3] = G(1, 0);
    //     GroptlibWriteArray[4] = G(1, 1);
    //     GroptlibWriteArray[5] = G(1, 2);
    //     GroptlibWriteArray[6] = G(2, 0);
    //     GroptlibWriteArray[7] = G(2, 1);
    //     GroptlibWriteArray[8] = G(2, 2);
    //     // Groptlib.Print("Groptlib"); //This was causing a loop
    //     ROFL_VAR1(xEigen);
    //     ROFL_VAR1(G);

    //     result->NewMemoryOnWrite();
    //     result->SetToZeros();
    //     *result = Groptlib;
    //     return *result;
    // };

    // void SampleSomProblem::RoptToEig(Vector x, Eigen::MatrixXf &xEigen) const
    // {
    //     Vector xT = x.GetTranspose(); // Eigen ADV init is row-major!!

    //     int totSz = xEigen.rows();

    //     const realdp *xArr = xT.ObtainWriteEntireData();
    //     for (int i = 0; i < totSz; ++i)
    //         xEigen(i) = xArr[i];
    // }

    void SampleSomProblem::RoptToEig(Vector x, Eigen::MatrixXd &xEigen) const
    {
        Vector xT = x.GetTranspose(); // Eigen ADV init is row-major!!

        int totSz = xEigen.rows();

        const realdp *xArr = xT.ObtainWriteEntireData();
        for (int i = 0; i < totSz; ++i)
            xEigen(i) = xArr[i];
    }

    void SampleSomProblem::getRi(const Variable &x, Eigen::MatrixXd &rOut, int i) const
    {
        Eigen::MatrixXd xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);

        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        int startId = i * rotSz;
        // int endId = (i+1) * rotSz;

        Eigen::MatrixXd rOutVec(rotSz, 1);
        rOutVec = xEigen.block(startId, 0, rotSz, 1);

        rOut = rOutVec.reshaped(sz_.p_, sz_.d_);
    }

    void SampleSomProblem::getRi(const Eigen::MatrixXd &xEig, Eigen::MatrixXd &rOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        int startId = i * rotSz;
        // int endId = (i+1) * rotSz;

        Eigen::MatrixXd rOutVec(rotSz, 1);
        rOutVec = xEig.block(startId, 0, rotSz, 1);

        rOut = rOutVec.reshaped(sz_.p_, sz_.d_);
    }

    void SampleSomProblem::getTi(const Variable &x, Eigen::MatrixXd &tOut, int i) const
    {
        Eigen::MatrixXd xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);

        // rOut already needs to have fixed size by here
        int rotSz = getRotSz();
        int translSz = getTranslSz();

        int startId = sz_.n_ * rotSz + i * translSz;
        // int endId = (i+1) * rotSz;

        ROFL_VAR1(startId);

        tOut = xEigen.block(startId, 0, translSz, 1);
    }


    void SampleSomProblem::getTi(const Eigen::MatrixXd &xEig, Eigen::MatrixXd &tOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz();
        int translSz = getTranslSz();

        int startId = sz_.n_ * rotSz + i * translSz;
        // int endId = (i+1) * rotSz;

        ROFL_VAR1(startId);

        tOut = xEig.block(startId, 0, translSz, 1);
    }

    int SampleSomProblem::getRotSz() const
    {
        return sz_.d_ * sz_.p_;
    }

    int SampleSomProblem::getTranslSz() const
    {
        return sz_.p_;
    }

    // Vector &SampleSomProblem::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
    // {
    //     //TODO: implement
    //     return *result;
    // };

}
