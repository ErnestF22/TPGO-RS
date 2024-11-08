#include "problems_SampleSom.h"

namespace ROPTLIB
{
    SampleSomProblem::SampleSomProblem() {}

    SampleSomProblem::SampleSomProblem(SomSize somSz, MatD &Tijs, Eigen::MatrixXi &edges)
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
        MatD xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);
        // ROFL_VAR1(x);

        realdp cost = 0.0;

        for (int e = 0; e < numEdges_; ++e)
        {
            MatD Ri(sz_.p_, sz_.d_);
            MatD Ti(sz_.p_, 1);
            MatD Tj(sz_.p_, 1);

            Eigen::VectorXd tij(sz_.d_);
            tij = Tijs_.col(e);

            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1
            getRi(xEigen, Ri, i);
            getTi(xEigen, Ti, i);
            getTi(xEigen, Tj, j);

            // ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

            double cost_e = (Ri * tij - Tj + Ti).squaredNorm();

            cost += cost_e;
        }

        ROFL_VAR1(cost);
        // ROFL_ASSERT(!std::isnan(corr));

        // Vector *resultEgrad;
        // *resultEgrad = Domain->RandominManifold();
        // EucGrad(x, resultEgrad);
        // x.AddToFields("EGrad", *resultEgrad); //x should have been const?? maybe only its reference

        return cost; // checked -> the - here should be OK
    };

    // Vector &SampleSomProblem::EucGrad(const Variable &x, Vector *result) const
    // {
    //     // result->NewMemoryOnWrite();
    //     // result->GetElement(0) = x.Field("B1x1D1");
    //     // result->GetElement(1) = x.Field("B2x2D2");
    //     // result->GetElement(2) = x.Field("B3x3D3");
    //     // Domain->ScalarTimesVector(x, 2, *result, result);

    //     // result->Print("Euc Grad result: printing it before NewMemoryOnWrite()");
    //     *result = Domain->RandominManifold();
    //     result->NewMemoryOnWrite();
    //     result->Print("Euc Grad result: printing it after NewMemoryOnWrite()");

    //     MatD xEig(fullSz_, 1);
    //     RoptToEig(x, xEig);

    //     VecMatD R(sz_.n_, MatD::Zero(sz_.p_, sz_.d_));
    //     getRotations(xEig, R);

    //     MatD T(MatD::Zero(sz_.p_, sz_.n_));
    //     getTranslations(xEig, T);

    //     // P = zeros(nrs, d*N);
    //     MatD P(MatD::Zero(sz_.p_, sz_.d_ * sz_.n_));
    //     double frct = 0.0;

    //     // LR = zeros(N,N);
    //     // PR = zeros(N,nrs);
    //     // BR_const = zeros(d,d);
    //     MatD Lr(MatD::Zero(sz_.n_, sz_.n_));
    //     MatD Pr(MatD::Zero(sz_.n_, sz_.p_));
    //     MatD Br(MatD::Zero(sz_.d_, sz_.d_));

    //     makePfrct(T, P, frct);
    //     makeLrPrBr(R, Lr, Pr, Br);

    //     VecMatD egR(sz_.n_, MatD::Zero(sz_.p_, sz_.d_));
    //     egradR(P, egR);

    //     MatD egT(MatD::Zero(sz_.p_, sz_.n_));
    //     egradT(T, Lr, Pr, egT);

    //     // result->NewMemoryOnWrite();
    //     // result = Domain->RandomInManifold();

    //     int rotSz = getRotSz();
    //     int translSz = getTranslSz();
    //     int gElemIdx = 0;

    //     // fill result with computed gradient values : R
    //     for (int i = 0; i < sz_.n_; ++i)
    //     {
    //         // ROFL_VAR1(gElemIdx);
    //         // ROFL_VAR2("\n", egR[gElemIdx]);
    //         result->GetElement(gElemIdx).SetToZeros(); // Ri
    //         // result->GetElement(gElemIdx).Print("Ri before assignment");

    //         Vector egRiVec(sz_.p_, sz_.d_);
    //         // egRiVec.Initialize();
    //         realdp *GroptlibWriteArray = egRiVec.ObtainWriteEntireData();
    //         for (int j = 0; j < rotSz; ++j)
    //         {
    //             // ROFL_VAR2(i, j);
    //             // egRiVec.Print("egRiVec before assignment");

    //             // ROFL_VAR1(egRiVec.GetElement(j, 0));

    //             GroptlibWriteArray[j] = egR[i].reshaped(sz_.d_ * sz_.p_, 1)(j);

    //             // ROFL_VAR1("");
    //             // egRiVec.Print("egRiVec after assignment");
    //         }
    //         egRiVec.CopyTo(result->GetElement(gElemIdx));
    //         // result->GetElement(gElemIdx).Print("grad Ri after assignment");
    //         gElemIdx++;
    //     }

    //     // fill result with computed gradient values : T

    //     Vector egTiVec(sz_.p_, sz_.n_);
    //     realdp *GroptlibWriteArray = egTiVec.ObtainWriteEntireData();
    //     for (int j = 0; j < sz_.p_ * sz_.n_; ++j)
    //     {
    //         // egTiVec.Print("egTiVec before assignment");

    //         // ROFL_VAR1(egRiVec.GetElement(j, 0));

    //         GroptlibWriteArray[j] = egT.reshaped(sz_.n_ * sz_.p_, 1)(j);

    //         // ROFL_VAR1("");
    //         // egTiVec.Print("egTiVec after assignment");
    //     }
    //     egTiVec.CopyTo(result->GetElement(gElemIdx));
    //     // result->GetElement(gElemIdx).Print("Eucl. grad Ti after assignment");

    //     // ROFL_VAR2("\n", egT);
    //     // ROFL_VAR1(gElemIdx);
    //     // result->GetElement(gElemIdx).Print();

    //     // result->NewMemoryOnWrite();
    //     // result->SetToZeros();
    //     // *result = Groptlib;

    //     result->Print("printing final result");

    //     // ROFL_ASSERT(0);

    //     x.AddToFields("EGrad", *result); // x should have been const?? maybe only its reference

    //     return *result;
    // };

    Vector &SampleSomProblem::RieGrad(const Variable &x, Vector *result) const
    {
        // result->NewMemoryOnWrite();
        // result->GetElement(0) = x.Field("B1x1D1");
        // result->GetElement(1) = x.Field("B2x2D2");
        // result->GetElement(2) = x.Field("B3x3D3");
        // Domain->ScalarTimesVector(x, 2, *result, result);

        result->Print("RieGrad: printing result at start of function (should be empty)");

        MatD xEig(fullSz_, 1);
        RoptToEig(x, xEig);

        VecMatD R(sz_.n_, MatD::Zero(sz_.p_, sz_.d_));
        getRotations(xEig, R);

        MatD T(MatD::Zero(sz_.p_, sz_.n_));
        getTranslations(xEig, T);

        // P = zeros(nrs, d*N);
        MatD P(MatD::Zero(sz_.p_, sz_.d_ * sz_.n_));
        double frct = 0.0;

        // LR = zeros(N,N);
        // PR = zeros(N,nrs);
        // BR_const = zeros(d,d);
        MatD Lr(MatD::Zero(sz_.n_, sz_.n_));
        MatD Pr(MatD::Zero(sz_.n_, sz_.p_));
        MatD Br(MatD::Zero(sz_.d_, sz_.d_));

        makePfrct(T, P, frct);
        makeLrPrBr(R, Lr, Pr, Br);

        VecMatD rgR(sz_.n_, MatD::Zero(sz_.p_, sz_.d_));
        rgradR(R, P, rgR);

        MatD rgT(MatD::Zero(sz_.p_, sz_.n_));
        rgradT(T, Lr, Pr, rgT);

        // result->NewMemoryOnWrite();
        // result = Domain->RandomInManifold();

        int rotSz = getRotSz();
        int translSz = getTranslSz();
        int gElemIdx = 0;

        // fill result with computed gradient values : R
        for (int i = 0; i < sz_.n_; ++i)
        {
            // ROFL_VAR1(gElemIdx);
            // ROFL_VAR2("\n", rgR[gElemIdx]);
            // result->GetElement(gElemIdx).SetToIdentity(); // Ri
            // result->GetElement(gElemIdx).Print("Ri before assignment");

            Vector rgRiVec(sz_.p_, sz_.d_);
            // rgRiVec.Initialize();
            realdp *GroptlibWriteArray = rgRiVec.ObtainWriteEntireData();
            for (int j = 0; j < rotSz; ++j)
            {
                // ROFL_VAR2(i, j);
                // rgRiVec.Print("rgRiVec before assignment");

                // ROFL_VAR1(rgRiVec.GetElement(j, 0));

                GroptlibWriteArray[j] = rgR[i].reshaped(sz_.d_ * sz_.p_, 1)(j);

                // ROFL_VAR1("");
                // rgRiVec.Print("rgRiVec after assignment");
            }
            rgRiVec.CopyTo(result->GetElement(gElemIdx));
            // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
            gElemIdx++;
        }

        // fill result with computed gradient values : T

        Vector rgTiVec(sz_.p_, sz_.n_);
        realdp *GroptlibWriteArray = rgTiVec.ObtainWriteEntireData();
        for (int j = 0; j < sz_.p_ * sz_.n_; ++j)
        {
            // rgTiVec.Print("rgTiVec before assignment");

            // ROFL_VAR1(rgRiVec.GetElement(j, 0));

            GroptlibWriteArray[j] = rgT.reshaped(sz_.n_ * sz_.p_, 1)(j);

            // ROFL_VAR1("");
            // rgTiVec.Print("rgTiVec after assignment");
        }
        rgTiVec.CopyTo(result->GetElement(gElemIdx));
        result->GetElement(gElemIdx).Print("grad Ti after assignment");

        // ROFL_VAR2("\n", rgT);
        // ROFL_VAR1(gElemIdx);
        // result->GetElement(gElemIdx).Print();

        // result->NewMemoryOnWrite();
        // result->SetToZeros();
        // *result = Groptlib;

        result->Print("printing final result");

        // ROFL_ASSERT(0);

        return *result;
    };

    // Vector &SampleSomProblem::Grad(const Variable &x, Vector *result) const
    // {
    //     // result->NewMemoryOnWrite();
    //     // result->GetElement(0) = x.Field("B1x1D1");
    //     // result->GetElement(1) = x.Field("B2x2D2");
    //     // result->GetElement(2) = x.Field("B3x3D3");
    //     // Domain->ScalarTimesVector(x, 2, *result, result);

    //     result->Print("Grad: printing it before anything happens result");

    //     MatD xEig(fullSz_, 1);
    //     RoptToEig(x, xEig);

    //     VecMatD R(sz_.n_, MatD::Zero(sz_.p_, sz_.d_));
    //     getRotations(xEig, R);

    //     MatD T(MatD::Zero(sz_.p_, sz_.n_));
    //     getTranslations(xEig, T);

    //     // P = zeros(nrs, d*N);
    //     MatD P(MatD::Zero(sz_.p_, sz_.d_ * sz_.n_));
    //     double frct = 0.0;

    //     // LR = zeros(N,N);
    //     // PR = zeros(N,nrs);
    //     // BR_const = zeros(d,d);
    //     MatD Lr(MatD::Zero(sz_.n_, sz_.n_));
    //     MatD Pr(MatD::Zero(sz_.n_, sz_.p_));
    //     MatD Br(MatD::Zero(sz_.d_, sz_.d_));

    //     makePfrct(T, P, frct);
    //     makeLrPrBr(R, Lr, Pr, Br);

    //     VecMatD egR(sz_.n_, MatD::Zero(sz_.p_, sz_.d_));
    //     egradR(P, egR);

    //     MatD rgT(MatD::Zero(sz_.p_, sz_.n_));
    //     rgradT(T, Lr, Pr, rgT);

    //     // result->NewMemoryOnWrite();
    //     // result = Domain->RandomInManifold();

    //     int rotSz = getRotSz();
    //     int translSz = getTranslSz();
    //     int gElemIdx = 0;

    //     // fill result with computed gradient values : R
    //     for (int i = 0; i < sz_.n_; ++i)
    //     {
    //         // ROFL_VAR1(gElemIdx);
    //         // ROFL_VAR2("\n", egR[gElemIdx]);
    //         result->GetElement(gElemIdx).SetToZeros(); // Ri
    //         // result->GetElement(gElemIdx).Print("Ri before assignment");

    //         Vector egRiVec(sz_.p_, sz_.d_);
    //         // egRiVec.Initialize();
    //         realdp *GroptlibWriteArray = egRiVec.ObtainWriteEntireData();
    //         for (int j = 0; j < rotSz; ++j)
    //         {
    //             // ROFL_VAR2(i, j);
    //             // egRiVec.Print("egRiVec before assignment");

    //             // ROFL_VAR1(egRiVec.GetElement(j, 0));

    //             GroptlibWriteArray[j] = egR[i].reshaped(sz_.d_ * sz_.p_, 1)(j);

    //             // ROFL_VAR1("");
    //             // egRiVec.Print("egRiVec after assignment");
    //         }

    //         { // eucgrad2rgrad scope!!

    //             class TmpProblem : public Problem
    //             {
    //             public:
    //                 TmpProblem(int p, int d) : p_(p), d_(d) {};

    //                 virtual ~TmpProblem() {};

    //                 virtual realdp f(const Variable &x) const { return 0; };

    //                 virtual Vector &Grad(const Variable &x, Vector *result) const { return *result; };

    //             private:
    //                 int d_;
    //                 int p_;
    //             };

    //             Vector rgRiVec(6, 1, 1);

    //             Stiefel maniSt(sz_.p_, sz_.d_);
    //             maniSt.ChooseParamsSet1();

    //             // Set the domain of the problem to be the product of Stiefel manifolds
    //             TmpProblem Prob(sz_.p_, sz_.d_);
    //             Prob.SetDomain(&maniSt);
    //             Domain->EucGradToGrad(x.GetElement(gElemIdx), egRiVec, &Prob, &rgRiVec);
    //             rgRiVec.CopyTo(result->GetElement(gElemIdx));
    //             rgRiVec.Print("scope rgRiVec");
    //             result->GetElement(gElemIdx).Print("scope result");
    //         }

    //         // result->GetElement(gElemIdx).Print("grad Ri after assignment");
    //         gElemIdx++;
    //     }

    //     // fill result with computed gradient values : T

    //     Vector rgTiVec(sz_.p_, sz_.n_);
    //     realdp *GroptlibWriteArray = rgTiVec.ObtainWriteEntireData();
    //     for (int j = 0; j < sz_.p_ * sz_.n_; ++j)
    //     {
    //         // rgTiVec.Print("rgTiVec before assignment");

    //         // ROFL_VAR1(rgRiVec.GetElement(j, 0));

    //         GroptlibWriteArray[j] = rgT.reshaped(sz_.n_ * sz_.p_, 1)(j);

    //         // ROFL_VAR1("");
    //         // rgTiVec.Print("rgTiVec after assignment");
    //     }
    //     rgTiVec.CopyTo(result->GetElement(gElemIdx));
    //     result->GetElement(gElemIdx).Print("grad Ti after assignment");

    //     // ROFL_VAR2("\n", rgT);
    //     // ROFL_VAR1(gElemIdx);
    //     // result->GetElement(gElemIdx).Print();

    //     // result->NewMemoryOnWrite();
    //     // result->SetToZeros();
    //     // *result = Groptlib;

    //     result->Print("printing final result");

    //     // ROFL_ASSERT(0);

    //     return *result;
    // };

    void SampleSomProblem::egradR(const MatD &P, VecMatD &egR) const
    {
        unStackH(P, egR, sz_.d_);
    }

    void SampleSomProblem::rgradR(const VecMatD &R, const MatD &P, VecMatD &rgR) const
    {
        VecMatD egR;
        egradR(P, egR);

        stiefelTangentProj(R, egR, rgR);
    }

    void SampleSomProblem::egradT(const MatD &T, const MatD &Lr, const MatD &Pr, MatD &egT) const
    {
        // ROFL_VAR7(egT, T.rows(), T.cols(), Lr.rows(), Lr.cols(), Pr.rows(), Pr.cols());
        egT = T * (Lr + Lr.transpose()) + Pr.transpose();
    }

    void SampleSomProblem::rgradT(const MatD &T, const MatD &Lr, const MatD &Pr, MatD &egT) const
    {
        egradT(T, Lr, Pr, egT);
    }

    // Vector &SampleSomProblem::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
    // {
    //     // TODO: implement

    //     result->NewMemoryOnWrite();

    //     // h.R = rsom_rhess_rot_stiefel(R, Rdot, problem_structs) + hrt;
    //     // h.T = rsom_rhess_transl_stiefel(T, Tdot, problem_structs) + htr;

    //     result->Print("EucHessianEta: printing it just after NewMemoryOnWrite()");

    //     return *result;
    // };

    Vector &SampleSomProblem::HessianEta(const Variable &x, const Vector &etax, Vector *result) const
    {
        // TODO: implement

        result->NewMemoryOnWrite();

        result->Print("HessianEta: printing it just after NewMemoryOnWrite()");

        return *result;
    };

    void SampleSomProblem::RoptToEig(Vector x, MatD &xEigen) const
    {
        Vector xT = x.GetTranspose(); // Eigen ADV init is row-major!!

        int totSz = xEigen.rows();

        const realdp *xArr = xT.ObtainWriteEntireData();
        for (int i = 0; i < totSz; ++i)
            xEigen(i) = xArr[i];
    }

    void SampleSomProblem::vstack(const VecMatD &in, MatD &out) const
    {
        // out has same number of columns, whereas number of rows is the product of the size of the other 2 dimensions

        int rowJump = sz_.p_; // TODO: generalize later
        for (int i = 0; i < in.size(); ++i)
        {
            out.block(rowJump * i, 0, rowJump, sz_.d_) = in[i];
        }
    }

    void SampleSomProblem::hstack(const VecMatD &in, MatD &out) const
    {
        // out has same number of columns, whereas number of rows is the product of the size of the other 2 dimensions

        int colJump = sz_.d_; // TODO: generalize later
        for (int i = 0; i < in.size(); ++i)
        {
            out.block(0, colJump * i, sz_.p_, colJump) = in[i];
        }
    }

    void SampleSomProblem::unStackV(const MatD &in, VecMatD &out, int rowsOut) const
    {
        int n = (int)in.rows() / rowsOut;

        int fixedSz = sz_.d_; // size that does not change in the 3D->2D transition (here, number of columns)

        // ROFL_VAR3(n, rowsOut, in.rows());
        ROFL_ASSERT(n * rowsOut == in.rows());

        out.clear();
        out.resize(n, MatD::Zero(rowsOut, in.cols()));

        for (int i = 0; i < n; ++i)
        {
            out[i] = in.block(rowsOut * i, 0, rowsOut, fixedSz);
        }
    }

    void SampleSomProblem::unStackH(const MatD &in, VecMatD &out, int colsOut) const
    {
        int n = (int)in.cols() / colsOut;

        int fixedSz = sz_.p_; // size that does not change in the 3D->2D transition (here, number of rows)

        ROFL_ASSERT(n * colsOut == in.cols());

        out.clear();
        out.resize(n, MatD::Zero(in.rows(), colsOut));

        for (int i = 0; i < n; ++i)
        {
            out[i] = in.block(0, colsOut * i, fixedSz, colsOut);
        }
    }

    void SampleSomProblem::getRi(const Variable &x, MatD &rOut, int i) const
    {
        MatD xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);

        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        int startId = i * rotSz;
        // int endId = (i+1) * rotSz;

        MatD rOutVec(rotSz, 1);
        rOutVec = xEigen.block(startId, 0, rotSz, 1);

        rOut = rOutVec.reshaped(sz_.p_, sz_.d_);
    }

    void SampleSomProblem::getRi(const MatD &xEig, MatD &rOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        int startId = i * rotSz;
        // int endId = (i+1) * rotSz;

        MatD rOutVec(rotSz, 1);
        rOutVec = xEig.block(startId, 0, rotSz, 1);

        rOut = rOutVec.reshaped(sz_.p_, sz_.d_);
    }

    void SampleSomProblem::getRotations(const MatD &xEig, VecMatD &rOut) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        // int endId = (i+1) * rotSz;

        rOut.clear();

        for (int i = 0; i < sz_.n_; ++i)
        {
            int startId = i * rotSz;
            MatD Ri(sz_.p_, sz_.d_);
            getRi(xEig, Ri, i); // TODO: this can probably be optimized better
            rOut.push_back(Ri);
        }
    }

    void SampleSomProblem::getTi(const Variable &x, MatD &tOut, int i) const
    {
        MatD xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);

        // rOut already needs to have fixed size by here
        int rotSz = getRotSz();
        int translSz = getTranslSz();

        int startId = sz_.n_ * rotSz + i * translSz;
        // int endId = (i+1) * rotSz;

        // ROFL_VAR1(startId);

        tOut = xEigen.block(startId, 0, translSz, 1);
    }

    void SampleSomProblem::getTi(const MatD &xEig, MatD &tOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz();
        int translSz = getTranslSz();

        int startId = sz_.n_ * rotSz + i * translSz;
        // int endId = (i+1) * rotSz;

        // ROFL_VAR1(startId);

        tOut = xEig.block(startId, 0, translSz, 1);
    }

    void SampleSomProblem::getTranslations(const MatD &xEig, MatD &tOut) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        // int endId = (i+1) * rotSz;

        for (int i = 0; i < sz_.n_; ++i)
        {
            int startId = i * rotSz;
            MatD Ti(sz_.p_, 1);
            getTi(xEig, Ti, i); // TODO: this can probably be optimized better
            tOut.col(i) = Ti;
        }
    }

    void SampleSomProblem::makePfrct(const MatD &T, MatD &P, double &frct) const
    {
        // P = zeros(nrs, d*N);

        // TODO: add some size checks at least for P

        P.setZero();
        frct = 0.0;

        for (int e = 0; e < numEdges_; ++e)
        {
            int ii = edges_(e, 0) - 1;
            int jj = edges_(e, 1) - 1;

            MatD T_j = T.col(jj);
            MatD T_i = T.col(ii);
            auto tij = Tijs_.col(e);
            auto Pe = 2 * (T_i * tij.transpose() - T_j * tij.transpose());
            P.block(0, ii * sz_.d_, sz_.p_, sz_.d_) += Pe;

            auto c = T_i * T_i.transpose() + T_j * T_j.transpose() - T_i * T_j.transpose() - T_j * T_i.transpose();
            auto d = tij * tij.transpose();
            frct += c.trace() + d.trace();
        }
    }

    void SampleSomProblem::makeLrPrBr(const VecMatD &R, MatD &Lr, MatD &Pr, MatD &Br) const
    {
        // LR = zeros(N,N);
        // PR = zeros(N,nrs);
        // BR_const = zeros(d,d);

        // TODO: add some size checks at least for Lr, Pr, Br

        Lr.setZero();
        Pr.setZero();
        Br.setZero();

        for (int e = 0; e < numEdges_; ++e)
        {
            int ii = edges_(e, 0) - 1;
            int jj = edges_(e, 1) - 1;
            auto Ri = R[ii];

            MatD bij(MatD::Zero(sz_.n_, 1));

            bij(ii) = 1;
            bij(jj) = -1;

            auto tij = Tijs_.col(e);

            Lr += bij * bij.transpose();

            Pr += 2 * bij * tij.transpose() * Ri.transpose();

            Br += tij * tij.transpose();
        }
    }

    int SampleSomProblem::getRotSz() const
    {
        return sz_.d_ * sz_.p_;
    }

    int SampleSomProblem::getTranslSz() const
    {
        return sz_.p_;
    }

    void SampleSomProblem::stiefelTangentProj(const VecMatD &Y, const VecMatD &Hin, VecMatD &Hout) const
    {
        ROFL_ASSERT(Y.size() == sz_.n_);

        MatD tmp = Y[0].transpose() * Hin[0];
        Hout.clear();
        Hout.resize(sz_.n_, MatD::Zero(tmp.rows(), tmp.cols()));

        for (int i = 0; i < sz_.n_; ++i)
        {
            stiefelTangentProj(Y[i], Hin[i], Hout[i]);
        }
    }

    void SampleSomProblem::stiefelTangentProj(const MatD &Y, const MatD &Hin, MatD &Hout) const
    {
        MatD tmp = Y.transpose() * Hin;

        ROFL_ASSERT(tmp.rows() == tmp.cols()); // kind of useless as the check is performed also in extractSymmetricPart()

        MatD sympart(MatD::Zero(tmp.rows(), tmp.cols()));
        extractSymmetricPart(Y.transpose() * Hin, sympart);
        Hout = Hin - Y * sympart;
    }

    void SampleSomProblem::extractSymmetricPart(const MatD &in, MatD &out) const
    {
        ROFL_ASSERT(in.rows() == in.cols());
        out = 0.5 * (in + in.transpose());
    }

} // end of namespace ROPTLIB
