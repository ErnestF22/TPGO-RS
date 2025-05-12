#include "problems_Ssom.h"

namespace ROPTLIB
{
    SsomProblem::SsomProblem() {}

    SsomProblem::SsomProblem(SomUtils::SomSize somSz, SomUtils::MatD &Tijs, Eigen::MatrixXi &edges)
    {
        sz_ = somSz;
        Tijs_ = Tijs;
        edges_ = edges;
        numEdges_ = Tijs_.cols();
        fullSz_ = sz_.d_ * sz_.p_ * sz_.n_ + sz_.p_ * sz_.n_;

        Rgt_.resize(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
        Tgt_ = SomUtils::MatD::Zero(sz_.d_, sz_.n_);

        src_ = 0; // TODO: src_ VS src (for sure in globalize, maybe also in other places)

        costCurr_ = 1e+10;
    }

    SsomProblem::SsomProblem(const SomUtils::SomSize somSz, const SomUtils::MatD &Tijs, const Eigen::MatrixXi &edges)
    {
        sz_ = somSz;
        Tijs_ = Tijs;
        edges_ = edges;
        numEdges_ = Tijs_.cols();
        fullSz_ = sz_.d_ * sz_.p_ * sz_.n_ + sz_.p_ * sz_.n_;

        Rgt_.resize(sz_.n_, SomUtils::MatD::Zero(sz_.d_, sz_.d_));
        Tgt_ = SomUtils::MatD::Zero(sz_.d_, sz_.n_);

        src_ = 0;

        costCurr_ = 1e+10;
    }

    SsomProblem::~SsomProblem() {};

    realdp SsomProblem::f(const Variable &x) const
    {
        SomUtils::MatD xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);
        // ROFL_VAR1(x);

        realdp cost = costEigenVec(xEigen);

        ROFL_VAR1(cost);
        // ROFL_ASSERT(!std::isnan(corr));

        // Vector *resultEgrad;
        // *resultEgrad = Domain->RandominManifold();
        // EucGrad(x, resultEgrad);
        // x.AddToFields("EGrad", *resultEgrad); //x should have been const?? maybe only its reference

        return cost; // checked -> the - here should be OK
    };

    double SsomProblem::costEigenVecSEdN(const SomUtils::MatD &xEigen) const
    {
        double cost = 0.0f;
        for (int e = 0; e < numEdges_; ++e)
        {
            SomUtils::MatD Ri(SomUtils::MatD::Zero(sz_.d_, sz_.d_));
            SomUtils::MatD Ti(SomUtils::MatD::Zero(sz_.d_, 1));
            SomUtils::MatD Tj(SomUtils::MatD::Zero(sz_.d_, 1));

            SomUtils::VecD tij(SomUtils::VecD::Zero(sz_.d_));
            tij = Tijs_.col(e);

            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1
            getRiSEdN(xEigen, Ri, i);
            getRiSEdN(xEigen, Ti, i);
            getRiSEdN(xEigen, Tj, j);

            // ROFL_VAR3(i, j, e);
            // ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

            double costE = (Ri * tij - Tj + Ti).norm(); // TODO: use squaredNorm() here directly
            // ROFL_VAR1(costE);

            cost += costE * costE;
        }
        return cost;
    }

    double SsomProblem::costEigen(const SomUtils::VecMatD &Reigen, const SomUtils::MatD &Teigen) const
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
            tij = Tijs_.col(e);

            // ROFL_VAR3(i, j, e);
            // ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

            double costEsq = (Ri * tij - Tj + Ti).squaredNorm(); // TODO: use squaredNorm() here directly
            // ROFL_VAR1(costE);

            cost += costEsq;
        }
        return cost;
    }

    double SsomProblem::costEigenVec(const SomUtils::MatD &xEigen) const
    {
        double cost = 0.0f;
        for (int e = 0; e < numEdges_; ++e)
        {
            SomUtils::MatD Ri(SomUtils::MatD::Zero(sz_.p_, sz_.d_));
            SomUtils::MatD Ti(SomUtils::MatD::Zero(sz_.p_, 1));
            SomUtils::MatD Tj(SomUtils::MatD::Zero(sz_.p_, 1));

            SomUtils::VecD tij(SomUtils::VecD::Zero(sz_.d_));
            tij = Tijs_.col(e);

            int i = edges_(e, 0) - 1; // !! -1
            int j = edges_(e, 1) - 1; // !! -1
            getRi(xEigen, Ri, i);
            getTi(xEigen, Ti, i);
            getTi(xEigen, Tj, j);

            // ROFL_VAR3(i, j, e);
            // ROFL_VAR4(Ri, tij.transpose(), Ti.transpose(), Tj.transpose());

            double costE = (Ri * tij - Tj + Ti).norm(); // TODO: use squaredNorm() here directly
            // ROFL_VAR1(costE);

            cost += costE * costE;
        }
        return cost;
    }

    // Vector &SsomProblem::EucGrad(const Variable &x, Vector *result) const
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

    //     SomUtils::MatD xEig(fullSz_, 1);
    //     RoptToEig(x, xEig);

    //     SomUtils::VecMatD R(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
    //     getRotations(xEig, R);

    //     SomUtils::MatD T(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
    //     getTranslations(xEig, T);

    //     // P = zeros(nrs, d*N);
    //     SomUtils::MatD P(SomUtils::MatD::Zero(sz_.p_, sz_.d_ * sz_.n_));
    //     double frct = 0.0;

    //     // LR = zeros(N,N);
    //     // PR = zeros(N,nrs);
    //     // BR_const = zeros(d,d);
    //     SomUtils::MatD Lr(SomUtils::MatD::Zero(sz_.n_, sz_.n_));
    //     SomUtils::MatD Pr(SomUtils::MatD::Zero(sz_.n_, sz_.p_));
    //     SomUtils::MatD Br(SomUtils::MatD::Zero(sz_.d_, sz_.d_));

    //     makePfrct(T, P, frct);
    //     makeLrPrBr(R, Lr, Pr, Br);

    //     SomUtils::VecMatD egR(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
    //     egradR(P, egR);

    //     SomUtils::MatD egT(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
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

    Vector &SsomProblem::RieGrad(const Variable &x, Vector *result) const
    {
        // result->NewMemoryOnWrite();
        // result->GetElement(0) = x.Field("B1x1D1");
        // result->GetElement(1) = x.Field("B2x2D2");
        // result->GetElement(2) = x.Field("B3x3D3");
        // Domain->ScalarTimesVector(x, 2, *result, result);

        // result->Print("RieGrad: printing result at start of function (should be empty)");

        SomUtils::MatD xEig(fullSz_, 1);
        RoptToEig(x, xEig);

        SomUtils::VecMatD R(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
        getRotations(xEig, R);

        SomUtils::MatD T(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
        getTranslations(xEig, T);

        // P = zeros(nrs, d*N);
        SomUtils::MatD P(SomUtils::MatD::Zero(sz_.p_, sz_.d_ * sz_.n_));
        double frct = 0.0;

        // LR = zeros(N,N);
        // PR = zeros(N,nrs);
        // BR_const = zeros(d,d);
        SomUtils::MatD Lr(SomUtils::MatD::Zero(sz_.n_, sz_.n_));
        SomUtils::MatD Pr(SomUtils::MatD::Zero(sz_.n_, sz_.p_));
        SomUtils::MatD Br(SomUtils::MatD::Zero(sz_.d_, sz_.d_));

        makePfrct(T, P, frct);
        // ROFL_VAR2(P, frct);
        makeLrPrBr(R, Lr, Pr, Br);
        // ROFL_VAR3(Lr, Pr, Br);

        SomUtils::VecMatD rgR(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
        // ROFL_VAR3(xEig.rows(), R.size(), T.size());
        // ROFL_VAR3(R.size(), P.size(), rgR.size());
        rgradR(R, P, rgR);
        // ROFL_VAR5(rgR[0],rgR[1],rgR[2],rgR[3],rgR[4]);

        SomUtils::MatD rgT(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
        rgradT(T, Lr, Pr, rgT);
        // ROFL_VAR1(rgT);

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
        // result->GetElement(gElemIdx).Print("grad Ti after assignment");

        // ROFL_VAR2("\n", rgT);
        // ROFL_VAR1(gElemIdx);
        // result->GetElement(gElemIdx).Print();

        // result->NewMemoryOnWrite();
        // result->SetToZeros();
        // *result = Groptlib;

        // result->Print("printing final result");

        // ROFL_ASSERT(0);

        return *result;
    };

    // Vector &SsomProblem::Grad(const Variable &x, Vector *result) const
    // {
    //     // result->NewMemoryOnWrite();
    //     // result->GetElement(0) = x.Field("B1x1D1");
    //     // result->GetElement(1) = x.Field("B2x2D2");
    //     // result->GetElement(2) = x.Field("B3x3D3");
    //     // Domain->ScalarTimesVector(x, 2, *result, result);

    //     result->Print("Grad: printing it before anything happens result");

    //     SomUtils::MatD xEig(fullSz_, 1);
    //     RoptToEig(x, xEig);

    //     SomUtils::VecMatD R(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
    //     getRotations(xEig, R);

    //     SomUtils::MatD T(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
    //     getTranslations(xEig, T);

    //     // P = zeros(nrs, d*N);
    //     SomUtils::MatD P(SomUtils::MatD::Zero(sz_.p_, sz_.d_ * sz_.n_));
    //     double frct = 0.0;

    //     // LR = zeros(N,N);
    //     // PR = zeros(N,nrs);
    //     // BR_const = zeros(d,d);
    //     SomUtils::MatD Lr(SomUtils::MatD::Zero(sz_.n_, sz_.n_));
    //     SomUtils::MatD Pr(SomUtils::MatD::Zero(sz_.n_, sz_.p_));
    //     SomUtils::MatD Br(SomUtils::MatD::Zero(sz_.d_, sz_.d_));

    //     makePfrct(T, P, frct);
    //     makeLrPrBr(R, Lr, Pr, Br);

    //     SomUtils::VecMatD egR(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
    //     egradR(P, egR);

    //     SomUtils::MatD rgT(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
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

    void SsomProblem::egradR(const SomUtils::MatD &P, SomUtils::VecMatD &egR) const
    {
        SomUtils::unStackH(P, egR, sz_.d_);
    }

    void SsomProblem::rgradR(const SomUtils::VecMatD &R, const SomUtils::MatD &P, SomUtils::VecMatD &rgR) const
    {
        SomUtils::VecMatD egR;
        egradR(P, egR);

        SomUtils::stiefelTangentProj(R, egR, rgR);
    }

    void SsomProblem::egradT(const SomUtils::MatD &T, const SomUtils::MatD &Lr, const SomUtils::MatD &Pr, SomUtils::MatD &egT) const
    {
        // ROFL_VAR7(egT, T.rows(), T.cols(), Lr.rows(), Lr.cols(), Pr.rows(), Pr.cols());
        egT = T * (Lr + Lr.transpose()) + Pr.transpose();
    }

    void SsomProblem::rgradT(const SomUtils::MatD &T, const SomUtils::MatD &Lr, const SomUtils::MatD &Pr, SomUtils::MatD &egT) const
    {
        egradT(T, Lr, Pr, egT);
    }

    // Vector &SsomProblem::EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const
    // {
    //     // TODO: implement

    //     result->NewMemoryOnWrite();

    //     // h.R = rsom_rhess_rot_stiefel(R, Rdot, problem_structs) + hrt;
    //     // h.T = rsom_rhess_transl_stiefel(T, Tdot, problem_structs) + htr;

    //     result->Print("EucHessianEta: printing it just after NewMemoryOnWrite()");

    //     return *result;
    // };

    void SsomProblem::hessGenprocEigen(const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR, const SomUtils::MatD &xT, const SomUtils::MatD &uT, SomUtils::VecMatD &rhR, SomUtils::MatD &rhT) const
    {
        // P = zeros(nrs, d*N);
        int staircaseStep = xT.rows();
        SomUtils::MatD P(SomUtils::MatD::Zero(staircaseStep, sz_.d_ * sz_.n_));
        double frct = 0.0;

        // LR = zeros(N,N);
        // PR = zeros(N,nrs);
        // BR_const = zeros(d,d);
        SomUtils::MatD Lr(SomUtils::MatD::Zero(sz_.n_, sz_.n_));
        SomUtils::MatD Pr(SomUtils::MatD::Zero(sz_.n_, staircaseStep));
        SomUtils::MatD Br(SomUtils::MatD::Zero(sz_.d_, sz_.d_));

        // SomUtils::VecMatD rhR(sz_.n_, SomUtils::MatD::Zero(staircaseStep, sz_.d_));
        // SomUtils::MatD rhT(SomUtils::MatD::Zero(staircaseStep, sz_.n_));
        makePfrct(xT, P, frct);
        // ROFL_VAR2(P, frct);
        makeLrPrBr(xR, Lr, Pr, Br);
        // ROFL_VAR3(Lr, Pr, Br);

        SomUtils::VecMatD rgR(sz_.n_, SomUtils::MatD::Zero(staircaseStep, sz_.d_));
        rgradR(xR, P, rgR);

        SomUtils::MatD rgT(SomUtils::MatD::Zero(staircaseStep, sz_.n_));
        rgradT(xT, Lr, Pr, rgT);

        // result->NewMemoryOnWrite();
        // result = Domain->RandomInManifold();

        // compute Hrr, Htt
        SomUtils::VecMatD hrr(sz_.n_, SomUtils::MatD::Zero(staircaseStep, sz_.d_));
        computeHrr(xR, uR, P, hrr);
        SomUtils::MatD htt(SomUtils::MatD::Zero(staircaseStep, sz_.n_));
        computeHtt(uT, Lr, htt);

        // h.R = rsom_rhess_rot_stiefel(R, Rdot, problem_structs) + hrt;
        // h.T = rsom_rhess_transl_stiefel(T, Tdot, problem_structs) + htr;
        SomUtils::VecMatD hrt(sz_.n_, SomUtils::MatD::Zero(staircaseStep, sz_.d_));
        computeHrt(xR, uT, hrt);
        SomUtils::MatD htr(SomUtils::MatD::Zero(staircaseStep, sz_.n_));
        computeHtr(uR, htr);

        // ROFL_VAR2(hrr[0], hrt[0]);
        // ROFL_VAR2(hrr[1], hrt[1]);
        // ROFL_VAR2(hrr[2], hrt[2]);
        // ROFL_VAR2(hrr[3], hrt[3]);
        // ROFL_VAR2(hrr[4], hrt[4]);
        // ROFL_VAR2(htr, htt);

        // TODO: use less memory

        // sum diag component with antidiag one
        std::transform(hrr.begin(), hrr.end(), hrt.begin(), rhR.begin(), std::plus<SomUtils::MatD>());
        // ROFL_VAR5(rhR[0], rhR[1], rhR[2], rhR[3], rhR[4]);
        rhT = htt + htr;
        // ROFL_VAR1(rhT);
    }

    void SsomProblem::hessGenprocEigenShifted(
        const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR, const SomUtils::MatD &xT, const SomUtils::MatD &uT,
        double mu, SomUtils::VecMatD &rhR, SomUtils::MatD &rhT) const
    {
        hessGenprocEigen(xR, uR, xT, uT, rhR, rhT);
        for (int i = 0; i < sz_.n_; ++i)
        {
            rhR[i] -= mu * uR[i];
        };
        // ROFL_VAR5(rhR[0], rhR[1], rhR[2], rhR[3], rhR[4]);

        // rhT = htt + htr;
        rhT -= mu * uT; // shift!
        // ROFL_VAR1(rhT);
    }

    Vector &SsomProblem::RieHessianEta(const Variable &x, const Vector &etax, Vector *result) const
    {
        // TODO: implement

        result->NewMemoryOnWrite();

        // result->Print("RieHessianEta: printing it just after NewMemoryOnWrite()");

        // x.Print("x inside RieHessianEta");
        // etax.Print("etax inside RieHessianEta");

        SomUtils::MatD xEig(fullSz_, 1);
        RoptToEig(x, xEig);

        SomUtils::MatD xEtaEig(fullSz_, 1);
        RoptToEig(etax, xEtaEig);

        SomUtils::VecMatD R(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
        getRotations(xEig, R);

        SomUtils::VecMatD uR(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
        getRotations(xEtaEig, uR);

        SomUtils::MatD T(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
        getTranslations(xEig, T);

        SomUtils::MatD uT(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
        getTranslations(xEtaEig, uT);

        // ROFL_VAR2(R[0], uR[0]);
        // ROFL_VAR2(T, uT);

        SomUtils::VecMatD rhR(sz_.n_, SomUtils::MatD::Zero(sz_.p_, sz_.d_));
        SomUtils::MatD rhT(SomUtils::MatD::Zero(sz_.p_, sz_.n_));
        hessGenprocEigen(R, uR, T, uT, rhR, rhT);

        int rotSz = getRotSz();
        int translSz = getTranslSz();

        int gElemIdx = 0;
        // fill result with computed gradient values : R
        for (int i = 0; i < sz_.n_; ++i)
        {
            // ROFL_VAR1(gElemIdx);
            // ROFL_VAR2("\n", rhR[gElemIdx]);
            // result->GetElement(gElemIdx).SetToIdentity(); // Ri
            // result->GetElement(gElemIdx).Print("Ri before assignment");

            Vector rhRiVec(sz_.p_, sz_.d_);
            // rhRiVec.Initialize();
            realdp *GroptlibWriteArray = rhRiVec.ObtainWriteEntireData();
            for (int j = 0; j < rotSz; ++j)
            {
                // ROFL_VAR2(i, j);
                // rhRiVec.Print("rhRiVec before assignment");

                // ROFL_VAR1(rhRiVec.GetElement(j, 0));

                GroptlibWriteArray[j] = rhR[i].reshaped(sz_.d_ * sz_.p_, 1)(j);

                // ROFL_VAR1("");
                // rhRiVec.Print("rhRiVec after assignment");
            }
            rhRiVec.CopyTo(result->GetElement(gElemIdx));
            // result->GetElement(gElemIdx).Print("Riem. grad Ri after assignment");
            gElemIdx++;
        }

        // fill result with computed gradient values : T

        Vector rhTiVec(sz_.p_, sz_.n_);
        realdp *GroptlibWriteArray = rhTiVec.ObtainWriteEntireData();
        for (int j = 0; j < sz_.p_ * sz_.n_; ++j)
        {
            // rhTiVec.Print("rhTiVec before assignment");

            // ROFL_VAR1(rhRiVec.GetElement(j, 0));

            GroptlibWriteArray[j] = rhT.reshaped(sz_.n_ * sz_.p_, 1)(j);

            // ROFL_VAR1("");
            // rhTiVec.Print("rhTiVec after assignment");
        }
        rhTiVec.CopyTo(result->GetElement(gElemIdx));
        // result->GetElement(gElemIdx).Print("RieHess T after assignment");

        // result->Print("RieHessianEta: printing just before end of function");

        return *result;
    };

    void SsomProblem::RoptToEig(Vector x, SomUtils::MatD &xEigen) const
    {
        Vector xT = x.GetTranspose(); // Eigen ADV init is row-major!!

        int totSz = xEigen.rows();

        const realdp *xArr = xT.ObtainWriteEntireData();
        for (int i = 0; i < totSz; ++i)
            xEigen(i) = xArr[i];
    }

    void SsomProblem::getRi(const Variable &x, SomUtils::MatD &rOut, int i) const
    {
        SomUtils::MatD xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);

        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        int startId = i * rotSz;
        // int endId = (i+1) * rotSz;

        SomUtils::MatD rOutVec(SomUtils::MatD::Zero(rotSz, 1));
        rOutVec = xEigen.block(startId, 0, rotSz, 1);

        rOut = rOutVec.reshaped(sz_.p_, sz_.d_);
    }

    void SsomProblem::getRi(const SomUtils::MatD &xEig, SomUtils::MatD &rOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        int startId = i * rotSz;
        // int endId = (i+1) * rotSz;

        SomUtils::MatD rOutVec(SomUtils::MatD::Zero(rotSz, 1));
        rOutVec = xEig.block(startId, 0, rotSz, 1);

        rOut = rOutVec.reshaped(sz_.p_, sz_.d_);
    }

    void SsomProblem::getRiSEdN(const SomUtils::MatD &xEig, SomUtils::MatD &rOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = sz_.d_ * sz_.d_; // as a vector

        int startId = i * rotSz;
        // int endId = (i+1) * rotSz;

        SomUtils::MatD rOutVec(SomUtils::MatD::Zero(rotSz, 1));
        rOutVec = xEig.block(startId, 0, rotSz, 1);

        rOut = rOutVec.reshaped(sz_.d_, sz_.d_);
    }

    void SsomProblem::getRotations(const SomUtils::MatD &xEig, SomUtils::VecMatD &rOut) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // vectorized

        ROFL_ASSERT(rotSz * sz_.n_ + sz_.p_ * sz_.n_ == xEig.rows())

        // int endId = (i+1) * rotSz;

        std::for_each(rOut.begin(), rOut.end(), [](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
            x.setZero();
        });

        for (int i = 0; i < sz_.n_; ++i)
        {
            int startId = i * rotSz;
            SomUtils::MatD Ri(SomUtils::MatD::Zero(sz_.p_, sz_.d_));
            getRi(xEig, Ri, i); // TODO: this can probably be optimized better
            rOut[i] = Ri;
            // ROFL_VAR3(i, rOut[i].rows(), rOut[i].cols())
        }
    }

    void SsomProblem::getTi(const Variable &x, SomUtils::MatD &tOut, int i) const
    {
        SomUtils::MatD xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);

        // rOut already needs to have fixed size by here
        int rotSz = getRotSz();
        int translSz = getTranslSz();

        int startId = sz_.n_ * rotSz + i * translSz;
        // int endId = (i+1) * rotSz;

        // ROFL_VAR1(startId);

        tOut.setZero(); // TODO: remove later
        tOut = xEigen.block(startId, 0, translSz, 1);
    }

    void SsomProblem::getTi(const SomUtils::MatD &xEig, SomUtils::MatD &tOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz();
        int translSz = getTranslSz();

        int startId = sz_.n_ * rotSz + i * translSz;
        // int endId = (i+1) * rotSz;

        // ROFL_VAR1(startId);
        tOut.setZero(); // TODO: remove later
        tOut = xEig.block(startId, 0, translSz, 1);
    }

    void SsomProblem::getTiSEdN(const SomUtils::MatD &xEig, SomUtils::MatD &tOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = sz_.d_ * sz_.d_;
        int translSz = sz_.d_;

        int startId = sz_.n_ * rotSz + i * translSz;
        // int endId = (i+1) * rotSz;

        // ROFL_VAR1(startId);
        tOut.setZero(); // TODO: remove later
        tOut = xEig.block(startId, 0, translSz, 1);
    }

    void SsomProblem::getTranslations(const SomUtils::MatD &xEig, SomUtils::MatD &tOut) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        // int endId = (i+1) * rotSz;

        for (int i = 0; i < sz_.n_; ++i)
        {
            int startId = i * rotSz;
            SomUtils::MatD Ti(sz_.p_, 1);
            getTi(xEig, Ti, i); // TODO: this can probably be optimized better
            tOut.col(i) = Ti;
        }
    }

    void SsomProblem::makePfrct(const SomUtils::MatD &T, SomUtils::MatD &P, double &frct) const
    {
        // P = zeros(nrs, d*N);

        // TODO: add some size checks at least for P

        P.setZero();
        frct = 0.0;

        for (int e = 0; e < numEdges_; ++e)
        {
            int ii = edges_(e, 0) - 1;
            int jj = edges_(e, 1) - 1;

            SomUtils::MatD T_j = T.col(jj);
            SomUtils::MatD T_i = T.col(ii);
            SomUtils::MatD tij = Tijs_.col(e);
            SomUtils::MatD Pe = 2 * (T_i * tij.transpose() - T_j * tij.transpose());
            P.block(0, ii * sz_.d_, T.rows(), sz_.d_) += Pe;

            SomUtils::MatD c = T_i * T_i.transpose() + T_j * T_j.transpose() - T_i * T_j.transpose() - T_j * T_i.transpose();
            SomUtils::MatD d = tij * tij.transpose();
            frct += c.trace() + d.trace();
        }
    }

    void SsomProblem::makeLrPrBr(const SomUtils::VecMatD &R, SomUtils::MatD &Lr, SomUtils::MatD &Pr, SomUtils::MatD &Br) const
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
            SomUtils::MatD Ri = R[ii];

            SomUtils::MatD bij(SomUtils::MatD::Zero(sz_.n_, 1));

            bij(ii) = 1;
            bij(jj) = -1;

            SomUtils::MatD tij = Tijs_.col(e);

            Lr += bij * bij.transpose();

            Pr += 2 * bij * tij.transpose() * Ri.transpose();

            Br += tij * tij.transpose();
        }
    }

    void SsomProblem::computeHrr(const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR, const SomUtils::MatD &P, SomUtils::VecMatD &hRR) const
    {
        int staircaseStep = xR[0].rows();
        SomUtils::VecMatD egR(sz_.n_, SomUtils::MatD::Zero(staircaseStep, sz_.d_));
        egradR(P, egR);
        // ROFL_VAR2(P, egR[0]);

        SomUtils::VecMatD term_1(sz_.n_, SomUtils::MatD::Zero(staircaseStep, sz_.d_));
        SomUtils::VecMatD term_2(sz_.n_, SomUtils::MatD::Zero(staircaseStep, sz_.d_));
        SomUtils::VecMatD DGf(sz_.n_, SomUtils::MatD::Zero(staircaseStep, sz_.d_));

        // TODO: improve efficiency of this
        //  Original code was:
        //  G = rsom_egrad_rot_stiefel(x, problem);
        //  term_1 = multiprod(u, 0.5*multiprod(multitransp(x), G) + 0.5*multiprod(multitransp(G), x));
        //  term_2 = multiprod(x, 0.5*multiprod(multitransp(u), G) + 0.5*multiprod(multitransp(G), u));
        //  DGf = - term_1 - term_2;
        //  h = stiefel_tangentProj(x, DGf);

        for (int i = 0; i < sz_.n_; ++i)
        {
            term_1[i] = uR[i] * (0.5 * (xR[i].transpose() * egR[i]) + 0.5 * (egR[i].transpose() * xR[i]));
            term_2[i] = xR[i] * (0.5 * (uR[i].transpose() * egR[i]) + 0.5 * (egR[i].transpose() * uR[i]));
            DGf[i] = -term_1[i] - term_2[i];
            // ROFL_VAR2(i, egR[i]);
            // ROFL_VAR2(xR[i], uR[i]);
            // ROFL_VAR3(term_1[i], term_2[i], DGf[i]);
        }

        std::for_each(hRR.begin(), hRR.end(), [](SomUtils::MatD &x) { //^^^ take argument by reference
            x.setZero();
        });
        SomUtils::stiefelTangentProj(xR, DGf, hRR);
        // for (int i = 0; i < sz_.n_; ++i)
        // {
        //     ROFL_VAR2(i, hRR[i]);
        // }
    }

    void SsomProblem::computeHtt(const SomUtils::MatD &uT, const SomUtils::MatD &LR, SomUtils::MatD &hTT) const
    {
        hTT = uT * (LR.transpose() + LR);
    }

    void SsomProblem::computeHrt(const SomUtils::VecMatD &xR, const SomUtils::MatD uT, SomUtils::VecMatD &hrt) const
    {
        int staircaseStep = xR[0].rows();
        SomUtils::VecMatD W(sz_.n_, SomUtils::MatD::Zero(staircaseStep, sz_.d_));

        for (int e = 0; e < numEdges_; ++e)
        {
            int ii = edges_(e, 0) - 1;
            int jj = edges_(e, 1) - 1;
            SomUtils::MatD uTi(SomUtils::MatD::Zero(staircaseStep, 1));
            uTi = uT.col(ii);
            SomUtils::MatD uTj(SomUtils::MatD::Zero(staircaseStep, 1));
            uTj = uT.col(jj);

            SomUtils::MatD tij = Tijs_.col(e);
            SomUtils::MatD wij = 2 * (uTi - uTj) * tij.transpose();
            W[ii] += wij;
        }

        std::for_each(hrt.begin(), hrt.end(), [](SomUtils::MatD &x) { //^^^ take argument by reference
            x.setZero();
        });

        SomUtils::stiefelTangentProj(xR, W, hrt);
        // ROFL_VAR1(uT);
        // for (int i = 0; i < sz_.n_; ++i)
        // {
        //     ROFL_VAR2(i, hrt[i]);
        // }
    }

    void SsomProblem::computeHtr(const SomUtils::VecMatD &uR, SomUtils::MatD &htr) const
    {
        htr.setZero();

        for (int e = 0; e < numEdges_; ++e)
        {
            int ii = edges_(e, 0) - 1;
            int jj = edges_(e, 1) - 1;
            SomUtils::MatD uRi = uR[ii];
            // uRj = uR[jj];

            SomUtils::MatD bij(SomUtils::MatD::Zero(sz_.n_, 1));
            bij(ii, 0) = 1;
            bij(jj, 0) = -1;

            SomUtils::MatD tij = Tijs_.col(e);
            SomUtils::MatD wij = 2 * bij * tij.transpose() * uRi.transpose();
            htr += wij.transpose();
            // ROFL_VAR5(bij.transpose(), tij.transpose(), uRi, wij, htr);
        }
    }

    int SsomProblem::getRotSz() const
    {
        return sz_.d_ * sz_.p_;
    }

    int SsomProblem::getTranslSz() const
    {
        return sz_.p_;
    }

    void SsomProblem::vectorizeR(const SomUtils::VecMatD &R, SomUtils::MatD &RvecOut) const
    {
        // int fullRotsSz = sz_.p_ * sz_.d_ * sz_.n_;

        // for (int i=0; i<fullRotsSz; ++i) {
        // }

        int n = R.size();

        int fullIdx = 0;
        for (int i = 0; i < n; ++i)
        {
            int ric = R[i].cols();
            int rir = R[i].rows();
            for (int j = 0; j < ric; ++j)
            {
                for (int k = 0; k < rir; ++k)
                {
                    RvecOut(fullIdx, 0) = R[i](k, j);
                    fullIdx++;
                    // ROFL_VAR4(i, j, k, fullIdx);
                }
            }
        }
        ROFL_ASSERT(fullIdx == RvecOut.rows())
    }

    void SsomProblem::vectorizeRT(const SomUtils::VecMatD &R, const SomUtils::MatD &T, SomUtils::MatD &XvecOut) const
    {
        // int fullRotsSz = sz_.p_ * sz_.d_ * sz_.n_;

        // for (int i=0; i<fullRotsSz; ++i) {
        // }

        int n = R.size();

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

    void SsomProblem::RoptToEigStiefel(Vector x, SomUtils::MatD &xEigen) const
    {
        Vector xT = x.GetTranspose(); // Eigen ADV init is row-major!!

        int totSz = xEigen.rows(); // xEigen is supposed to be a vectorized matrix

        const realdp *xArr = xT.ObtainWriteEntireData();
        for (int i = 0; i < totSz; ++i)
            xEigen(i) = xArr[i];
    }

    void SsomProblem::makeAdjMatFromEdges(Eigen::MatrixXi &adjMat) const
    {
        ROFL_ASSERT(adjMat.rows() == sz_.n_)
        ROFL_ASSERT(adjMat.cols() == sz_.n_)

        adjMat.setZero();
        for (int k = 0; k < numEdges_; ++k)
        {
            int ii = edges_(k, 0) - 1;
            int jj = edges_(k, 1) - 1;
            adjMat(ii, jj) = 1;
        }
    }

    void SsomProblem::setGtR(const SomUtils::VecMatD &R)
    {
        Rgt_ = R;
    }

    void SsomProblem::setGtT(const SomUtils::MatD &T)
    {
        Tgt_ = T;
    }

    void SsomProblem::setGt(const SomUtils::VecMatD &R, const SomUtils::MatD &T)
    {
        Rgt_ = R;
        Tgt_ = T;
    }

    bool SsomProblem::getRsRecoverySuccess() const
    {
        return rsRecoverySuccess_;
    }

    void SsomProblem::setCostCurr(double cc)
    {
        costCurr_ = cc;
    }

    void SsomProblem::makeHmat(const SomUtils::MatD &XvecNext, const SomUtils::SomSize &szNext, SomUtils::MatD &Hmat) const
    {
        int staircaseStepLevel = szNext.p_; // TODO: it can maybe be deducted fron XvecNext.size()?

        ROFL_VAR1(staircaseStepLevel)

        ROPTLIB::SsomProblem ProbNextLocal(szNext, Tijs_, edges_);

        SomUtils::VecMatD xRi(szNext.n_, SomUtils::MatD::Zero(staircaseStepLevel, szNext.d_));
        SomUtils::MatD xTi(SomUtils::MatD::Zero(staircaseStepLevel, szNext.n_));
        ROFL_VAR1("Calling getRotations()")
        ProbNextLocal.getRotations(XvecNext, xRi);
        ProbNextLocal.getTranslations(XvecNext, xTi);

        int vecsz = XvecNext.rows();
        ROFL_VAR1(vecsz)

        for (int i = 0; i < vecsz; ++i)
        {
            SomUtils::MatD eI(SomUtils::MatD::Zero(vecsz, 1));
            eI(i) = 1;
            SomUtils::VecMatD uRi(szNext.n_, SomUtils::MatD::Zero(staircaseStepLevel, szNext.d_));
            SomUtils::MatD uTi(SomUtils::MatD::Zero(staircaseStepLevel, szNext.n_));

            ROFL_VAR1("Calling getRotations()")
            ProbNextLocal.getRotations(eI, uRi);
            ProbNextLocal.getTranslations(eI, uTi);

            SomUtils::VecMatD rhrI(szNext.n_, SomUtils::MatD::Zero(staircaseStepLevel, szNext.d_));
            SomUtils::MatD rhtI(SomUtils::MatD::Zero(staircaseStepLevel, szNext.n_));

            ROFL_VAR2(i, vecsz)
            ROFL_VAR1(Tijs_)
            ROFL_VAR3(szNext.n_, xRi[szNext.n_ - 1].rows(), xRi[szNext.n_ - 1].cols())
            ROFL_VAR2(uRi[0], uTi)
            ProbNextLocal.hessGenprocEigen(xRi, uRi, xTi, uTi, rhrI, rhtI);
            ROFL_VAR2(rhrI[0], rhtI)

            SomUtils::MatD rhVecI(SomUtils::MatD::Zero(vecsz, 1));
            ProbNextLocal.vectorizeRT(rhrI, rhtI, rhVecI);
            Hmat.col(i) = rhVecI;
            ROFL_VAR2(i, rhVecI.transpose())
        }
    }

    double runSsom(ROPTLIB::SsomProblem &Prob,
                   const ROPTLIB::Vector &startX,
                   int src,
                   SomUtils::VecMatD &Rout,
                   SomUtils::MatD &Tout,
                   int &staircaseStepIdx)
    {
        ROFL_VAR1("Start of runRsomRS()")
        ROFL_VAR1(Prob.costEigen(Prob.Rgt_, Prob.Tgt_));

        // output the parameters of the manifold of domain
        ROPTLIB::RTRNewton *RTRNewtonSolver = new ROPTLIB::RTRNewton(&Prob, &startX); // USE INITGUESS HERE!
        RTRNewtonSolver->Verbose = ROPTLIB::ITERRESULT;
        // RTRNewtonSolver->Max_Iteration = 500;
        // RTRNewtonSolver->Max_Inner_Iter = 500;
        // ROPTLIB::PARAMSMAP solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 10)};
        // RTRNewtonSolver->SetParams(solverParams);
        RTRNewtonSolver->CheckParams();

        // % Solve.
        // [x, xcost, info, options] = trustregions(problem);
        RTRNewtonSolver->Run();
        // Numerically check gradient consistency (optional).
        auto Xopt = RTRNewtonSolver->GetXopt();
        auto XoptCost = RTRNewtonSolver->Getfinalfun();

        // Prob.CheckGradHessian(Xopt);

        // std::cout << "Prob.GetUseGrad() " << Prob.GetUseGrad() << std::endl;
        // std::cout << "Prob.GetUseHess() " << Prob.GetUseHess() << std::endl;
        // std::cout << "Prob.GetNumGradHess() " << Prob.GetNumGradHess() << std::endl;

        // Outputs
        Xopt.Print("Xopt");
        std::cout << "XoptCost " << XoptCost << std::endl; // x cost

        delete RTRNewtonSolver;

        // RS
        int d = Prob.sz_.d_;
        int n = Prob.sz_.n_;
        int r0 = d + 1;
        int e = Prob.sz_.e_;

        integer numoftypes = 3; // 2 i.e. (3D) Stiefel + Euclidean
        integer numofmani1 = n; // num of Stiefel manifolds
        integer numofmani2 = 1;
        integer numofmani3 = 1;

        ROPTLIB::Stiefel mani1(d, d);
        mani1.ChooseParamsSet2();
        ROPTLIB::Euclidean mani2(d, n);
        ROPTLIB::Euclidean mani3(e);
        ROPTLIB::ProductManifold ProdManiSsom(numoftypes, 
            &mani1, numofmani1, &mani2, numofmani2, &mani3, numofmani3);

        SomUtils::MatD XoptEigVec(SomUtils::MatD::Zero(d * d * n + d * n + e, 1));
        Prob.RoptToEig(Xopt, XoptEigVec);
        ROFL_VAR1(XoptEigVec.transpose())

        double costLast = XoptCost;
        auto ProbPrev = Prob;
        // int staircaseStepIdx;
        SomUtils::VecMatD RmanoptOutEig(n, SomUtils::MatD::Zero(d, d));
        SomUtils::MatD TmanoptOutEig(SomUtils::MatD::Zero(d, n));

        auto rGt = Prob.Rgt_;
        auto tGt = Prob.Tgt_;

        // !! COMMENTING OUT THE STAIRCASE LOOP AND THE RECOVERY PROCEDURE FOR NOW
        // for (staircaseStepIdx = r0; staircaseStepIdx <= d * n + 1; ++staircaseStepIdx)
        // {
        //     ROFL_VAR1(staircaseStepIdx)
        //     ROFL_VAR1(costLast)

        //     SomUtils::VecMatD R(n, SomUtils::MatD::Zero(staircaseStepIdx - 1, d));
        //     SomUtils::MatD T(SomUtils::MatD::Zero(staircaseStepIdx - 1, n));
        //     {
        //         SomUtils::SomSize somSzScope(staircaseStepIdx - 1, d, n); // TODO: improve getRotations() and getTranslations() and avoid local scope
        //         ROPTLIB::SsomProblem ProbScope(somSzScope, Prob.Tijs_, Prob.edges_);
        //         ROFL_VAR1("Calling getRotations()")
        //         ProbScope.getRotations(XoptEigVec, R);
        //         ProbScope.getTranslations(XoptEigVec, T);
        //     }

        //     // SomUtils::VecMatD Rnext(n, SomUtils::MatD::Zero(staircaseStepIdx, d));
        //     // SomUtils::MatD Tnext(SomUtils::MatD::Zero(staircaseStepIdx, n));

        //     // Prob.catZeroRow3dArray(R, Rnext);
        //     // Prob.catZeroRow(T, Tnext);

        //     SomUtils::SomSize somSzNext(staircaseStepIdx, d, n);
        //     ROPTLIB::SsomProblem ProbNext(somSzNext, Prob.Tijs_, Prob.edges_);

        //     ProbNext.setGt(rGt, tGt);
        //     ROFL_VAR1(ProbNext.costEigen(ProbNext.Rgt_, ProbNext.Tgt_));
        //     ROFL_VAR1(ProbNext.costEigen(R, T));

        //     ROPTLIB::Stiefel mani1next(somSzNext.p_, somSzNext.d_);
        //     mani1next.ChooseParamsSet2();
        //     ROPTLIB::Euclidean mani2next(somSzNext.p_, somSzNext.n_);
        //     ROPTLIB::ProductManifold ProdManiNext(numoftypes, &mani1next, numofmani1, &mani2next, numofmani2);
        //     ROPTLIB::Vector Y0;
        //     double lambda;
        //     SomUtils::VecMatD vLambdaR(n, SomUtils::MatD::Zero(somSzNext.p_, somSzNext.d_));
        //     SomUtils::MatD vLambdaT(SomUtils::MatD::Zero(somSzNext.p_, somSzNext.n_));
        //     ProbPrev.setCostCurr(costLast);
        //     for (auto &Rm : R)
        //         ROFL_VAR1(ProbPrev.checkIsOnStiefel(Rm))
        //     ROFL_VAR1("Calling ProbPrev.rsomEscapeHessianGenprocEigen()")
        //     ProbPrev.rsomEscapeHessianGenprocEigen(R, T, Y0, lambda, vLambdaR, vLambdaT);

        //     if (lambda > -1e-2)
        //     {
        //         ROFL_VAR2(lambda, "R, T eigenvals > 0: exiting staircase")
        //         // staircaseStepSkipped = 0;
        //         RmanoptOutEig = R;
        //         TmanoptOutEig = T;
        //         break;
        //     }

        //     double costNewStart = ProbNext.f(Y0);
        //     ROFL_VAR1(costNewStart)

        //     // Run next step of staircase with found initial guess

        //     Y0.Print("Y0");

        //     // Set Prob params
        //     ProbNext.SetDomain(&ProdManiNext);
        //     ProbNext.SetUseGrad(true);
        //     ProbNext.SetUseHess(true);

        //     ROPTLIB::RTRNewton *RTRNewtonSolverNext = new ROPTLIB::RTRNewton(&ProbNext, &Y0); // USE INITGUESS HERE!
        //     RTRNewtonSolverNext->Verbose = ROPTLIB::ITERRESULT;
        //     // RTRNewtonSolverNext->Max_Iteration = 500;
        //     // RTRNewtonSolverNext->Max_Inner_Iter = 500;
        //     // ROPTLIB::PARAMSMAP solverParams = {std::pair<std::string, double>("Max_Inner_Iter", 10)};
        //     // RTRNewtonSolverNext->SetParams(solverParams);
        //     RTRNewtonSolverNext->CheckParams();

        //     // % Solve.
        //     // [x, xcost, info, options] = trustregions(problem);
        //     RTRNewtonSolverNext->Run();
        //     auto XoptNext = RTRNewtonSolverNext->GetXopt();
        //     realdp XoptNextCost = RTRNewtonSolverNext->Getfinalfun();
        //     // Numerically check gradient consistency (optional).
        //     ProbNext.CheckGradHessian(XoptNext);

        //     costLast = XoptNextCost;

        //     ROFL_VAR1(costLast)

        //     XoptEigVec.resize(staircaseStepIdx * d * n + staircaseStepIdx * n, 1);
        //     ProbNext.RoptToEig(XoptNext, XoptEigVec);
        //     XoptNext.Print("XoptNext");
        //     SomUtils::VecMatD XoptNextR(n, SomUtils::MatD::Zero(staircaseStepIdx, d));
        //     SomUtils::MatD XoptNextT(SomUtils::MatD::Zero(staircaseStepIdx, n));
        //     ROFL_VAR1("Calling getRotations()")
        //     ProbNext.getRotations(XoptEigVec, XoptNextR);
        //     ProbNext.getTranslations(XoptEigVec, XoptNextT);

        //     // std::cout << "Prob.GetUseGrad() " << Prob.GetUseGrad() << std::endl;
        //     // std::cout << "Prob.GetUseHess() " << Prob.GetUseHess() << std::endl;
        //     // std::cout << "Prob.GetNumGradHess() " << Prob.GetNumGradHess() << std::endl;

        //     // Outputs
        //     XoptNext.Print("XoptNext");
        //     std::cout << "XoptNextCost " << XoptNextCost << std::endl; // x cost

        //     ProbNext.RoptToEig(XoptNext, XoptEigVec);
        //     ROFL_VAR1(XoptEigVec.transpose())

        //     ProbPrev = ProbNext;

        //     // save output
        //     for (int i = 0; i < n; ++i)
        //     {
        //         RmanoptOutEig[i].resize(staircaseStepIdx, d);
        //     }
        //     TmanoptOutEig.resize(staircaseStepIdx, n);

        //     RmanoptOutEig = XoptNextR;
        //     TmanoptOutEig = XoptNextT;

        //     // Rank stopping condition
        //     ROFL_VAR1(staircaseStepIdx)
        //     SomUtils::MatD XoutRhSt(SomUtils::MatD::Zero(staircaseStepIdx, d * n));
        //     SomUtils::hstack(RmanoptOutEig, XoutRhSt);
        //     Eigen::FullPivLU<SomUtils::MatD> luDecomp(XoutRhSt);
        //     auto rank = luDecomp.rank();
        //     if (rank < staircaseStepIdx)
        //     {
        //         staircaseStepIdx++;
        //         ROFL_VAR1("Rank stopping condition reached -> Exiting RS");
        //         break;
        //     }

        //     delete RTRNewtonSolverNext;

        //     // break; // TODO: For now, just 1 RS step allowed -> remove it later after fixing linesearch
        // }

        // // Recovery procedure

        // ROFL_VAR1("Running recovery procedure")

        // // back to SE(d)^N
        // SomUtils::VecMatD Rrecovered(n, SomUtils::MatD::Zero(d, d));
        // SomUtils::MatD Trecovered(SomUtils::MatD::Zero(d, n));
        // bool recSEDNsuccess = ProbPrev.recoverySEdN(staircaseStepIdx,
        //                                             RmanoptOutEig, TmanoptOutEig,
        //                                             Rrecovered, Trecovered);

        // ROFL_VAR1("Printing R, T SE(d)^N")
        // for (auto &m : Rrecovered)
        //     ROFL_VAR1(m)
        // ROFL_VAR1(Trecovered)

        // ROFL_VAR1(recSEDNsuccess)

        // // globalize

        // ROFL_VAR1("Running globalization procedure")

        // Rout.resize(n, SomUtils::MatD::Zero(d, d));
        // Tout.resize(d, n);
        // bool globalRecoverySuccess = ProbPrev.globalize(src, Rrecovered, Trecovered,
        //                                                 Rout, Tout);
        // ROFL_VAR1(globalRecoverySuccess)

        return costLast;
    }

} // end of namespace ROPTLIB
