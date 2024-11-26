#include "problems_SampleSom.h"

namespace ROPTLIB
{
    SampleSomProblem::SampleSomProblem() {}

    SampleSomProblem::SampleSomProblem(SomUtils::SomSize somSz, SomUtils::MatD &Tijs, Eigen::MatrixXi &edges)
    {
        sz_ = somSz;
        Tijs_ = Tijs;
        edges_ = edges;
        numEdges_ = Tijs_.cols();
        fullSz_ = sz_.d_ * sz_.p_ * sz_.n_ + sz_.p_ * sz_.n_;
    }

    SampleSomProblem::SampleSomProblem(const SomUtils::SomSize somSz, const SomUtils::MatD &Tijs, const Eigen::MatrixXi &edges)
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
        SomUtils::MatD xEigen(fullSz_, 1);
        RoptToEig(x, xEigen);
        // ROFL_VAR1(x);

        realdp cost = costEigen(xEigen);

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

    Vector &SampleSomProblem::RieGrad(const Variable &x, Vector *result) const
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

    void SampleSomProblem::egradR(const SomUtils::MatD &P, SomUtils::VecMatD &egR) const
    {
        unStackH(P, egR, sz_.d_);
    }

    void SampleSomProblem::rgradR(const SomUtils::VecMatD &R, const SomUtils::MatD &P, SomUtils::VecMatD &rgR) const
    {
        SomUtils::VecMatD egR;
        egradR(P, egR);

        stiefelTangentProj(R, egR, rgR);
    }

    void SampleSomProblem::egradT(const SomUtils::MatD &T, const SomUtils::MatD &Lr, const SomUtils::MatD &Pr, SomUtils::MatD &egT) const
    {
        // ROFL_VAR7(egT, T.rows(), T.cols(), Lr.rows(), Lr.cols(), Pr.rows(), Pr.cols());
        egT = T * (Lr + Lr.transpose()) + Pr.transpose();
    }

    void SampleSomProblem::rgradT(const SomUtils::MatD &T, const SomUtils::MatD &Lr, const SomUtils::MatD &Pr, SomUtils::MatD &egT) const
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

    void SampleSomProblem::hessGenprocEigen(const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR, const SomUtils::MatD &xT, const SomUtils::MatD &uT, SomUtils::VecMatD &rhR, SomUtils::MatD &rhT) const
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

    void SampleSomProblem::hessGenprocEigenShifted(
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

    Vector &SampleSomProblem::RieHessianEta(const Variable &x, const Vector &etax, Vector *result) const
    {
        // TODO: implement

        result->NewMemoryOnWrite();

        // result->Print("RieHessianEta: printing it just after NewMemoryOnWrite()");

        x.Print("x inside RieHessianEta");
        etax.Print("etax inside RieHessianEta");

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
        result->GetElement(gElemIdx).Print("RieHess T after assignment");

        result->Print("RieHessianEta: printing just before end of function");

        return *result;
    };

    void SampleSomProblem::RoptToEig(Vector x, SomUtils::MatD &xEigen) const
    {
        Vector xT = x.GetTranspose(); // Eigen ADV init is row-major!!

        int totSz = xEigen.rows();

        const realdp *xArr = xT.ObtainWriteEntireData();
        for (int i = 0; i < totSz; ++i)
            xEigen(i) = xArr[i];
    }

    void SampleSomProblem::vstack(const SomUtils::VecMatD &in, SomUtils::MatD &out) const
    {
        ROFL_ASSERT(out.cols() == in[0].cols() && out.rows() == in[0].rows() * in.size());

        // out has same number of columns, whereas number of rows is the product of the size of the other 2 dimensions

        int rowJump = in[0].rows();
        for (int i = 0; i < in.size(); ++i)
        {
            out.block(rowJump * i, 0, rowJump, out.cols()) = in[i];
        }
    }

    void SampleSomProblem::hstack(const SomUtils::VecMatD &in, SomUtils::MatD &out) const
    {
        // out has same number of columns, whereas number of rows is the product of the size of the other 2 dimensions

        ROFL_ASSERT_VAR5(out.rows() == in[0].rows() && out.cols() == in[0].cols() * in.size(), out.rows(), out.cols(), in[0].rows(), in[0].cols(), in.size());

        int colJump = in[0].cols();
        for (int i = 0; i < in.size(); ++i)
        {
            // ROFL_VAR2(out.block(0, colJump * i, out.rows(), colJump), in[i]);
            out.block(0, colJump * i, out.rows(), colJump) = in[i];
        }
    }

    void SampleSomProblem::unStackV(const SomUtils::MatD &in, SomUtils::VecMatD &out, int rowsOut) const
    {
        int n = (int)in.rows() / rowsOut;

        int fixedSz = sz_.d_; // size that does not change in the 3D->2D transition (here, number of columns)

        // ROFL_VAR3(n, rowsOut, in.rows());
        ROFL_ASSERT(n * rowsOut == in.rows());

        out.clear();
        out.resize(n, SomUtils::MatD::Zero(rowsOut, in.cols()));

        for (int i = 0; i < n; ++i)
        {
            out[i] = in.block(rowsOut * i, 0, rowsOut, fixedSz);
            ROFL_ASSERT(out[i].cols() == in.cols())
        }
    }

    void SampleSomProblem::unStackH(const SomUtils::MatD &in, SomUtils::VecMatD &out, int colsOut) const
    {
        int n = (int)in.cols() / colsOut;

        int fixedSz = in.rows(); // size that does not change in the 3D->2D transition (here, number of rows)

        ROFL_ASSERT(n * colsOut == in.cols());

        out.clear();
        out.resize(n, SomUtils::MatD::Zero(in.rows(), colsOut));

        for (int i = 0; i < n; ++i)
        {
            out[i] = in.block(0, colsOut * i, fixedSz, colsOut);
            ROFL_ASSERT(out[i].rows() == in.rows())
        }
    }

    void SampleSomProblem::getRi(const Variable &x, SomUtils::MatD &rOut, int i) const
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

    void SampleSomProblem::getRi(const SomUtils::MatD &xEig, SomUtils::MatD &rOut, int i) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

        int startId = i * rotSz;
        // int endId = (i+1) * rotSz;

        SomUtils::MatD rOutVec(SomUtils::MatD::Zero(rotSz, 1));
        rOutVec = xEig.block(startId, 0, rotSz, 1);

        rOut = rOutVec.reshaped(sz_.p_, sz_.d_);
    }

    void SampleSomProblem::getRotations(const SomUtils::MatD &xEig, SomUtils::VecMatD &rOut) const
    {
        // rOut already needs to have fixed size by here
        int rotSz = getRotSz(); // as a vector

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
        }
    }

    void SampleSomProblem::getTi(const Variable &x, SomUtils::MatD &tOut, int i) const
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

    void SampleSomProblem::getTi(const SomUtils::MatD &xEig, SomUtils::MatD &tOut, int i) const
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

    void SampleSomProblem::getTranslations(const SomUtils::MatD &xEig, SomUtils::MatD &tOut) const
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

    void SampleSomProblem::makePfrct(const SomUtils::MatD &T, SomUtils::MatD &P, double &frct) const
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

    void SampleSomProblem::makeLrPrBr(const SomUtils::VecMatD &R, SomUtils::MatD &Lr, SomUtils::MatD &Pr, SomUtils::MatD &Br) const
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

    void SampleSomProblem::computeHrr(const SomUtils::VecMatD &xR, const SomUtils::VecMatD &uR, const SomUtils::MatD &P, SomUtils::VecMatD &hRR) const
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
        stiefelTangentProj(xR, DGf, hRR);
        // for (int i = 0; i < sz_.n_; ++i)
        // {
        //     ROFL_VAR2(i, hRR[i]);
        // }
    }

    void SampleSomProblem::computeHtt(const SomUtils::MatD &uT, const SomUtils::MatD &LR, SomUtils::MatD &hTT) const
    {
        hTT = uT * (LR.transpose() + LR);
    }

    void SampleSomProblem::computeHrt(const SomUtils::VecMatD &xR, const SomUtils::MatD uT, SomUtils::VecMatD &hrt) const
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

        stiefelTangentProj(xR, W, hrt);
        // ROFL_VAR1(uT);
        // for (int i = 0; i < sz_.n_; ++i)
        // {
        //     ROFL_VAR2(i, hrt[i]);
        // }
    }

    void SampleSomProblem::computeHtr(const SomUtils::VecMatD &uR, SomUtils::MatD &htr) const
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

    int SampleSomProblem::getRotSz() const
    {
        return sz_.d_ * sz_.p_;
    }

    int SampleSomProblem::getTranslSz() const
    {
        return sz_.p_;
    }

    void SampleSomProblem::stiefelTangentProj(const SomUtils::VecMatD &Y, const SomUtils::VecMatD &Hin, SomUtils::VecMatD &Hout) const
    {
        ROFL_ASSERT_VAR1(Y.size() == sz_.n_, Y.size());

        SomUtils::MatD tmp = Y[0].transpose() * Hin[0];
        Hout.clear();
        Hout.resize(sz_.n_, SomUtils::MatD::Zero(tmp.rows(), tmp.cols()));

        for (int i = 0; i < sz_.n_; ++i)
        {
            stiefelTangentProj(Y[i], Hin[i], Hout[i]);
        }
    }

    void SampleSomProblem::stiefelTangentProj(const SomUtils::MatD &Y, const SomUtils::MatD &Hin, SomUtils::MatD &Hout) const
    {
        SomUtils::MatD tmp = Y.transpose() * Hin;

        ROFL_ASSERT(tmp.rows() == tmp.cols()); // kind of useless as the check is performed also in extractSymmetricPart()

        SomUtils::MatD sympart(SomUtils::MatD::Zero(tmp.rows(), tmp.cols()));
        extractSymmetricPart(Y.transpose() * Hin, sympart);
        Hout = Hin - Y * sympart;
    }

    void SampleSomProblem::extractSymmetricPart(const SomUtils::MatD &in, SomUtils::MatD &out) const
    {
        ROFL_ASSERT(in.rows() == in.cols());
        out = 0.5 * (in + in.transpose());
    }

} // end of namespace ROPTLIB
