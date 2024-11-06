#include "manifolds_Rotations.h"
#include "manifolds_Stiefel.h"

#include <vector>

#include <eigen3/Eigen/Dense>

#include <rofl/common/io.h>
#include <rofl/common/macros.h>

// #include <ars3d/RshRotation.h>

// #include <ars3d_test/rsh_utils.h>

#include "problems_Problem.h"
#include "others_def.h"

struct SomSize
{
    int p_;
    int d_;
    int n_;

    SomSize() : p_(1), d_(1), n_(1) {}

    SomSize(int p, int d, int n) : p_(p), d_(d), n_(n)
    {
    }
};

namespace ROPTLIB
{

    class SampleSomProblem : public Problem
    {
    public:
        SampleSomProblem();

        SampleSomProblem(SomSize somSz, Eigen::MatrixXd &Tijs, Eigen::MatrixXi &edges);

        virtual ~SampleSomProblem();

        virtual realdp f(const Variable &x) const;

        virtual Vector &EucGrad(const Variable &x, Vector *result) const;

        virtual Vector &RieGrad(const Variable &x, Vector *result) const;

        virtual Vector &Grad(const Variable &x, Vector *result) const;

        void egradR(const Eigen::MatrixXd &P, std::vector<Eigen::MatrixXd> &egR) const;

        void rgradR(const std::vector<Eigen::MatrixXd> &R, const Eigen::MatrixXd &P, std::vector<Eigen::MatrixXd> &rgR) const;

        void egradT(const Eigen::MatrixXd &T, const Eigen::MatrixXd &Lr, const Eigen::MatrixXd &Pr, Eigen::MatrixXd &egT) const;

        void rgradT(const Eigen::MatrixXd &T, const Eigen::MatrixXd &Lr, const Eigen::MatrixXd &Pr, Eigen::MatrixXd &egT) const;

        // virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

        // void RoptToEig(Vector x, Eigen::MatrixXf &xEigen) const;

        void RoptToEig(Vector x, Eigen::MatrixXd &xEigen) const;

        void vstack(const std::vector<Eigen::MatrixXd> &in, Eigen::MatrixXd &out) const;

        void hstack(const std::vector<Eigen::MatrixXd> &in, Eigen::MatrixXd &out) const;

        /**
         * Unstack a vertically stacked "3D" array
         */
        void unStackV(const Eigen::MatrixXd &in, std::vector<Eigen::MatrixXd> &out, int rowsOut = 3) const;

        void unStackH(const Eigen::MatrixXd &in, std::vector<Eigen::MatrixXd> &out, int colsOut = 3) const;

        void getRi(const Variable &x, Eigen::MatrixXd &rOut, int i) const;

        void getRi(const Eigen::MatrixXd &xEig, Eigen::MatrixXd &rOut, int i) const;

        void getRotations(const Eigen::MatrixXd &xEig, std::vector<Eigen::MatrixXd> &rOut) const;

        void getTi(const Variable &x, Eigen::MatrixXd &rOut, int i) const;

        void getTi(const Eigen::MatrixXd &xEig, Eigen::MatrixXd &tOut, int i) const;

        void getTranslations(const Eigen::MatrixXd &xEig, Eigen::MatrixXd &tOut) const;

        void makePfrct(const Eigen::MatrixXd &T, Eigen::MatrixXd &P, double &frct) const;

        void makeLrPrBr(const std::vector<Eigen::MatrixXd> &R, Eigen::MatrixXd &Lr, Eigen::MatrixXd &Pr, Eigen::MatrixXd &Br) const;

        int getRotSz() const;

        int getTranslSz() const;

        void stiefelTangentProj(const std::vector<Eigen::MatrixXd> &Y, const std::vector<Eigen::MatrixXd> &Hin, std::vector<Eigen::MatrixXd> &Hout) const
        {
            ROFL_ASSERT(Y.size() == sz_.n_);

            Eigen::MatrixXd tmp = Y[0].transpose() * Hin[0];
            Hout.clear();
            Hout.resize(sz_.n_, Eigen::MatrixXd::Zero(tmp.rows(), tmp.cols()));

            for (int i = 0; i < sz_.n_; ++i)
            {
                stiefelTangentProj(Y[i], Hin[i], Hout[i]);
            }
        }

        void stiefelTangentProj(const Eigen::MatrixXd &Y, const Eigen::MatrixXd &Hin, Eigen::MatrixXd &Hout) const
        {
            Eigen::MatrixXd tmp = Y.transpose() * Hin;

            ROFL_ASSERT(tmp.rows() == tmp.cols()); // kind of useless as the check is performed also in extractSymmetricPart()

            Eigen::MatrixXd sympart(Eigen::MatrixXd::Zero(tmp.rows(), tmp.cols()));
            extractSymmetricPart(Y.transpose() * Hin, sympart);
            Hout = Hin - Y * sympart;
        }

        void extractSymmetricPart(const Eigen::MatrixXd &in, Eigen::MatrixXd &out) const
        {
            ROFL_ASSERT(in.rows() == in.cols());
            out = 0.5 * (in + in.transpose());
        }

        // private:
        SomSize sz_;
        Eigen::MatrixXd Tijs_;
        Eigen::MatrixXi edges_;
        int numEdges_; // num edges
        int fullSz_;
    };

} // end of namespace