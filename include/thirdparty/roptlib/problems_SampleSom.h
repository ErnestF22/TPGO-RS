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
        // virtual Vector &EucGrad(const Variable &x, Vector *result) const;
        // virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

        // void RoptToEig(Vector x, Eigen::MatrixXf &xEigen) const;

        void RoptToEig(Vector x, Eigen::MatrixXd &xEigen) const;

        void getRi(const Variable &x, Eigen::MatrixXd &rOut, int i) const;

        void getRi(const Eigen::MatrixXd &xEig, Eigen::MatrixXd &rOut, int i) const;

        void getTi(const Variable &x, Eigen::MatrixXd &rOut, int i) const;

        void getTi(const Eigen::MatrixXd &xEig, Eigen::MatrixXd &tOut, int i) const;

        int getRotSz() const;

        int getTranslSz() const;

        // private:
        SomSize sz_;
        Eigen::MatrixXd Tijs_;
        Eigen::MatrixXi edges_;
        int numEdges_; // num edges
        int fullSz_;
    };

} // end of namespace