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

namespace ROPTLIB
{

    class SampleSomProblem : public Problem
    {
    public:
        SampleSomProblem(std::vector<double> &rshSrc, std::vector<double> &rshDst, int arsLmax);

        virtual ~SampleSomProblem();
        virtual realdp f(const Variable &x) const;
        virtual Vector &EucGrad(const Variable &x, Vector *result) const;
        // virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

        void RoptToEig(Vector x, Eigen::Matrix3d &xEigen) const;

        void rsh_corr_grad(const Eigen::Vector3d &ea,
                           double &corr,
                           Eigen::Vector3d &grad,
                           Eigen::Matrix3d &Hessian,
                           std::vector<Eigen::MatrixXd> &Mf) const;

        void rsh_corr_grad(const Eigen::Vector3d &ea,
                           double &corr,
                           Eigen::Vector3d &grad,
                           Eigen::Matrix3d &Hessian) const;

        // private:
        std::vector<double> rshSrc_;
        std::vector<double> rshDst_;
        int lmax_;
    };

} // end of namespace