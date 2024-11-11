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

using MatD = Eigen::MatrixXd;
using VecMatD = std::vector<MatD>;

namespace ROPTLIB
{

    class SampleSomProblem : public Problem
    {
    public:
        /**
         * Default (empty) constructor
         */
        SampleSomProblem();

        /**
         * Standard constructor with problem data inputs (most commonly used)
         */
        SampleSomProblem(SomSize somSz, MatD &Tijs, Eigen::MatrixXi &edges);

        /**
         * Standard (empty) destructor
         */
        virtual ~SampleSomProblem();

        /**
         * Compute cost function, given input x
         * OBS. ROPT to Eig conversion performed here internally
         */
        virtual realdp f(const Variable &x) const;

        /**
         * Computation of Euclidean Gradient of cost function
         */
        // virtual Vector &EucGrad(const Variable &x, Vector *result) const;

        /**
         * Riemannian Gradient (backup: not used, as ROPTLIB wants Grad() directly)
         */
        virtual Vector &RieGrad(const Variable &x, Vector *result) const;

        /**
         * Computation of Riemannian Gradient of cost function
         */
        // virtual Vector &Grad(const Variable &x, Vector *result) const;

        /**
         * Function that computes Euclidean gradient of Translation estimation cost (with Eigen inputs/outputs)
         */
        void egradR(const MatD &P, VecMatD &egR) const;

        /**
         * Function that computes Riemannian gradient of Rotation estimation cost (with Eigen inputs/outputs)
         */
        void rgradR(const VecMatD &R, const MatD &P, VecMatD &rgR) const;

        /**
         * Function that computes Euclidean gradient of Translation estimation cost (with Eigen inputs/outputs)
         */
        void egradT(const MatD &T, const MatD &Lr, const MatD &Pr, MatD &egT) const;

        /**
         * Function that computes Riemannian gradient of Translation estimation cost (with Eigen inputs/outputs)
         */
        void rgradT(const MatD &T, const MatD &Lr, const MatD &Pr, MatD &egT) const;

        /**
         * Euclidean Hessian action i.e., *result = H(x)[etax]
         */
        // virtual Vector &EucHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

        /**
         * Hessian action i.e., *result = H(x)[etax]
         */
        virtual Vector &RieHessianEta(const Variable &x, const Vector &etax, Vector *result) const;

        /**
         * Convert ROPTLIB Vector into Eigen equivalent
         * TOCHECK: for now seems to work only with vectors (not matrices, tensors)
         */
        void RoptToEig(Vector x, MatD &xEigen) const;

        /**
         * Vertically stack input vector
         */
        void vstack(const VecMatD &in, MatD &out) const;

        /**
         * Horizontally stack input vector
         */
        void hstack(const VecMatD &in, MatD &out) const;

        /**
         * Unstack a vertically stacked "3D" array
         */
        void unStackV(const MatD &in, VecMatD &out, int rowsOut = 3) const;

        /**
         * Unstack a horizontally stacked "3D" array
         */
        void unStackH(const MatD &in, VecMatD &out, int colsOut = 3) const;

        /**
         * Get i-th Rotation from ROPTLIB variable x
         */
        void getRi(const Variable &x, MatD &rOut, int i) const;

        /**
         * Get i-th Rotation from Eigen-converted variable xEig
         */
        void getRi(const MatD &xEig, MatD &rOut, int i) const;

        /**
         * Get all roatations in vector of pxd matrices tOut (with n elements)
         */
        void getRotations(const MatD &xEig, VecMatD &rOut) const;

        /**
         * Get i-th Translation from ROPTLIB variable x
         */
        void getTi(const Variable &x, MatD &rOut, int i) const;

        /**
         * Get i-th Translation from Eigen-converted variable xEig
         */
        void getTi(const MatD &xEig, MatD &tOut, int i) const;

        /**
         * Get all translations in p * n matrix tOut
         */
        void getTranslations(const MatD &xEig, MatD &tOut) const;

        /**
         * Computes matrices used in rotation estimation cost
         */
        void makePfrct(const MatD &T, MatD &P, double &frct) const;

        /**
         * Computes matrices used in translation estimation cost
         */
        void makeLrPrBr(const VecMatD &R, MatD &Lr, MatD &Pr, MatD &Br) const;

        void computeHrr(const VecMatD &xR, const VecMatD &uR, const MatD &P, VecMatD &hRR) const;

        void computeHtt(const MatD &uT, const MatD &LR, MatD &hTT) const;

        void computeHrt(const VecMatD &xR, const MatD uT, VecMatD &hrt) const;

        void computeHtr(const VecMatD &uR, MatD &htr) const;

        /**
         * Return p * d (size of a sigle Stiefel-rotation)
         */
        int getRotSz() const;

        /**
         * Return p (size of a sigle Stiefel-translation)
         */
        int getTranslSz() const;

        /**
         * Project Hin onto tangent space at Y (3D Stiefel)
         */
        void stiefelTangentProj(const VecMatD &Y, const VecMatD &Hin, VecMatD &Hout) const;

        /**
         * Project Hin onto tangent space at Y (single Stiefel matrix)
         */
        void stiefelTangentProj(const MatD &Y, const MatD &Hin, MatD &Hout) const;

        /**
         * Extract symmetric part of input (square) matrix
         */
        void extractSymmetricPart(const MatD &in, MatD &out) const;

        // private: //TODO: separate public from private members

        /**
         * Struct that contains problem size info
         */
        SomSize sz_;

        /**
         * (Estimated) Relative translations between nodes
         * Size: d x e
         */
        MatD Tijs_;

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
    };

} // end of namespace ROPTLIB