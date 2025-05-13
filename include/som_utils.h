
#ifndef SOM_UTILS_H_
#define SOM_UTILS_H_

#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include <rofl/common/io.h>
#include <rofl/common/macros.h>

#include <eigen3/Eigen/Dense>

#include <boost/date_time/posix_time/posix_time.hpp>

#include "thirdparty/roptlib/manifolds_Element.h"
#include "thirdparty/roptlib/problems_Problem.h"
#include "thirdparty/roptlib/solvers_Solvers.h"

namespace SomUtils
{
    struct SomSize
    {
        int p_;
        int d_;
        int n_;
        int e_;

        SomSize() : p_(1), d_(1), n_(1), e_(0) {}

        SomSize(int p, int d, int n) : p_(p), d_(d), n_(n), e_(0)
        {
        }

        SomSize(int p, int d, int n, int e) : p_(p), d_(d), n_(n), e_(e)
        {
        }
    };

    using MatD = Eigen::MatrixXd;
    using VecMatD = std::vector<MatD>;
    using VecD = Eigen::VectorXd;

    /**
     * Called from readCsvInitguess
     */
    void deserializeRow(const std::string &row, double &pt);

    /**
     * Read Initguess from csv file in vectorized form (no size checks)
     * and save it into ROPTLIB vector
     * TODO: rename it as this method is not necessarily to be used only for reading initguesses
     * Return success boolean
     */
    bool readCsvInitguess(std::string fname, ROPTLIB::Vector &csvVec);

    /**
     * Read Initguess from csv file in vectorized form (no size checks)
     * and save it into Eigen MatrixXd
     * Return success boolean
     */
    bool readCsvEigen(std::string fname, SomUtils::MatD &csvEig);

    /**
     * Called from readCsvTijs
     */
    void deserializeRowTijs(const std::string &row, Eigen::VectorXd &tij);

    /**
     * Read Tijs from csv file in vectorized form (no size checks)
     * Need numEdges
     * d is unused ATM but might be needed for more complex cases
     */
    bool readCsvTijs(std::string fname, Eigen::MatrixXd &Tijs, int d, int numEdges);

    /**
     * Called from readCsvEdges
     */
    void deserializeRowEdges(const std::string &row, Eigen::Vector2i &edges);

    /**
     * Read Edges from rows in csv files (no size checks)
     */
    bool readCsvEdges(std::string fname, Eigen::MatrixXi &edges);

    // template <typename M>
    // M load_csv(const std::string &path)
    // {
    //     std::ifstream indata;
    //     indata.open(path);
    //     std::string line;
    //     std::vector<double> values;
    //     uint rows = 0;
    //     while (std::getline(indata, line))
    //     {
    //         std::stringstream lineStream(line);
    //         std::string cell;
    //         while (std::getline(lineStream, cell, ','))
    //         {
    //             values.push_back(std::stod(cell));
    //         }
    //         ++rows;
    //     }
    //     return Eigen::Map<const Eigen::Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, Eigen::ColMajor>>(values.data(), rows, values.size() / rows);
    // }

    /* User - specified linesearch algorithm */
    double LinesearchInput(integer iter, const ROPTLIB::Variable &x1, const ROPTLIB::Vector &exeta1, realdp initialstepsize, realdp initialslope, const ROPTLIB::Problem *prob, const ROPTLIB::Solvers *solver);

    /**
     * Used to compute minimum distance inside DijkstraAlgo
     */
    int miniDist(const std::vector<int> &distance, const std::vector<bool> &Tset);

    /**
     * Dijkstra Algorithm (shortest path) as per https://favtutor.com/blogs/dijkstras-algorithm-cpp
     */
    void DijkstraAlgo(const Eigen::MatrixXi &adjMat, int src); // adjacency matrix

    /**
     * @brief No spaces after comma delimiter
     * Return success boolean
     */
    bool readMatlabCsvEdges(std::string fname, Eigen::MatrixXi &edges);

    /**
     * @brief No spaces after comma delimiter
     * Return success boolean
     */
    bool readMatlabCsvTijs(std::string fname, Eigen::MatrixXd &Tijs, int d, int numEdges);

    /**
     * @brief Read a single integer from a csv file and save it to reference param out
     * Return success boolean
     */
    bool readSingleIntCsv(std::string fname, int &out);

    /**
     * @brief Read a single double from a csv file and save it to reference param out
     * Return success boolean
     */
    bool readSingleDoubleCsv(std::string fname, double &out);

    // function d=rot_distSingle(R1,R2)
    // switch size(R1,1)
    //     case 2
    //         theta1=rot2ToAngle(R1);
    //         theta2=rot2ToAngle(R2);
    //         d=abs(modAngle(theta1-theta2));
    //     case 3
    //
    //     otherwise
    //         error('Not implemented yet')
    // end
    // %[d(r1,r2),w] = angleaxis(dc2quat(R1(:,:,r1)'*R2(:,:,r2)));
    // %d(r1,r2)=acos(max(min((trace(R1(:,:,r1)'*R2(:,:,r2))-1)/2,1),-1));
    // %d(r1,r2)=-trace(hat(logrot(R1(:,:,r1)'*R2(:,:,r2)))^2)/2;

    /**
     * Compute the geodesic distance between the rotations R1 and R2 in SO(3)
     */
    double rotDistSingle(const Eigen::Matrix3d &R1, const Eigen::Matrix3d &R2);

    /**
     * Compute the geodesic distance between the rotations R1 and R2 in SO(2)
     */
    double rotDistSingle(const Eigen::Matrix2d &R1, const Eigen::Matrix2d &R2);

    // function theta=rot2ToAngle(R)
    /**
     * Return angle corresponding to input 2D rotation matrix
     */
    double rot2ToAngle(const Eigen::Matrix2d &R);

    // function a=modAngle(a)
    // a=mod(a+pi,2*pi)-pi;
    /**
     * Maps any scalar angle (in radians) to the equivalent between -pi and pi
     */
    double modAngle(double a);

    /**
     * Compute the error between translations T1 and T2 in R^3
     * // [tij,lambdaij]=cnormalize(gij(1:3,4));
     * // [tijtruth,lambdaijtruth]=cnormalize(gijtruth(1:3,4));
     * // translErr=acos(max(-1,min(1,tij'*tijtruth)));
     * !! this has part of the code for computing scale error!
     */
    double translErr(const Eigen::Vector3d &T1, const Eigen::Vector3d &T2);

    /**
     *  %function [xn,normx] = cnormalize(x,normx)
     *  %Normalize each column of x. The argument normx represents
     *  %a vector with the norms of the columns of x.
     *  %The function does not normalize vectors with norms less than 1e-15.
     *  !! Here normx is not computed, but instead just passed (useful when input x may be scaled arbitrarily)
     * */
    void cnormalizeLambdas(const MatD &x, const Eigen::RowVectorXd &normx, MatD &xn);

    /**
     *  %function [xn,normx] = cnormalize(x,normx)
     *  %Normalize each column of x. The argument normx represents
     *  %a vector with the norms of the columns of x.
     *  In this version, the norm of each column is also computed and returned as reference.
     *  %The function does not normalize vectors with norms less than 1e-15.
     * */
    void cnormalize(const MatD &x, MatD &xn, MatD &normx);

    /**
     * Compute the inverse of a rigid body transformation represented as a 4x4 matrix
     */
    void invg(const SomUtils::MatD &gIn, SomUtils::MatD &gOut);

    /**
     * Compute relative pose from two absolute poses
     */
    void computeRelativePose(const Eigen::MatrixXd &g1, const Eigen::MatrixXd &g2, Eigen::MatrixXd &pose);

    /**
     * Generate and return a string containing current timestamp
     */
    std::string generateStampedString(const std::string prefix, const std::string postfix);

    /**
     * Read a csv file with a single column and save it to reference @param out
     */
    bool readCsvVecEigen(const std::string &filenameIn, Eigen::MatrixXd &out);

    /**
     * Vertically stack input vector
     */
    void vstack(const SomUtils::VecMatD &in, SomUtils::MatD &out);

    /**
     * Horizontally stack input vector
     */
    void hstack(const SomUtils::VecMatD &in, SomUtils::MatD &out);

    /**
     * Unstack a vertically stacked "3D" array
     */
    void unStackV(const SomUtils::MatD &in, SomUtils::VecMatD &out, int rowsOut = 3);

    /**
     * Unstack a horizontally stacked "3D" array
     */
    void unStackH(const SomUtils::MatD &in, SomUtils::VecMatD &out, int colsOut = 3);

    /**
     * Project Hin onto tangent space at Y (3D Stiefel)
     */
    void stiefelTangentProj(const SomUtils::VecMatD &Y, const SomUtils::VecMatD &Hin, SomUtils::VecMatD &Hout);

    /**
     * Project Hin onto tangent space at Y (single Stiefel matrix)
     */
    void stiefelTangentProj(const SomUtils::MatD &Y, const SomUtils::MatD &Hin, SomUtils::MatD &Hout);

    /**
     * Extract symmetric part of input (square) matrix
     */
    void extractSymmetricPart(const SomUtils::MatD &in, SomUtils::MatD &out);

    /**
     * @brief Return @param mOut as a copy of @param mIn but with an extra zero-row at the bottom of the matrix
     */
    void catZeroRow(const SomUtils::MatD &mIn, SomUtils::MatD &mOut);

    /**
     * @brief 3D version of catZeroRow where catZeroRow is applied to each pair of elements
     * in @param muIn, @param muOut
     */
    void catZeroRow3dArray(const SomUtils::VecMatD &mIn, SomUtils::VecMatD &mOut);

    /**
     * @brief Return whether floating point elements matrices @param a, @param b are equal
     * up to a difference of @param thr
     * @return true is they are equal, @return false otherwise
     */
    bool isEqualFloats(const SomUtils::MatD &a, const SomUtils::MatD &b, double thr = 1e-5);

    /**
     * @brief 3D version of isEqualFloats() comparing all element pairs in @param a and @param b
     * @return true if they are all equal, @return false if at least one pair is not (floaty-)equal
     */
    bool isEqualFloats(const SomUtils::VecMatD &a, const SomUtils::VecMatD &b, double thr = 1e-5);

    /**
     * @brief Return whether double-precision floating point inputs @param a, @param b are equal
     * up to a difference of @param thr
     * @return true is they are equal, @return false otherwise
     */
    bool isEqualDoubles(double a, double b, double thr = 1e-6);

    /**
     * @brief Compute determinants of each element in input @param a3d
     * and return them in output reference vector @param dets
     */
    void multidet(const SomUtils::VecMatD &a3d, std::vector<double> &dets);

    /**
     * @brief Normalize @param mIn matrix in 2D Euclidean space
     * returning the normalized output in reference @param mOut
     */
    void normalizeEucl(const SomUtils::MatD &mIn, SomUtils::MatD &mOut);

    /**
     * @brief Apply 2D version of normalizeEucl to each pair of elements in @param mIn, @param mOut
     */
    void normalizeEucl(const SomUtils::VecMatD &mIn, SomUtils::VecMatD &mOut);

    /**
     * Return mean of input vector of doubles @param v
     */
    double stlVecDoublesMean(const std::vector<double> &v);

    /**
     * Compute errors of a single run of RSOM methods (RS, ICP, Procrustes)
     * and return them in the rotErrs and translErrs vectors (reference params).
     */
    void computeErrorsSingleRsom(const Eigen::MatrixXi &edges,
                                 const SomUtils::VecMatD &R, const SomUtils::MatD &T,
                                 const SomUtils::VecMatD &Rgt, const SomUtils::MatD &Tgt,
                                 std::vector<double> &rotErrs, std::vector<double> &translErrs);

} // end of namespace SomUtils

#endif /*SOM_UTILS_H_*/