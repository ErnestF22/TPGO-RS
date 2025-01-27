
#ifndef SOM_UTILS_H_
#define SOM_UTILS_H_

#include <string>
#include <sstream>
#include <fstream>
#include <vector>

#include <rofl/common/io.h>
#include <rofl/common/macros.h>

#include <eigen3/Eigen/Dense>

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

        SomSize() : p_(1), d_(1), n_(1) {}

        SomSize(int p, int d, int n) : p_(p), d_(d), n_(n)
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
     */
    void readCsvInitguess(std::string fname, ROPTLIB::Vector &csvVec);

    /**
     * Read Initguess from csv file in vectorized form (no size checks)
     * and save it into Eigen MatrixXd
     */
    void readCsvEigen(std::string fname, SomUtils::MatD &csvEig);

    /**
     * Called from readCsvTijs
     */
    void deserializeRowTijs(const std::string &row, Eigen::VectorXd &tij);

    /**
     * Read Tijs from csv file in vectorized form (no size checks)
     * Need numEdges
     * d is unused ATM but might be needed for more complex cases
     */
    void readCsvTijs(std::string fname, Eigen::MatrixXd &Tijs, int d, int numEdges);

    /**
     * Called from readCsvEdges
     */
    void deserializeRowEdges(const std::string &row, Eigen::Vector2i &edges);

    /**
     * Read Edges from rows in csv files (no size checks)
     */
    void readCsvEdges(std::string fname, Eigen::MatrixXi &edges);

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
     */
    void readMatlabCsvEdges(std::string fname, Eigen::MatrixXi &edges);

    /**
     * @brief No spaces after comma delimiter
     */
    void readMatlabCsvTijs(std::string fname, Eigen::MatrixXd &Tijs, int d, int numEdges);

    /**
     * @brief Read a single integer from a csv file and save it to reference param out
     */
    void readSingleIntCsv(std::string fname, int &out);

    /**
     * @brief Read a single double from a csv file and save it to reference param out
     */
    void readSingleDoubleCsv(std::string fname, double &out);

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

} // end of namespace SomUtils

#endif /*SOM_UTILS_H_*/