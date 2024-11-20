
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
     */
    void readCsvInitguess(std::string fname, ROPTLIB::Vector &csvVec);

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

} // end of namespace SomUtils

#endif /*SOM_UTILS_H_*/