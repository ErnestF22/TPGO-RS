
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
     * Read Initguess as a row (no size checks)
     */
    void readCsvInitguess(std::string fname, ROPTLIB::Vector &csvVec);

} // end of namespace SomUtils

#endif /*SOM_UTILS_H_*/