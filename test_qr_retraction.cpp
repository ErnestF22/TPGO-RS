#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include <filesystem>
#include <boost/lexical_cast.hpp>

#include <rofl/common/macros.h>

#include "thirdparty/qr_unique_sizeless/main.h"

#include "thirdparty/roptlib/problems_SampleSom.h"



// function Y = retraction_qr_stiefel(X, U, t)
//     % It is necessary to call qr_unique rather than simply qr to ensure
//     % this is a retraction, to avoid spurious column sign flips.
//     if nargin < 3
//         Y = qr_unique(X + U);
//     else
//         Y = qr_unique(X + t*U);
//     end
// end


int main(int argc, char **argv)
{

    SomUtils::MatD X1(4, 3);
    X1 << -0.247290964128394, -0.709429062803934, 0.121926006253623,
        -0.952190182267795, 0.242099566705116, 0.138797585524395,
        0.0716506430449076, 0.660194541380901, 0.00964725947025160,
        0.164460394030062, 0.0473431323550654, 0.982739136107011;

    SomUtils::MatD U1(4, 3);
    U1 << -0.450304594672209, 0.190944660134435, 0.197705433478300,
        0.175317099774566, 0.656906194094820, 0.0569067264599740,
        0.433825084892155, -0.0446722018118250, -0.164969390406632,
        0.148942349700657, 0.124988141333314, -0.0309466129108940;

    SomUtils::MatD Y1mlab(4, 3);

    Y1mlab << -0.580544344017758, -0.414199930254277, 0.265923014231035,
        -0.646519703818813, 0.750269857223769, 0.0604195917632961,
        0.420660755634195, 0.496932079523241, -0.131098259878416,
        0.260816153396013, 0.136352731763924, 0.953125211968653;

    SomUtils::MatD Y1(4, 3);
    SomUtils::MatD tmp(3, 3);
    main_qr_unique(X1 + U1, Y1, tmp);

    std::cout << "Y1mlab" << std::endl
              << Y1mlab << std::endl
              << "Y1" << std::endl
              << Y1 << std::endl;

    return 0;
}
