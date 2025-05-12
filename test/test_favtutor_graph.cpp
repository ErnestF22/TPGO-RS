#include "som_utils.h"

int main(int argc, char **argv)
{
    Eigen::MatrixXi graph(Eigen::MatrixXi(6, 6));
    graph << 0, 1, 2, 0, 0, 0,
        1, 0, 0, 5, 1, 0,
        2, 0, 0, 2, 3, 0,
        0, 5, 2, 0, 2, 2,
        0, 1, 3, 2, 0, 1,
        0, 0, 0, 2, 1, 0;
    SomUtils::DijkstraAlgo(graph, 0);
    return 0;
}