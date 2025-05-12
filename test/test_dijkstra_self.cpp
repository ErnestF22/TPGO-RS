// From https://www.geeksforgeeks.org/dijkstras-shortest-path-algorithm-greedy-algo-7/
// C++ program for Dijkstra's single source shortest path
// algorithm. The program is for adjacency matrix
// representation of the graph
#include <iostream>

#include <limits.h>
#include <vector>
#include <eigen3/Eigen/Dense>

#include "thirdparty/roptlib/problems_SampleSom.h"

#include <rofl/common/macros.h>



int main(int argc, char **argv)
{

    /* Let us create the example graph discussed above */
    int n = 9;
    Eigen::MatrixXi adjmat(n, n);
    adjmat << 0, 4, 0, 0, 0, 0, 0, 8, 0,
        4, 0, 8, 0, 0, 0, 0, 11, 0,
        0, 8, 0, 7, 0, 4, 0, 0, 2,
        0, 0, 7, 0, 9, 14, 0, 0, 0,
        0, 0, 0, 9, 0, 10, 0, 0, 0,
        0, 0, 4, 14, 10, 0, 2, 0, 0,
        0, 0, 0, 0, 0, 2, 0, 1, 6,
        8, 11, 0, 0, 0, 0, 1, 0, 7,
        0, 0, 2, 0, 0, 0, 6, 7, 0;

    // Function call

    Eigen::MatrixXi edges(27, 2);
    edges << 1, 2,
        1, 8,
        2, 1,
        2, 3,
        2, 8,
        3, 2,
        3, 4,
        3, 9,
        4, 3,
        4, 5,
        4, 6,
        5, 3,
        5, 6,
        6, 3,
        6, 4,
        6, 5,
        6, 7,
        7, 6,
        7, 8,
        7, 9,
        8, 1,
        8, 2,
        8, 7,
        8, 9,
        9, 3,
        9, 7,
        9, 8;

    ROPTLIB::SampleSomProblem Prob;

    std::vector<double> dist;
    std::vector<int> prev;

    int src = 0;
    Prob.dijkstraSP(n, src, adjmat, dist, prev);


    ROFL_VAR1("Dijkstra self implementation output")

    for (int i = 0; i < dist.size(); ++i)
    {
        ROFL_VAR1(dist[i]);
    }

    // Vertex      Distance from Source
    // 0                 0
    // 1                 4
    // 2                 12
    // 3                 19
    // 4                 21
    // 5                 11
    // 6                 9
    // 7                 8
    // 8                 14

    std::vector<std::vector<int>> listNodes;

    Prob.dijkstraBT(src, n, prev, listNodes);

    std::cout << std::endl;

    ROFL_VAR1("list nodes output")
    for (int i = 0; i < listNodes.size(); ++i)
    {
        for (int j = 0; j < listNodes[i].size(); ++j)
            ROFL_VAR2(i, listNodes[i][j]);
        std::cout << std::endl;
    }

    // The minimum distance from 0 to 2 = 12. 0->1->2
    // The minimum distance from 0 to 3 = 19. 0->1->2->3
    // The minimum distance from 0 to 4 = 21. 0->7->6->5->4
    // The minimum distance from 0 to 5 = 11. 0->7->6->5
    // The minimum distance from 0 to 6 = 9. 0->7->6
    // The minimum distance from 0 to 7 = 8. 0->7
    // The minimum distance from 0 to 8 = 14. 0->1->2->8

    std::vector<std::vector<int>> listEdges;
    Prob.dijkstraBTedges(src, n, prev, edges, listEdges);

    std::cout << std::endl;

    ROFL_VAR1("list edges output")
    for (int i = 0; i < listEdges.size(); ++i)
    {
        for (int j = 0; j < listEdges[i].size(); ++j)
            ROFL_VAR2(i, listEdges[i][j]);
        std::cout << std::endl;
    }

    return 0;
}
