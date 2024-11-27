// From https://www.w3schools.com/dsa/dsa_ref_traveling_salesman.php

// def nearest_neighbor_tsp(distances):
//     n = len(distances)
//     visited = [False] * n
//     route = [0]
//     visited[0] = True
//     total_distance = 0

//     for _ in range(1, n):
//         last = route[-1]
//         nearest = None
//         min_dist = float('inf')
//         for i in range(n):
//             if not visited[i] and distances[last][i] < min_dist:
//                 min_dist = distances[last][i]
//                 nearest = i
//         route.append(nearest)
//         visited[nearest] = True
//         total_distance += min_dist

//     total_distance += distances[route[-1]][0]
//     route.append(0)
//     return route, total_distance

#include <iostream>
#include <numeric>
#include <bits/stdc++.h>

#include <rofl/common/io.h>
#include <rofl/common/macros.h>

#include <eigen3/Eigen/Dense>

int calculate_distance(std::vector<int> &route, const Eigen::MatrixXi &distances)
{
    int totalDistance = 0;
    for (int i = 0; i < route.size() - 1; ++i)
        totalDistance += distances(route[i], route[i + 1]);
    totalDistance += distances(route.back(), route.front());
    return totalDistance;
}

void brute_force_tsp(const Eigen::MatrixXi &distances, std::vector<int> &shortestRoute, int &minDistance)
{
    int n = distances.rows();
    ROFL_ASSERT(n == distances.cols())

    std::vector<int> cities(n);
    std::iota(std::begin(cities), std::end(cities), 0);
    shortestRoute.clear();
    minDistance = 1e+6;

    // for perm in permutations(cities)
    // {
    //     current_route = [0] + list(perm)
    //     current_distance = calculate_distance(current_route, distances)

    //     if (current_distance < min_distance)
    //     {
    //         min_distance = current_distance
    //         shortest_route = current_route
    //     }
    // }

    for (int i = 0; i < cities.size(); ++i)
    {
        ROFL_VAR2(i, cities[i])
    }

    // Generate all permuatation of an array
    int permNum = 0;
    do
    {
        std::vector<int> currentRoute(cities);
        for (int i = 0; i < currentRoute.size(); ++i)
        {
            ROFL_VAR3(permNum, i, currentRoute[i])
        }
        // currentRoute.insert(currentRoute.begin(), 0);
        int currentDistance = calculate_distance(currentRoute, distances);

        if (currentDistance < minDistance)
        {
            minDistance = currentDistance;
            shortestRoute = currentRoute;
        }
        permNum ++;
    } while (std::next_permutation(cities.begin(), cities.end()));

    shortestRoute.push_back(0);
    // return shortest_route, min_distance;
}

int main(int argc, char **argv)
{
    Eigen::MatrixXi distances(Eigen::MatrixXi(6, 6));
    distances << 0, 2, 2, 5, 9, 3,
        2, 0, 4, 6, 7, 8,
        2, 4, 0, 8, 6, 3,
        5, 6, 8, 0, 4, 9,
        9, 7, 6, 4, 0, 10,
        3, 8, 3, 9, 10, 0;

    // route, total_distance = nearest_neighbor_tsp(distances);
    std::vector<int> shortestRoute;
    int minDistance = 1e+6;
    brute_force_tsp(distances, shortestRoute, minDistance);

    // print("Route:", route);
    // print("Total distance:", total_distance);
    for (int i = 0; i < shortestRoute.size(); ++i)
    {
        ROFL_VAR2(i, shortestRoute[i])
    }
    ROFL_VAR1(minDistance);

    return 0;
}