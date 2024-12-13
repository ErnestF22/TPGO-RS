#include "som_utils.h"

namespace SomUtils
{
    void deserializeRow(const std::string &row, double &pt)
    {
        std::stringstream ss(row);
        std::vector<std::string> ptCoordStrings;

        // ROFL_VAR1("\n");
        for (std::string strI; ss >> strI;)
        {
            ptCoordStrings.push_back(strI);

            if (ss.peek() == ',')
                ss.ignore();

            strI.erase(std::remove(strI.begin(), strI.end(), ','), strI.end()); // this should do nothing if line separator is null
            // ROFL_VAR1(strI);

            double ptCoord = std::stod(strI);

            // ROFL_VAR1(ptCoord);
            pt = ptCoord;
        }
        // ROFL_VAR1(pt.transpose());
    }

    void readCsvInitguess(std::string fname, ROPTLIB::Vector &csvVec)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }

        std::string line;
        // getline(fout, csvVec.header, '\n');
        // ROFL_VAR1(csvVec.header);
        csvVec.ObtainWriteEntireData();

        /////
        // Vector rgTiVec(sz_.p_, sz_.n_);
        realdp *GroptlibWriteArray = csvVec.ObtainWriteEntireData();

        /////

        int j = 0;
        while (getline(fout, line, '\n'))
        {
            // ROFL_VAR1(line);

            // add all the column data
            // of a row to a vector

            deserializeRow(line, GroptlibWriteArray[j]);
            j++;
            // ROFL_VAR2(j, line);

            // csvVec.pts.push_back(pt);
        }

        // rgTiVec.CopyTo(result->GetElement(gElemIdx));
        // csvVec.Print("csv read result");

        fout.close();
    }

    void deserializeRowTijs(const std::string &row, Eigen::RowVectorXd &tij)
    {
        std::stringstream ss(row);
        std::vector<std::string> ptCoordStrings;

        // ROFL_VAR1("\n");
        int idx = 0;
        for (std::string strI; ss >> strI;)
        {
            ptCoordStrings.push_back(strI);

            if (ss.peek() == ',')
                ss.ignore();

            strI.erase(std::remove(strI.begin(), strI.end(), ','), strI.end()); // this should nothing
            // ROFL_VAR2(idx,strI);
            tij(idx) = stod(strI);
            idx++; // idx should not go higher that d-1
        }
    }

    void readCsvTijs(std::string fname, Eigen::MatrixXd &Tijs, int d, int numEdges)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }

        std::string line;
        // getline(fout, csvVec.header, '\n');
        // ROFL_VAR1(csvVec.header);

        /////

        int j = 0;
        while (getline(fout, line, '\n'))
        {
            // ROFL_VAR1(line);

            // add all the column data
            // of a row to a vector

            Eigen::RowVectorXd tijI(numEdges);
            deserializeRowTijs(line, tijI);
            ROFL_VAR1(tijI)
            Tijs.row(j) = tijI;
            j++;
            // ROFL_VAR1(line);

            // csvVec.pts.push_back(pt);
        }

        // rgTiVec.CopyTo(result->GetElement(gElemIdx));
        // csvVec.Print("csv read result");

        fout.close();
    }

    void deserializeRowEdges(const std::string &row, Eigen::Vector2i &edgeI)
    {
        std::stringstream ss(row);
        std::vector<std::string> ptCoordStrings;

        // ROFL_VAR1("\n");
        int idx = 0;
        for (std::string strI; ss >> strI;)
        {
            ptCoordStrings.push_back(strI);

            if (ss.peek() == ',')
                ss.ignore();

            strI.erase(std::remove(strI.begin(), strI.end(), ','), strI.end()); // this should nothing
            // ROFL_VAR2(idx,strI);
            edgeI(idx) = stoi(strI);
            idx++; // idx should not go higher that 1
        }
    }

    void readCsvEdges(std::string fname, Eigen::MatrixXi &edges)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }

        std::string line;
        // getline(fout, csvVec.header, '\n');
        // ROFL_VAR1(csvVec.header);

        /////

        int j = 0;
        while (getline(fout, line, '\n'))
        {
            // ROFL_VAR1(line);

            // add all the column data
            // of a row to a vector

            Eigen::Vector2i edgeI;
            deserializeRowEdges(line, edgeI);
            edges.row(j) = edgeI;
            j++;
            // ROFL_VAR1(line);

            // csvVec.pts.push_back(pt);
        }

        // rgTiVec.CopyTo(result->GetElement(gElemIdx));
        // csvVec.Print("csv read result");

        fout.close();
    }

    double LinesearchInput(integer iter, const ROPTLIB::Variable &x1, const ROPTLIB::Vector &exeta1, realdp initialstepsize, realdp initialslope, const ROPTLIB::Problem *prob, const ROPTLIB::Solvers *solver)
    {
        return 1;
    };

    int miniDist(const std::vector<int> &distance, const std::vector<bool> &Tset) // finding minimum distance
    {
        int minimum = INT_MAX, ind;

        int n = distance.size();
        ROFL_ASSERT(n == Tset.size())

        for (int k = 0; k < n; k++)
        {
            if (Tset[k] == false && distance[k] <= minimum)
            {
                minimum = distance[k];
                ind = k;
            }
        }
        return ind;
    }

    void DijkstraAlgo(const Eigen::MatrixXi &adjMat, int src) // adjacency matrix
    {
        int n = adjMat.rows();

        ROFL_ASSERT(n == adjMat.cols())
        std::vector<int> distance(adjMat.rows()); // // array to calculate the minimum distance for each node
        std::vector<bool> Tset(adjMat.cols());    // boolean array to mark visited and unvisited for each node

        for (int k = 0; k < 6; k++)
        {
            distance[k] = INT_MAX;
            Tset[k] = false;
        }

        distance[src] = 0; // Source vertex distance is set 0

        for (int k = 0; k < n; k++)
        {
            int m = miniDist(distance, Tset);
            Tset[m] = true;
            for (int k = 0; k < n; k++)
            {
                // updating the distance of neighbouring vertex
                if (!Tset[k] && adjMat(m, k) && distance[m] != INT_MAX && distance[m] + adjMat(m, k) < distance[k])
                    distance[k] = distance[m] + adjMat(m, k);
            }
        }
        std::cout << "Vertex\t\tDistance from source vertex" << std::endl;
        for (int k = 0; k < n; k++)
        {
            char str = 65 + k;
            std::cout << str << "\t\t\t" << distance[k] << std::endl;
        }
    }

    void readMatlabCsvTijs(std::string fname, Eigen::MatrixXd &Tijs, int d, int numEdges)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }

        std::string line;
        // getline(fout, csvVec.header, '\n');
        // ROFL_VAR1(csvVec.header);

        /////

        int j = 0;
        while (getline(fout, line, '\n'))
        {
            std::stringstream ss(line);
            std::vector<std::string> result;

            while (ss.good())
            {
                std::string substr;
                getline(ss, substr, ',');
                result.push_back(substr);
                // ROFL_VAR1(substr)
            }

            Eigen::VectorXd rowI(numEdges, 1);

            for (int idx = 0; idx < result.size(); ++idx)
                rowI(idx, 0) = stod(result[idx]);
            Tijs.row(j) = rowI;
            j++;
        }

        fout.close();
    }

    void readMatlabCsvEdges(std::string fname, Eigen::MatrixXi &edges)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }

        std::string line;
        // getline(fout, csvVec.header, '\n');
        // ROFL_VAR1(csvVec.header);

        /////

        int j = 0;
        while (getline(fout, line, '\n'))
        {
            std::stringstream ss(line);
            std::vector<std::string> result;

            while (ss.good())
            {
                std::string substr;
                getline(ss, substr, ',');
                result.push_back(substr);
                // ROFL_VAR1(substr)
            }

            Eigen::Vector2i edgeI;

            for (int e = 0; e < result.size(); ++e)
                edgeI(e, 0) = stoi(result[e]);
            edges.row(j) = edgeI;
            j++;
        }

        fout.close();
    }

    void readSingleIntCsv(std::string fname, int &out)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return;
        }

        std::string line;
        getline(fout, line, '\n');

        out = stoi(line);
    }

} // end of namespace SomUtils