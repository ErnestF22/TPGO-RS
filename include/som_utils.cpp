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

    void readCsvEigen(std::string fname, SomUtils::MatD &csvEig)
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
        // Vector rgTiVec(sz_.p_, sz_.n_);
        realdp Grealdp;

        /////

        int j = 0;
        while (getline(fout, line, '\n'))
        {
            // ROFL_VAR1(line);

            // add all the column data
            // of a row to a vector

            deserializeRow(line, csvEig.data()[j]);
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

    void readSingleDoubleCsv(std::string fname, double &out)
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

        out = stod(line);
    }

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

    double rotDistSingle(const Eigen::Matrix3d &R1, const Eigen::Matrix3d &R2)
    {
        //         R=R1'*R2;
        //         s1=R(6)-R(8);
        //         s2=R(7)-R(3);
        //         s3=R(2)-R(4);
        //         d=atan2(sqrt(s1*s1+s2*s2+s3*s3),R(1)+R(5)+R(9)-1);
        Eigen::Matrix3d R = R1.transpose() * R2;
        auto Rcolmaj = R.reshaped<Eigen::ColMajor>(9, 1);
        double s1 = Rcolmaj(5, 0) - Rcolmaj(7, 0);
        double s2 = Rcolmaj(6, 0) - Rcolmaj(2, 0);
        double s3 = Rcolmaj(1, 0) - Rcolmaj(3, 0);

        double d = atan2(sqrt(s1 * s1 + s2 * s2 + s3 * s3), Rcolmaj(0, 0) + Rcolmaj(4, 0) + Rcolmaj(8, 0) - 1);
        return d;
    }

    double rotDistSingle(const Eigen::Matrix2d &R1, const Eigen::Matrix2d &R2)
    {
        double theta1 = rot2ToAngle(R1);
        double theta2 = rot2ToAngle(R2);
        double d = abs(modAngle(theta1 - theta2));
        return d;
    }

    // function theta=rot2ToAngle(R)
    double rot2ToAngle(const Eigen::Matrix2d &R)
    {
        // theta=atan2(R(2,1)-R(1,2),trace(R));
        double theta = atan2(R(1, 0) - R(0, 1), R.trace());
        return theta;
    }

    double modAngle(double a)
    {
        // function a=modAngle(a)
        // a=mod(a+pi,2*pi)-pi;
        return (std::fmod(a + M_PI, 2 * M_PI)) - M_PI;
    }

    double translErr(const Eigen::Vector3d &T1, const Eigen::Vector3d &T2)
    {
        // [tij,lambdaij]=cnormalize(gij(1:3,4));
        // [tijtruth,lambdaijtruth]=cnormalize(gijtruth(1:3,4));
        // translErr=acos(max(-1,min(1,tij'*tijtruth)));

        Eigen::MatrixXd tij(T1.rows(), T1.cols()), tijtruth(T2.rows(), T2.cols());
        Eigen::MatrixXd lambdaij, lambdaijtruth; // when computing errors these are just scalars
        cnormalize(T1, tij, lambdaij);           // lambda!
        cnormalize(T2, tijtruth, lambdaijtruth); // lambda!

        auto tmp = tij.transpose() * tijtruth;

        return acos(std::max<double>(-1, std::min<double>(1, tmp(0, 0))));
    }

    void cnormalizeLambdas(const MatD &x, const Eigen::RowVectorXd &normx, MatD &xn)
    {
        // function [xn,normx] = cnormalize(x,normx)
        // d = size(x,1);
        // if ~exist('normx','var')
        //     normx = sqrt(sum(x.^2,1));
        // end
        // xn=zeros(size(x));

        ROFL_ASSERT(normx.cols() == x.cols())

        int r = x.rows();
        xn.resizeLike(x);
        xn.setZero();

        // idxNotZero=normx>1e-14;
        // for k=1:size(x,3)
        //     xn(:,idxNotZero(1,:,k),k) =  x(:,idxNotZero(1,:,k),k) ./ (repmat(normx(1,idxNotZero(1,:,k),k),d,1));
        // end
        for (int i = 0; i < normx.cols(); ++i)
        {
            if (normx(0, i) > 1e-14)
                xn.col(i) = x.col(i) / normx(0, i);
        }
    }

    void cnormalize(const MatD &x, MatD &xn, MatD &normx)
    {
        // function [xn,normx] = cnormalize(x,normx)
        // d = size(x,1);
        // if ~exist('normx','var')
        //     normx = sqrt(sum(x.^2,1));
        // end
        // xn=zeros(size(x));

        ROFL_ASSERT(xn.rows() == x.rows() && xn.cols() == x.cols())

        int r = x.rows();
        normx.resize(1, x.cols());
        xn.setZero();
        normx.setZero();

        // idxNotZero=normx>1e-14;
        // for k=1:size(x,3)
        //     xn(:,idxNotZero(1,:,k),k) =  x(:,idxNotZero(1,:,k),k) ./ (repmat(normx(1,idxNotZero(1,:,k),k),d,1));
        // end
        for (int i = 0; i < normx.cols(); ++i)
        {
            normx(0, i) = x.col(i).norm();
            if (normx(0, i) > 1e-14)
                xn.col(i) = x.col(i) / normx(0, i);
        }
    }

    void invg(const SomUtils::MatD &gIn, SomUtils::MatD &gOut)
    {
        // g1=g;
        // d=size(g,2)-1;
        // R=g(1:d,1:d,:);
        // for ig=1:size(g,3)
        //     g1(1:d,1:d,ig)=R(:,:,ig)';
        //     g1(1:d,d+1,ig)=-R(:,:,ig)'*g(1:d,d+1,ig);
        ROFL_ASSERT(gIn.rows() == 4 && gIn.cols() == 4 && gOut.rows() == 4 && gOut.cols() == 4)

        int d = 3;
        auto R = gIn.block(0, 0, d, d);

        gOut.block(0, 0, d, d) = R.transpose();
        gOut.block(0, d, d, 1) = -R.transpose() * gIn.block(0, d, d, 1);
    }

    void computeRelativePose(const Eigen::MatrixXd &g1, const Eigen::MatrixXd &g2, Eigen::MatrixXd &pose)
    {
        // pose.resize(4, 4);
        // pose.setIdentity();

        SomUtils::MatD g1inv(4, 4);
        invg(g1, g1inv);

        pose = g1inv * g2;
    }

} // end of namespace SomUtils