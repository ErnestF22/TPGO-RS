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

    bool readCsvInitguess(std::string fname, ROPTLIB::Vector &csvVec)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return false;
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
        return true;
    }

    bool readCsvEigen(std::string fname, SomUtils::MatD &csvEig)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return false;
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
        return true;
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

    bool readCsvTijs(std::string fname, Eigen::MatrixXd &Tijs, int d, int numEdges)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return false;
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
        return true;
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

    bool readCsvEdges(std::string fname, Eigen::MatrixXi &edges)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return false;
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
        return true;
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

    bool readMatlabCsvTijs(std::string fname, Eigen::MatrixXd &Tijs, int d, int numEdges)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return false;
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
            {
                // ROFL_VAR1(result[idx])
                rowI(idx, 0) = stod(result[idx]);
            }
            Tijs.row(j) = rowI;
            j++;
        }

        fout.close();
        return true;
    }

    bool readMatlabCsvEdges(std::string fname, Eigen::MatrixXi &edges)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return false;
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
        return true;
    }

    bool readSingleIntCsv(std::string fname, int &out)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return false;
        }

        std::string line;
        getline(fout, line, '\n');

        out = stoi(line);

        return true;
    }

    bool readSingleDoubleCsv(std::string fname, double &out)
    {
        std::fstream fout;
        fout.open(fname, std::ios::in);

        if (!fout.is_open())
        {
            std::cerr << "Error opening file!" << std::endl;
            return false;
        }

        std::string line;
        getline(fout, line, '\n');

        out = stod(line);
        return true;
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

    std::string generateStampedString(const std::string prefix, const std::string postfix)
    {
        boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
        std::ostringstream formatter;
        std::string formatstring = prefix + "%Y%m%d_%H%M_%S" + postfix;
        formatter.imbue(std::locale(std::cout.getloc(), new boost::posix_time::time_facet(formatstring.c_str())));
        formatter << now;
        return formatter.str();
    }

    bool readCsvVecEigen(const std::string &filenameIn, Eigen::MatrixXd &out)
    {
        std::string line;
        std::ifstream fileIn(filenameIn);
        if (!fileIn.is_open())
        {
            ROFL_ERR("Error opening file")
            ROFL_VAR1(filenameIn)
            return false;
        }
        int ctr = 0;
        while (std::getline(fileIn, line))
        {
            // ROFL_VAR1(line);
            double val = std::stod(line);
            out(ctr, 0) = val;
            ctr++;
        }
        return true;
    }

    void vstack(const SomUtils::VecMatD &in, SomUtils::MatD &out)
    {
        ROFL_ASSERT(out.cols() == in[0].cols() && out.rows() == in[0].rows() * in.size());

        // out has same number of columns, whereas number of rows is the product of the size of the other 2 dimensions

        int rowJump = in[0].rows();
        for (int i = 0; i < in.size(); ++i)
        {
            out.block(rowJump * i, 0, rowJump, out.cols()) = in[i];
        }
    }

    void hstack(const SomUtils::VecMatD &in, SomUtils::MatD &out)
    {
        // out has same number of columns, whereas number of rows is the product of the size of the other 2 dimensions

        ROFL_ASSERT_VAR5(out.rows() == in[0].rows() && out.cols() == in[0].cols() * in.size(),
                         out.rows(), out.cols(), in[0].rows(), in[0].cols(), in.size());

        int colJump = in[0].cols();
        for (int i = 0; i < in.size(); ++i)
        {
            // ROFL_VAR2(out.block(0, colJump * i, out.rows(), colJump), in[i]);
            out.block(0, colJump * i, out.rows(), colJump) = in[i];
        }
    }

    void unStackH(const SomUtils::MatD &in, SomUtils::VecMatD &out, int colsOut)
    {
        int n = (int)in.cols() / colsOut;

        int fixedSz = in.rows(); // size that does not change in the 3D->2D transition (here, number of rows)

        ROFL_ASSERT(n * colsOut == in.cols());

        out.clear();
        out.resize(n, SomUtils::MatD::Zero(in.rows(), colsOut));

        for (int i = 0; i < n; ++i)
        {
            out[i] = in.block(0, colsOut * i, fixedSz, colsOut);
            ROFL_ASSERT(out[i].rows() == in.rows())
        }
    }

    void unStackV(const SomUtils::MatD &in, SomUtils::VecMatD &out, int rowsOut)
    {
        int n = (int)in.rows() / rowsOut;

        int fixedSz = rowsOut; // size that does not change in the 3D->2D transition (here, number of columns)

        // ROFL_VAR3(n, rowsOut, in.rows());
        ROFL_ASSERT(n * rowsOut == in.rows());

        out.clear();
        out.resize(n, SomUtils::MatD::Zero(rowsOut, in.cols()));

        for (int i = 0; i < n; ++i)
        {
            out[i] = in.block(rowsOut * i, 0, rowsOut, fixedSz);
            ROFL_ASSERT(out[i].cols() == in.cols())
        }
    }

    void stiefelTangentProj(const SomUtils::VecMatD &Y, const SomUtils::VecMatD &Hin, SomUtils::VecMatD &Hout)
    {
        int n = Hin.size();
        ROFL_ASSERT_VAR2(Y.size() == n, Y.size(), n);

        SomUtils::MatD tmp = Y[0].transpose() * Hin[0];
        Hout.clear();
        Hout.resize(n, SomUtils::MatD::Zero(tmp.rows(), tmp.cols()));

        for (int i = 0; i < n; ++i)
        {
            stiefelTangentProj(Y[i], Hin[i], Hout[i]);
        }
    }

    void stiefelTangentProj(const SomUtils::MatD &Y, const SomUtils::MatD &Hin, SomUtils::MatD &Hout)
    {
        SomUtils::MatD tmp = Y.transpose() * Hin;

        ROFL_ASSERT(tmp.rows() == tmp.cols()); // kind of useless as the check is performed also in extractSymmetricPart()

        SomUtils::MatD sympart(SomUtils::MatD::Zero(tmp.rows(), tmp.cols()));
        extractSymmetricPart(Y.transpose() * Hin, sympart);
        Hout = Hin - Y * sympart;
    }

    void extractSymmetricPart(const SomUtils::MatD &in, SomUtils::MatD &out)
    {
        ROFL_ASSERT(in.rows() == in.cols());
        out = 0.5 * (in + in.transpose());
    }

    void catZeroRow(const SomUtils::MatD &mIn, SomUtils::MatD &mOut)
    {
        ROFL_ASSERT(mOut.rows() == mIn.rows() + 1);
        ROFL_ASSERT(mOut.cols() == mIn.cols());

        mOut.setZero();
        mOut.block(0, 0, mIn.rows(), mIn.cols()) = mIn;
    }

    void catZeroRow3dArray(const SomUtils::VecMatD &mIn, SomUtils::VecMatD &mOut)
    {
        ROFL_ASSERT(mIn.size() == mOut.size())

        int n = mIn.size();
        for (int i = 0; i < n; ++i)
        {
            catZeroRow(mIn[i], mOut[i]);
        }
    }

    bool isEqualFloats(const SomUtils::MatD &a, const SomUtils::MatD &b, double thr)
    {
        ROFL_ASSERT(a.rows() == b.rows() && a.cols() == b.cols())

        double val = (a - b).cwiseAbs().maxCoeff();

        if (val > thr)
            return false;

        return true;
    }

    bool isEqualFloats(const SomUtils::VecMatD &a, const SomUtils::VecMatD &b, double thr)
    {
        int n = a.size();
        ROFL_ASSERT(n == b.size())
        bool retval = true;
        for (int i = 0; i < n; ++i)
        {
            if (!isEqualFloats(a, b, thr))
                return false;
        }
        return true;
    }

    bool isEqualDoubles(double a, double b, double thr)
    {
        return abs(a - b) < thr;
    }

    void multidet(const SomUtils::VecMatD &a3d, std::vector<double> &dets)
    {
        int n = a3d.size();
        dets.clear();
        dets.assign(n, 0.0);
        ROFL_ASSERT(n == dets.size())

        for (int i = 0; i < n; ++i)
            dets[i] = a3d[i].determinant();
    }

    void normalizeEucl(const SomUtils::MatD &mIn, SomUtils::MatD &mOut)
    {
        ROFL_ASSERT(mIn.rows() == mOut.rows() && mIn.cols() == mOut.cols())

        mOut = mIn;
        double normF = mIn.norm(); // TODO: maybe use .normalized() directly?

        mOut /= normF;
        // ROFL_VAR3(mIn, normF, mOut);
    }

    void normalizeEucl(const SomUtils::VecMatD &mIn, SomUtils::VecMatD &mOut)
    {
        ROFL_ASSERT(mIn.size() == mOut.size())
        std::for_each(mOut.begin(), mOut.end(), [](SomUtils::MatD &x) { //^^^ take argument by reference: LAMBDA FUNCTION
            x.setZero();
        });

        int n = mIn.size();
        for (int i = 0; i < n; ++i)
            normalizeEucl(mIn[i], mOut[i]);
    }

    double stlVecDoublesMean(const std::vector<double> &v)
    {
        int sz = v.size();
        double total = 0.0;
        for (int i = 0; i < sz; ++i)
        {
            total += v[i];
        }
        return total / sz;
    }

    void computeErrorsSingleRsom(const Eigen::MatrixXi &edges,
                                 const SomUtils::VecMatD &R, const SomUtils::MatD &T,
                                 const SomUtils::VecMatD &Rgt, const SomUtils::MatD &Tgt,
                                 std::vector<double> &rotErrs, std::vector<double> &translErrs)
    {
        // Compute errors
        ROFL_VAR1("Printing R, T out")
        for (auto &m : R)
            ROFL_VAR1(m)
        ROFL_VAR1(T)

        int n = R.size();
        int d = R[0].cols();
        int numEdges = edges.rows();

        for (int e = 0; e < numEdges; ++e)
        {
            int i = edges(e, 0) - 1;
            int j = edges(e, 1) - 1;

            Eigen::Matrix3d ri = R[i].block(0, 0, d, d);
            Eigen::Matrix3d rigt = Rgt[i].block(0, 0, d, d);
            Eigen::Matrix3d rj = R[j].block(0, 0, d, d);
            Eigen::Matrix3d rjgt = Rgt[j].block(0, 0, d, d);
            Eigen::Vector3d ti = T.col(i);
            Eigen::Vector3d tigt = Tgt.col(i);
            Eigen::Vector3d tj = T.col(j);
            Eigen::Vector3d tjgt = Tgt.col(j);

            Eigen::Matrix4d transfI(Eigen::Matrix4d::Identity());
            transfI.block(0, 0, d, d) = ri;
            transfI.block(0, d, d, 1) = ti;
            Eigen::Matrix4d transfJ(Eigen::Matrix4d::Identity());
            transfJ.block(0, 0, d, d) = rj;
            transfJ.block(0, d, d, 1) = tj;

            Eigen::Matrix4d transfIgt(Eigen::Matrix4d::Identity());
            transfIgt.block(0, 0, d, d) = rigt;
            transfIgt.block(0, d, d, 1) = tigt;
            Eigen::Matrix4d transfJgt(Eigen::Matrix4d::Identity());
            transfJgt.block(0, 0, d, d) = rjgt;
            transfJgt.block(0, d, d, 1) = tjgt;

            Eigen::MatrixXd p(Eigen::MatrixXd::Identity(d + 1, d + 1));
            SomUtils::computeRelativePose(transfI, transfJ, p);

            Eigen::MatrixXd pGt(Eigen::MatrixXd::Identity(d + 1, d + 1));
            SomUtils::computeRelativePose(transfIgt, transfJgt, pGt);

            Eigen::Matrix3d pR = p.block(0, 0, d, d);
            Eigen::Matrix3d pRgt = pGt.block(0, 0, d, d);

            double rotDistEdge = SomUtils::rotDistSingle(pR, pRgt);
            ROFL_VAR2(e, rotDistEdge);

            Eigen::Vector3d pT = p.block(0, d, d, 1);
            Eigen::Vector3d pTgt = pGt.block(0, d, d, 1);

            double translDistEdge = SomUtils::translErr(ti, tigt);
            // ROFL_VAR2(ti.transpose(), tigt.transpose());
            ROFL_VAR2(e, translDistEdge);

            rotErrs[e] = rotDistEdge;
            translErrs[e] = translDistEdge;
        }
    }

    double RElU(double x)
    {
        return (x > 0) ? x : 0;
    }

} // end of namespace SomUtils