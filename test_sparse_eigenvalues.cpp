#include <iostream>
#include <fstream>

#include <rofl/common/param_map.h>
#include <rofl/common/profiler.h>

#include "som_utils.h"

#include <Spectra/SymEigsSolver.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SparseCore>


int main (int argc, char **argv)
{
    int vecsz = 80;
    SomUtils::MatD Hmat(vecsz * vecsz, 1);
    SomUtils::readCsvEigen("../data/hmat.csv", Hmat);
    

    Hmat = Hmat.reshaped<Eigen::ColMajor>(vecsz, vecsz);

    Eigen::SparseMatrix<double, Eigen::ColMajor> HmatSparse(Hmat.rows(), Hmat.cols());
    HmatSparse = Hmat.sparseView();
    // Eigen::EigenSolver<Eigen::SparseMatrix<double>> es;
    // es.compute(HmatSparse);
    // auto eigvals = es.eigenvalues();
    // auto eigvecs = es.eigenvectors();

    // Construct matrix operation object using the wrapper class DenseSymMatProd
    Spectra::SparseGenMatProd<double> op(HmatSparse);

    // Construct eigen solver object, requesting the largest three eigenvalues
    Spectra::GenEigsSolver<Spectra::SparseGenMatProd<double>> eigs(op, 3, 6);

    // Initialize and compute
    eigs.init();
    int nconv = eigs.compute(Spectra::SortRule::SmallestReal);

    // Retrieve results
    Eigen::VectorXcd evalues;
    if(eigs.info() == Spectra::CompInfo::Successful)
        evalues = eigs.eigenvalues();

    std::cout << "Eigenvalues found:\n" << evalues << std::endl;


    return 0;
}