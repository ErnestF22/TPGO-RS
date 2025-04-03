function test_hmat_sparse_eigenvalues
hmat = readmatrix("../data/hmat.csv");

hmat = reshape(hmat, 80, 80);

eig(hmat)

end %file function
