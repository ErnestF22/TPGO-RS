function test_hmat_sparse_eigenvalues
hmat = readmatrix("../data/hmat.csv");

hmat = reshape(hmat, 80, 80);

[~, eigvals] = eig(hmat);


disp(sort(diag(real(eigvals))))

end %file function
