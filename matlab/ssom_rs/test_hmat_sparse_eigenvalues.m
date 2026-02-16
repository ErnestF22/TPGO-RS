function test_hmat_sparse_eigenvalues
hmat = readmatrix("../data/hmat.csv");

hmat = reshape(hmat, 80, 80);

[~, eigvals] = eig(hmat);


disp("sort(diag(real(eigvals)), 'descend')")
disp(sort(diag(real(eigvals)), 'descend'))

end %file function
