function tijs_3d = make_tijs_3d(tijs_mat)
%MAKE_TIJS_3D function that take the Tijs matrix in as input and produces a
%3D array with every Tij present in tijs_mat ordered "blk-column"-wise

N = size(tijs_mat, 2);
d = size(tijs_mat, 1)/N;

tijs_3d = zeros(d, 1, N*N);

%TODO: version of this function without a for loop
ctr = 1;
for jj = 1:N
    for ii = 1:N
        tij = tijs_mat(ii*d-(d-1):ii*d, jj);
        tijs_3d(:, 1, ctr) = tij;
        ctr = ctr + 1;
    end
end

end %function