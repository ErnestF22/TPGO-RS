function ti_3d = make_t_globalframe_3d(ti_codemeta)
%MAKE_T_GLOBALFRAME_3D function that take the T matrix in as input and produces a
%3D array with every Tij present in tijs_mat ordered "blk-column"-wise
%Note: here T is considered as given by codemeta, i.e. of size dxN

N = size(ti_codemeta, 2);
d = size(ti_codemeta, 1);

ti_3d = zeros(d, 1, N);

%TODO: version of this function without a for loop
for ii = 1:N
    tij = ti_codemeta(:, ii);
    ti_3d(:, :, ii) = tij;
end

end %function