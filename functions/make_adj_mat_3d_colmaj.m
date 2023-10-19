function adj_mat_3d = make_adj_mat_3d_colmaj(adj_mat)
%MAKE_ADJ_MAT_3D_COLMAJ takes an adjacency matrix as input and transforms
%it into a 3d array, where to each index in the third dimension corresponds 
%the boolean stating whether or not an edge exists between the
%corresponding row and column index, taken in column-major ordering

N = size(adj_mat, 1); % = size(adj_mat, 2)

adj_mat_3d = zeros(1,1,N*N);

%TODO: rewrite this function without for loops
ctr = 1;
for jj = 1:N
    for ii = 1:N
        adj_mat_3d(:,:,ctr) = adj_mat(ii, jj);
        ctr = ctr + 1;
    end
end

end %function