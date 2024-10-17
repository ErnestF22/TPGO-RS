function A_stiefel = make_rand_stiefel_mat(num_rows_stiefel, d, N)
%MAKE_RAND_STIEFEL_MAT Return A_stiefel, a 3D array with (:,:,ii) elements
%on Stiefel manifold

eps = 1e-5;

A = randn(num_rows_stiefel, d, N);  % random iid ~N(0,1)
A_stiefel = zeros(num_rows_stiefel, d, N); 
for ii = 1:N
    A_i = A(:,:,ii);
    oA = orth(A_i);
%     oA = orth( A_i.' ).'; % orthogonal rows
%     nA = bsxfun( @rdivide, oA, sqrt( sum( oA.^2, 2 ) ) ); % normalize to unit length
    A_stiefel(:,:,ii) = oA;
    if (max(oA'*oA - eye(d), [], 'all') >= eps)
        fprintf("oA is NOT on Stiefel!\n");
        break;
    end
end


end
