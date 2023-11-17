function ret_val = make_rand_stiefel_3d_array(n, p, k, array_type)
%MAKE_RAND_STIEFEL_MAT Make random Stiefel array of size n x p x k

if ~exist('array_type','var')
    array_type='double';
end

ret_val =  qr_unique(randn(n, p, k, array_type));

end

