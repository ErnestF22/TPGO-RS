function transf = make_transf(rotMats, translMats)
% Take as input multiple rotation matrices and an equal number of
% translation matrices and composes all the corresponding transformation matrices
% This version works with arbitrary values of d;
%The built-in functions like rigidtform3d work only with fixed dimensions

%OBS. rot, transl were already MATLAB function names
%rotMats is d-D i.e. has size (d,d,N)
%translMats is given as already vectorized i.e. with size (d*N,1)
    d = size(rotMats, 2); %!! rotMats is passed as 3D!
    nrows_aff = size(rotMats, 3)*(d+1);
    N = size(rotMats, 3);
    transf = zeros(nrows_aff, d+1);
    for ii = 1:N
        transf_i = eye(d+1);
        transf_i(1:d, 1:d) = rotMats(:,:,ii);
        transl_i_aff = [translMats((ii-1)*d+1:(ii-1)*d+d, :); 1];
        transf_i(:, d+1) = transl_i_aff;
    
        transf((ii-1)*(d+1)+1:(ii-1)*(d+1)+(d+1), :) = transf_i; 
    end
end