function stiefel_samples = generate_stiefel_sampling(nrs, d, N, num_samples)
%GENERATE_STIEFEL_SAMPLING Return an array of samples all in Stiefel 3D
%manifold, with size as given in the function parameters.

% arr_out

if ~exist('num_samples','var')
    num_samples=100;
end

% A simple method of generating such samples is as follows. 
% Draw nm random samples from N(0,1) and arrange them into an n×m matrix X. 
% Then X(X^T X)^{−1/2} is a random matrix that follows the uniform distribution 
% on the Stiefel manifold V_m(\real^{n}) 
% (e.g., Theorem 2.2.1 in Chikuse, Y. (2003). Statistics on Special Manifolds).

stiefel_samples = zeros(nrs, d, N, num_samples);



% N = 1
X = rand(nrs,d);

stief_rand_mat = X * inv(sqrtm(X' * X));



stiefel_samples(:,:,:,1) = stief_rand_mat;



end %file function

