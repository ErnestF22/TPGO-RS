function [outputArg1] = GenerateRandomSymmetricBooleanMat(n)
%Generate and return random symmetric boolean matrix of size nxn
d = rand(n,1); % The diagonal values
t = triu(bsxfun(@min,d,d.').*rand(n),1); % The upper trianglar random values
M = diag(d)+t+t.'; % Put them together in a symmetric matrix

M_bool = M;
M_bool(M >= 0.1) = 1; %0.1 was supposed to be 0.5 initially
M_bool(M_bool < 1) = 0;

outputArg1 = M_bool;

end