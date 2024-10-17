function [outputArg1] = GenerateRandomSymmetricMat(n)
%Generate and return random symmetric matrix of size nxn
d = 100*rand(n,1); % The diagonal values
t = triu(bsxfun(@min,d,d.').*rand(n),1); % The upper trianglar random values
M = diag(d)+t+t.'; % Put them together in a symmetric matrix

outputArg1 = M;

end