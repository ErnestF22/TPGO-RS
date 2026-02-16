function test_qr_retraction
% Thin QR factorization ensuring diagonal of R is real, positive if possible.
%
% function [Q, R] = qr_unique(A)
%
% If A is a matrix, then Q, R are matrices such that A = QR where Q'*Q = I
% and R is upper triangular. If A is real, then so are Q and R.
% This is a thin QR factorization in the sense that if A has more rows than
% columns, than Q has the same size as A.
% 
% If A has full column rank, then R has positive reals on its diagonal.
% Otherwise, it may have zeros on its diagonal.
%
% This is equivalent to a call to Matlab's qr(A, 0), up to possible
% sign/phase changes of the columns of Q and of the rows of R to ensure
% the stated properties of the diagonal of R. If A has full column rank,
% this decomposition is unique, hence the name of the function.
%
% If A is a 3D array, then Q, R are also 3D arrays and this function is
% applied to each slice separately.

resetRands(0);
X1 = make_rand_stiefel_3d_array(4,3,1);
disp("X1")
disp(X1)

resetRands(0);
U1 = stiefel_randTangentNormVector(X1);
disp("U1")
disp(U1)

Y1 = retraction_qr_stiefel(X1, U1);
disp("Y1")
disp(Y1)


end %file function

function Y = retraction_qr_stiefel(X, U, t)
    % It is necessary to call qr_unique rather than simply qr to ensure
    % this is a retraction, to avoid spurious column sign flips.
    if nargin < 3
        Y = qr_unique(X + U);
    else
        Y = qr_unique(X + t*U);
    end
end

