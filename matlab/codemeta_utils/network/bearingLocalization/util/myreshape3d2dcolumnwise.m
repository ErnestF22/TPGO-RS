function C = myreshape3d2dcolumnwise(A)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
% Reshapes 3D array A 
%
% If A is n x n x k then
% C is nk x n
% C=[A(:,:,1); A(:,:,2);...;A(:,:,k)];

C = permute(A,[1 3 2]);
C = reshape(C,[],size(A,2),1);




end

