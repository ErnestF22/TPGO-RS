function A = myreshape3d2dcolumnwise_inv(C)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%
% Reshapes 3D array A 
%
% If A is n x n x k then
% C is nk x n
% C=[A(:,:,1); A(:,:,2);...;A(:,:,k)];

k = size(C,1)/size(C,2);
A = zeros( size(C,2),size(C,2),k);
for i=1:k
   A(:,:,i) = C( size(C,2)*(i-1)+1:size(C,2)*i ,:); 
end
 




end
