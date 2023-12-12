function [discs, maxEigenvalue] = gershgorinDisc(A)
% Return the Gershgorin Discs of the complex square matrix A
% NOTE: radii are sums of rows
% INPUTS:
%   A := a [nxn] complex matrix
% OUTPUTS:
%   discs := [nx2] vector of [center, radius]
%   maxEigenvalue := largest real eigenvalue encompassed by a Gershgorin
%       Disc

B = abs(A); % Take the absolute value of all entries

for i = 1:length(A)
    discs(i,:) = [A(i,i), sum(B(i,:))-B(i,i)];
end

maxEigenvalue = max(sum(discs,2));
end