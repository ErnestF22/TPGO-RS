%Project a matrix on the set of PSD matrices
function X=projectPSD(Z)
[V,S]=eig((Z+Z')/2);
X=V*diag(max(diag(S),0))*V';
