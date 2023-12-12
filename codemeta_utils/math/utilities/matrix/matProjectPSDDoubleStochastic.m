%Project a matrix on the intersection of the sets of PSD and double stochastic matrices
%function X=projectPSDDoubleStochastic(Z)
function X=matProjectPSDDoubleStochastic(Z,varargin)
X=matProjectPSD(matProjectDoubleStochastic(Z,varargin{:}));