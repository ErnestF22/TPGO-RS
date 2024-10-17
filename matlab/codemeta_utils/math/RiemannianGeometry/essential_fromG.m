%Compute a point in the QREM corresponding to a pair of rigid poses
%function Q=essential_fromG(G1,G2,varargin)
%See essential_fromRT for options
function Q=essential_fromG(G1,G2,varargin)
Q=essential_fromRT(G2R(G1),G2T(G1),G2R(G2),G2T(G2),varargin{:});
