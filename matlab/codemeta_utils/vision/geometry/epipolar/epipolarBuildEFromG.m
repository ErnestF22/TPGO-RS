%Construct essential matrix from poses
%function E=epipolarBuildEFromG(G1,G2,varargin)
%See epipolarBuildEFromRT
function E=epipolarBuildEFromG(G1,G2,varargin)
[R1,T1]=G2RT(G1);
if exist('G2','var') && ~isempty(G2)
    [R2,T2]=G2RT(G2);
else
    R2=[];
    T2=[];
end
E=epipolarBuildEFromRT(R1,T1,R2,T2,varargin{:});
