%Find a low-rank SDP approximation for a matrix of pair-wise rotations
%function GEst=sfm_rawRotationsLowRankSDP(G,varargin)
%Optional Inputs
%   'spectrahedral'     Add constraints given by convex hull of SO(3)
function GEst=sfm_rawRotationsLowRankSDP(G,varargin)
NRotations=size(G,1)/3;
idxRotation=reshape(1:3*NRotations,[3 NRotations]);
flagSpectrahedral = false;

ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'spectrahedral'
            flagSpectrahedral = true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if flagSpectrahedral
    B=generateSpectrahedralConstraints();
end

%solve optimization
cvx_begin SDP quiet
    variable GEst(3*NRotations,3*NRotations) symmetric 
    maximize trace(G'*GEst)
    subject to
        GEst >= 0;
        for iRotation=1:NRotations
            iidx=idxRotation(:,iRotation);
            GEst(iidx,iidx) == eye(3);
            if flagSpectrahedral
                for jRotation=iRotation+1:NRotations
                    jidx=idxRotation(:,jRotation);
                    reshape(B*vec(GEst(iidx,jidx)),4,4) + eye(4) >= 0;
                end
            end
        end
cvx_end
GEst=full(GEst);


function B=generateSpectrahedralConstraints()
A{1} = diag([1 1 -1 -1]);
A{2} = zeros(4,4); A{2}(1,4) = -1; A{2}(2,3) = 1; A{2} = A{2} + A{2}';
A{3} = zeros(4,4); A{3}(1,3) = 1; A{3}(2,4) = 1; A{3} = A{3} + A{3}';
A{4} = zeros(4,4); A{4}(1,4) = 1; A{4}(2,3) = 1; A{4} = A{4} + A{4}';
A{5} = diag([1 -1 1 -1]);
A{6} = zeros(4,4); A{6}(1,2) = -1; A{6}(3,4) = 1; A{6} = A{6} + A{6}';
A{7} = zeros(4,4); A{7}(1,3) = -1; A{7}(2,4) = 1; A{7} = A{7} + A{7}';
A{8} = zeros(4,4); A{8}(1,2) = 1; A{8}(3,4) = 1; A{8} = A{8} + A{8}';
A{9} = diag([1 -1 -1 1]);
B = zeros(16,9);
for i = 1:9
    B(:,i) = A{i}(:);
end
