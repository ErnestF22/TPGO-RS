%Generate 3-D points lying on planes
%function [X,NVec,idxX]=homFlowDatasetStructure(NPlanes,varargin)
%Optional inputs
%   'NPlanes',n         number of planes (can be 1 or 2)
%   'NPointsPlane',n    number of points to generate for each plane
%Outputs
%   X       [3 x NPoints] array of 3-D coordinates
%   G       camera pose in front of the points
%   NVec    [4 x NPlanes] array of plane normals and distances
%   idxX    [1 x NPoints] vector with point-plane memberships
function [X,G,NVec,idxX]=homFlowDatasetStructure(varargin)
NPointsPlane=100;
NPlanes=1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'nplanes'
            ivarargin=ivarargin+1;
            NPlanes=varargin{ivarargin};
        case 'npointsplane'
            ivarargin=ivarargin+1;
            NPointsPlane=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NVec(:,1)=planeNdToNVec(sphere_randn(cnormalize([0;-1;1]),0.2),0);
if NPlanes>1
    NVec(:,2)=planeNdToNVec(sphere_randn(cnormalize([0;0;1]),0.2),0);
end

offsetPlanes=[zeros(3,1) [0;2;0]];

X=[];
idxX=[];
for iPlane=1:NPlanes
    X=[X planeProject(NVec(:,iPlane),rigidTransform(eye(3),offsetPlanes(:,iPlane),2*(rand(3,NPointsPlane)-0.5)))];
    idxX=[idxX iPlane*ones(1,NPointsPlane)];
end

G=RT2G(eye(3),[0;1;-3]);
