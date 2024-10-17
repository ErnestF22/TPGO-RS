%Add pose of each camera to the data
%function data=sfm_addPose(data,G,varargin)
%Adds field poseTruth to data structure. If G is omitted, initialize all to the
%identity. The poses are in the 'reference' interpretation.
%Optional inputs
%   memberName  name of the field to which the pose should be added
function data=sfm_poseAdd(data,G,varargin)
flagGivenG=exist('G','var') && ~isempty(G);
memberName='poseTruth';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'membername'
            ivarargin=ivarargin+1;
            memberName=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ~flagGivenG
    NImg=size(data.imgSize,2);
    G=eye(4,4,NImg);
end
data.(memberName)=G;
