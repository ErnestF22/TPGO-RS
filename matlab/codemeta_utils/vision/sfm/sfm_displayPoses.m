function sfm_displayPoses(data,varargin)
memberPoses='poseEstimated';

methodAbsolutePoses='reference';
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'memberposes'
            ivarargin=ivarargin+1;
            memberPoses=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

testNetworkDisplay(data.(memberPoses),'methodAbsolutePoses',methodAbsolutePoses)

