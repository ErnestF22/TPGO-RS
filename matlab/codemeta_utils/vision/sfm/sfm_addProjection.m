%function data=sfm_addProjectionMatrix(data)
%Requires fields calibration and pose, which are used to compute the
%projection matrices. Adds the new field projectionMatrix
function data=sfm_addProjection(data,varargin)
memberNamePose='pose';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'membernamepose'
            ivarargin=ivarargin+1;
            memberNamePose=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

data.projection=GK2P(data.(memberNamePose),data.calibration,'References');
