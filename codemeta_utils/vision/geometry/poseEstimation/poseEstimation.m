function GEst=poseEstimation(X,x,varargin)
methodAbsolutePoses='pose';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'references'
            methodAbsolutePoses='reference';
        case 'poses'
            methodAbsolutePoses='pose';
        case 'methodabsoluteposes'
            ivarargin=ivarargin+1;
            methodAbsolutePoses=lower(varargin{ivarargin});
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end


GEstQuad=poseEstimationQuadratic(X,x,'methodAbsolutePoses',methodAbsolutePoses);
GEst=poseEstimationRefineFromG(GEstQuad,X,x(1:2,:),'methodAbsolutePoses',methodAbsolutePoses);
