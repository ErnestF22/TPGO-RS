function alpha=fieldEstimationLocal(fieldData,idxList,kList,funEstimate,varargin)
dAlpha=8;            %dimension of the parameter vectors (is a constant)

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'dimalpha'
            ivarargin=ivarargin+1;
            dAlpha=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NX=size(fieldData,2);
alpha=zeros(dAlpha,NX);

%initialize alpha from local estimates
for iX=1:NX
    alpha(:,iX)=funEstimate(fieldData(:,idxList(iX,1:kList(iX))));
end
