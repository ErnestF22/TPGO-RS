%Fit essential matrix using a subset of image points and detect inliers
%function [E,residuals,flagInlier]=epipolarEstimateESampleAndOutliers(x1,x2,sampleIdx,varargin)
%Inputs
%   x1,x2       [2 x Nx] arrays with the image points
%   sampleIdx   [1 x Nx] boolean array or [1 x NSample] index array
%               indicating which image points should be used for estimation
%Optional arguments
%   'threshold',t           Threshold for computing the flagInlier output
%   'methodResiduals', m    Name of the method used for computing the
%                           residuals (see epipolarResiduals)
%   'methodE'               Method used for estimation
%       '5pt'                   5-point algorithm
%       '5pt+validation'        5-point algorithm + use the remaining
%            points in the sample and the threshold to validate each
%            solution
%       '8pt'                   linear 8-point algorithm (actually uses all
%             the points provided for the solution)
%If 'methodE' is not provided or empty, the method is chosen depending on
%the number of points. If the number of points in flagSample is 5,6, or 7,
%use the 5-point algorithm plus validation. Otherwise, use the 8-point
%algorithm.
%Outputs
%   E           estimated essential matrix(ces). Beware that different methods
%       might return different numbers of solutions, varying from none
%       (5pt+validation) to one ('8pt') to many (5pt and 5pt+validation).
%   residuals   computed residuals for each image pair
%   flagInlier  boolean vectors indicating if each image pair is an inlier
function [E,residuals,flagInlier]=epipolarEstimateESampleAndOutliers(x1,x2,sampleIdx,varargin)
threshold=1e-2;
methodResiduals='sampsonabs';
methodE=[];

%detect if flagSample contains booleans or indeces and count number of
%points to use for estimation
if isempty(sampleIdx)
    sampleIdx=true(1,size(x1,2));
end

if length(unique([0 1 sampleIdx]))==2
    NSample=sum(sampleIdx);
    sampleIdx=find(sampleIdx);
else
    NSample=length(sampleIdx);
end

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'threshold'
            ivarargin=ivarargin+1;
            threshold=lower(varargin{ivarargin});
        case 'methodresiduals'
            ivarargin=ivarargin+1;
            methodResiduals=lower(varargin{ivarargin});
        case 'methode'
            ivarargin=ivarargin+1;
            methodE=lower(varargin{ivarargin});
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if NSample<5
    error('Number of points in sampleIdx should be at least 5')
end

if isempty(methodE)
    if NSample==5
        methodE='5pt';
    elseif NSample>5 && NSample<8
        methodE='5pt+validation';
    elseif NSample>=8
        methodE='8pt';
    end
end

switch lower(methodE)
    case '5pt'
        E=epipolarNormalize(solveE(homogeneous(x1(:,sampleIdx),3), homogeneous(x2(:,sampleIdx),3)));
    case '5pt+validation'
        sampleIdxEst=sampleIdx(1:5);
        sampleIdxVal=sampleIdx(6:end);
        E=epipolarNormalize(solveE(homogeneous(x1(:,sampleIdxEst),3), homogeneous(x2(:,sampleIdxEst),3)));
        if ~isempty(sampleIdxVal)
            residuals=epipolarResiduals(E,x1(:,sampleIdxVal),x2(:,sampleIdxVal),methodResiduals);
            flagValid=all(residuals<threshold,2);
            E=E(:,:,flagValid);
        end
    case '8pt'
        E=epipolarEstimateE8pt(x1(:,sampleIdx),x2(:,sampleIdx));
    otherwise
        error('MethodE invalid')
end

residuals = epipolarResiduals(E,x1,x2,methodResiduals);
flagInlier = residuals < threshold;
