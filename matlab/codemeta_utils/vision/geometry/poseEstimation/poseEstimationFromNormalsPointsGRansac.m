function [G,output]=poseEstimationFromNormalsPointsGRansac(n1,d1,x22D,threshold,varargin)
optsRansac={};
methodRefine='none';
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'optsransac'
            ivarargin=ivarargin+1;
            optsRansac=varargin{ivarargin};
        case 'refineall'
            methodRefine='all';
        case 'refineend'
            methodRefine='end';
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

switch methodRefine
    case {'none','end'}
        flagRefine=false;
    case 'all'
        flagRefine=true;
    otherwise
        error('methodRefine not valid')
end

data=[n1;d1;x22D];
[G,output]=ransac(data,@(subData) funModelEstimation(subData,flagRefine),9,@funResiduals,threshold,optsRansac{:});
if strcmp(methodRefine,'end')
    G=funModelEstimation(data(:,output.flagInliers),true);
end

function G=funModelEstimation(subData,flagRefine)
n1=subData(1:3,:);
d1=subData(4,:);
x22D=subData(5:6,:);
G=poseEstimationFromNormalsPointsG(n1,d1,x22D);
if flagRefine
    G=poseEstimationFromNormalsPointsRefineG(n1,d1,x22D,G);
end

function e=funResiduals(G,data)
n1=data(1:3,:);
d1=data(4,:);
x22D=data(5:6,:);
e=poseEstimationFromNormalsPointsResidualsG(G,n1,d1,x22D);
