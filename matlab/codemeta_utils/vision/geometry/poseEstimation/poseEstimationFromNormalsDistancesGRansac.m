function [G,output]=poseEstimationFromNormalsDistancesGRansac(n1,d1,n2,d2,threshold,varargin)
optsRansac={};
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'optsransac'
            ivarargin=ivarargin+1;
            optsRansac=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

data=[n1;d1;n2;d2];
[G,output]=ransac(data,@funModelEstimation,3,@funResiduals,threshold,optsRansac{:});

function G=funModelEstimation(subData)
n1=subData(1:3,:);
d1=subData(4,:);
n2=subData(5:7,:);
d2=subData(8,:);
G=poseEstimationFromNormalsDistancesG(n1,d1,n2,d2);

function e=funResiduals(G,data)
n1=data(1:3,:);
d1=data(4,:);
n2=data(5:7,:);
d2=data(8,:);
[eRot,eTransl]=poseEstimationFromNormalsDistancesResidualsG(G,n1,d1,n2,d2);
e=sqrt(eRot.^2+eTransl.^2);
