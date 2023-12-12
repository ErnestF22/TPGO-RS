%Generate synthetic datasets for testing QuickShift and QuickMatch
%function [X,membershipPrior,NPoints]=quickshift_test_datasets(datasetName,varargin)
%Input arguments
%   datasetName     Use 'clustering' for a clustering problem, and
%                   'matching' for a matching problem
%Optional arguments
%   'NPointsClass'    Number of points in each class
%   
function [X,membershipPrior,NPoints]=quickshift_test_datasets(datasetName,varargin)
if ~exist('datasetName','var') || isempty(datasetName)
    datasetName='clustering';
end
NPointsCluster=50;
NOutliersCluster=0;
sigmaCorrespondences=0.01;
NClusters=4;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'npointsclass'
            ivarargin=ivarargin+1;
            NPointsCluster=varargin{ivarargin};
        case 'nclasses'
            ivarargin=ivarargin+1;
            NClusters=varargin{ivarargin};
        case 'noutliersclass'
            ivarargin=ivarargin+1;
            NOutliersCluster=varargin{ivarargin};
        case 'sigmacorrespondences'
            ivarargin=ivarargin+1;
            sigmaCorrespondences=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

switch lower(datasetName)
    case 'clustering'
        if NClusters~=4
            error('This function needs to be extended to allow other numbers of Gaussian clusters')
        end
        %generate points in 4 clusters
        X=[randn(2,NPointsCluster) randn(2,NPointsCluster)+4 ...
            [randn(1,NPointsCluster)+4;0.5*randn(1,NPointsCluster)] ...
            [0.5*randn(1,NPointsCluster);randn(1,NPointsCluster)+4]];
        membershipPrior=reshape(repmat(1:4,NPointsCluster,1),1,[]);
        
    case 'matching'
        X0=rand(2,NPointsCluster);
        X=[];
        for iCluster=1:NClusters
            X=[X X0+sigmaCorrespondences*randn(size(X0))];
        end
        membershipPrior=reshape(repmat(1:NClusters,NPointsCluster,1),1,[]);
end        

if NOutliersCluster>0
    minX=min(X,[],2);
    maxX=max(X,[],2);
    widthX=maxX-minX;
    XOutliers=(1.4*widthX*ones(1,NClusters*NOutliersCluster)).*rand(2,NClusters*NOutliersCluster)+(minX-0.2*widthX)*ones(1,NClusters*NOutliersCluster);
    membershipPriorOutliers=reshape(repmat(1:4,NOutliersCluster,1),1,[]);
    
    X=[X XOutliers];
    membershipPrior=[membershipPrior membershipPriorOutliers];
end

NPoints=NClusters*NPointsCluster;
