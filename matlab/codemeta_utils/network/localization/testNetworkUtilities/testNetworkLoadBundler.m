%function [t_node,X]=testNetworkLoadBundler(dirName, varargin)
%Load the results from bundler into a t_node structure. By default, the
%loaded poses are in the 'reference' interpretation
%Optional parameters
%   'methodAbsolutePoses',method    same as testNetworkAddGroundTruth

%%AUTORIGHTS%%

function [t_node,X,br]=testNetworkLoadBundler(dirName,varargin)
methodAbsolutePoses='reference';
thresholdCorrespCount=30;
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
            methodAbsolutePoses=varargin{ivarargin};
        case 'thresholdcount'
            ivarargin=ivarargin+1;
            thresholdCorrespCount=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

[G,X,br]=testNetworkLoadBundlerEstimate(dirName);

br.dirName=dirName;
flagHasPairs=exist(fullfile(dirName,'pairsKeysOnly'),'dir');
%prepare adjacency matrix
A=zeros(br.NCameras);
if flagHasPairs
    [pairPoses,E1,E2]=bundlerGetPairsPoses(br,... 
        'thresholdcount',thresholdCorrespCount);
    A(sub2ind(br.NCameras*[1 1],E1,E2))=1;
    idxNoPairs=getCamerasNoPairs(br);
    if length(idxNoPairs)>0
        warning('%d cameras without any edge',length(idxNoPairs));
    end
end

idxNoCameras=getCamerasInvalidPoses(G);
if length(idxNoCameras)>0
    warning('%d cameras not localized by bundler');
end

%load everything into the structure
t_node=testNetworkCreateStruct(A);
t_node=testNetworkAddGroundTruth(t_node,G,'methodabsoluteposes',methodAbsolutePoses);
if exist(fullfile(dirName,'pairsKeysOnly'),'dir')
    t_node=testNetworkAddMeasurements(t_node,'method','given',pairPoses);
end

function idx=getCamerasInvalidPoses(G)
NCameras=size(G,3);
idx=[];
for iCamera=1:NCameras
    if det(G(:,:,iCamera))==0
        idx=[idx iCamera];
    end
end

function idx=getCamerasNoPairs(br)
idx=(sum(br.correspondenceCounts)==0);

    