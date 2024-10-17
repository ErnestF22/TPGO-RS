%Prepare network for averaging poses
%function t_node=multiCalibration_networkPrepare(poses,c,flagUseCovariances)
%Convert relative pose estimates and correspondences into a network data
%structure to be used in localization_MLE_rigid
%The function proceeds as follows:
%   1)  Get list of node names from the correspondence data using
%       multiCalibration_correspondencesToNodeNames
%   2)  Create list of edges from the 'nodeNames' field in the
%       correspondence data
%   3)  Initialize the newtork data structure
%   4)  Fill poses for the edge list using the arguments poses and c
%   5)  Add these poses to the network data structure
function [t_node,nodeNames]=multiCalibration_networkPrepare(poses,c,flagUseCovariances)
nodeNames=multiCalibration_correspondencesToNodeNames(c);
NNodes=length(nodeNames);

%create edge list from correspondences
NEdges=length(c);
E=zeros(NEdges,2);
for iEdge=1:NEdges
    for iEndPoint=1:2
        [~,E(iEdge,iEndPoint)]=ismember(c{iEdge}.nodeNames{iEndPoint},nodeNames);
    end
end
E=[E;fliplr(E)];

%initialize network data structure
t_node=testNetworkCreateStruct(E,'Edges',NNodes,'Type','Single');
E=t_node.E; 
%Note: the edges in the structure do not align with those in the correspondences

%fill poses using correspondences and edge list
NEdges=size(E,1);
gij=zeros(4,4,NEdges);
Gammaij=zeros(6,6,NEdges);
for iCorrespondence=1:length(c)
    %get node index from name
    [~,iNode]=ismember(c{iCorrespondence}.nodeNames{1},nodeNames);
    [~,jNode]=ismember(c{iCorrespondence}.nodeNames{2},nodeNames);
    %get edge index from node indeces and fill in the pose
    iEdge=find(E(:,1)==iNode & E(:,2)==jNode);
    gij(:,:,iEdge)=poses(iCorrespondence).G;
    Gammaij(:,:,iEdge)=poses(iCorrespondence).Sigma;
    %get reverse edge indeces and fill in the pose
    iEdgeReverse=find(E(:,1)==jNode & E(:,2)==iNode);
    gij(:,:,iEdgeReverse)=invg(poses(iCorrespondence).G);
    Gammaij(:,:,iEdgeReverse)=poses(iCorrespondence).Sigma;
end

%add relative measurements to the network data structure
t_node=testNetworkAddMeasurements(t_node,'Method','Given',gij);
if ~flagUseCovariances
    t_node=testNetworkAddDispersionMatricesRT(t_node,'identity');
else
    t_node=testNetworkAddDispersionMatricesRT(t_node,'given',Gammaij);
end

