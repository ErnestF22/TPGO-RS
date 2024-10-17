%Evaluate ADMM medoids Lagrangian
%Inputs
%   nodeData    structure array as produced by admmMedoids_initNodeData

function L=admmMedoids_Lagrangian(nodeData,rho)
nbNodes=length(nodeData);
allL=nan(1,nbNodes);
for iNode=1:nbNodes
    lambdaNode=admmMedoids_centerAggregateVariable(nodeData,iNode,'lambda');
    lambdaNode=[lambdaNode(:,:,:,1) lambdaNode(:,:,:,2)];
    sz=size(lambdaNode);
    lambdaNode=squeeze(mat2cell(lambdaNode,sz(1),sz(2),ones(1,sz(3))));
    
    zNode=admmMedoids_centerAggregateVariable(nodeData,iNode,'z');
    sz=size(zNode);
    zNode=squeeze(mat2cell(zNode,sz(1),sz(2),ones(1,sz(3))));
    rhoNode=cellfun(@(x) rho*ones(1,size(x,2)),zNode,...
        'UniformOutput',false);
    
    xNode=nodeData(iNode).x;
    muNode=nodeData(iNode).mu;
    allL(iNode)=medoids_centersAssignmentCost(xNode,muNode,...
        'bias',lambdaNode,'priorCenters',zNode,rhoNode);
end
L=sum(allL);
